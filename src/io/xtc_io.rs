// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of functions for reading and writing xtc files.

use std::marker::PhantomData;
use std::os::raw::{c_float, c_int};
use std::path::Path;

use crate::errors::{ReadXdrError, WriteXdrError, XdrError};
use crate::io::xdrfile::{
    self, OpenMode, XdrFile, XdrFrameData, XdrGroupWriter, XdrRangeReader, XdrReader, XdrWriter,
};
use crate::iterators::AtomIterator;
use crate::structures::{group::Group, vector3d::Vector3D};
use crate::system::general::System;

/**************************/
/*       READING XTC      */
/**************************/

/// Iterator over an xtc file.
pub struct XtcReader<'a> {
    system: *mut System,
    xtc: XdrFile,
    phantom: PhantomData<&'a mut System>,
}

impl<'a> XdrReader<'a> for XtcReader<'a> {
    /// Create an iterator over an xtc file.
    fn new(system: &'a mut System, filename: impl AsRef<Path>) -> Result<Self, ReadXdrError> {
        let n_atoms = system.get_n_atoms();

        // sanity check the number of atoms
        match XdrFile::check_xtc(filename.as_ref(), n_atoms) {
            Err(e) => return Err(e),
            Ok(false) => {
                return Err(ReadXdrError::AtomsNumberMismatch(Box::from(
                    filename.as_ref(),
                )))
            }
            Ok(true) => (),
        };

        // open the xtc file and save the handle to it
        let xtc = match XdrFile::open_xdr(filename.as_ref(), OpenMode::Read) {
            Ok(x) => x,
            Err(XdrError::FileNotFound(x)) => return Err(ReadXdrError::FileNotFound(x)),
            Err(XdrError::InvalidPath(x)) => return Err(ReadXdrError::InvalidPath(x)),
        };

        Ok(XtcReader {
            system,
            xtc,
            phantom: PhantomData,
        })
    }

    /// Transform `XtcReader` iterator into `XdrRangeReader<XtcReader>` iterator with the specified time range.
    /// Time (`start_time` and `end_time`) must be provided in picoseconds.
    ///
    /// ## Details
    /// The frames with time lower than `start_time` are skipped over, i.e., the coordinates of the atoms
    /// are not read at all, making the `XdrRangeReader` very efficient. Note however that calling this function nonetheless
    /// involves some initial overhead, as it needs to locate the starting point of the iteration.
    ///
    /// Iteration is ended at the frame with time corresponding to `end_time` or once the end of the xtc file is reached.
    /// The range is inclusive on both ends, i.e., frames with `time = start_time` and `time = end_time` will be included in the iteration.
    ///
    /// ## Returns
    /// `XdrRangeReader<XtcReader>` if the the specified time range is valid. Else returns `ReadXdrError`.
    ///
    /// ## Example
    /// Creating an `XdrRangeReader` iterator over an xtc file.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// // load system from file
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// // create the iterator for partial reading of an xtc file
    /// let range_iterator = system
    ///     .xtc_iter("trajectory.xtc")
    ///     .unwrap()
    ///     .with_range(10_000.0, 100_000.0)
    ///     .unwrap();
    /// ```
    ///
    /// ## Notes
    /// - If the frame corresponding to the `start_time` doesn't exist in the xtc file,
    /// the iterator starts at the frame closest in time to but greater than the `start_time`.
    /// - If either the `start_time` or the `end_time` is negative, it results in a `ReadXdrError::TimeRangeNegative` error.
    /// - If the `start_time` is greater than the `end_time`, it results in a `ReadXdrError::InvalidTimeRange` error.
    /// - If the `start_time` exceeds the time of any frame in the xtc file, it results in a `ReadXdrError::StartNotFound` error.
    fn with_range(
        self,
        start_time: f32,
        end_time: f32,
    ) -> Result<XdrRangeReader<'a, XtcReader<'a>>, ReadXdrError> {
        xdrfile::sanity_check_timerange(start_time, end_time)?;

        let reader = XdrRangeReader::new(self, start_time, end_time);

        // jump to the start of iteration
        unsafe {
            if xdrfile::xtc_jump_to_start(reader.xdrreader.xtc.handle, reader.start_time as c_float)
                != 0 as c_int
            {
                return Err(ReadXdrError::StartNotFound(start_time.to_string()));
            }
        }

        Ok(reader)
    }
}

impl<'a> Iterator for XtcReader<'a> {
    type Item = Result<&'a mut System, ReadXdrError>;

    /// Read next frame in an xtc file.
    ///
    /// ## Returns
    /// `None` in case the file has been fully read.
    /// `Some(&mut System)` in case the frame has been read successfully.
    /// `Some(ReadXdrError)` in case an error occured while reading.
    fn next(&mut self) -> Option<Self::Item> {
        unsafe {
            match XtcFrameData::from_frame(&mut self.xtc, (*self.system).get_n_atoms()) {
                None => None,
                Some(Err(e)) => Some(Err(e)),
                Some(Ok(data)) => {
                    data.update_system(&mut *self.system);
                    Some(Ok(&mut *self.system))
                }
            }
        }
    }
}

impl<'a> Iterator for XdrRangeReader<'a, XtcReader<'a>> {
    type Item = Result<&'a mut System, ReadXdrError>;

    /// Read the next frame in the specified time range of the xtc file.
    fn next(&mut self) -> Option<Self::Item> {
        unsafe {
            let system = self.xdrreader.system;

            match XtcFrameData::from_frame(&mut self.xdrreader.xtc, (*system).get_n_atoms()) {
                None => None,
                Some(Err(e)) => Some(Err(e)),
                Some(Ok(data)) => {
                    if (*system).get_simulation_time() >= self.end_time {
                        None
                    } else {
                        data.update_system(&mut *system);
                        Some(Ok(&mut *system))
                    }
                }
            }
        }
    }
}

/// ## Methods for reading xtc files.
impl System {
    /// Create an `XtcReader` structure which is an iterator over an xtc file.
    ///
    /// ## Returns
    /// `XtcReader` if the xtc file exists and matches the structure file.
    /// Else returns `ReadXdrError`.
    ///
    /// ## Examples
    /// Iterating through an xtc trajectory and calculating
    /// and printing the current center of geometry of the system.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// // load system from file
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// // create the `XtcReader`
    /// let mut trajectory = match system.xtc_iter("trajectory.xtc") {
    ///     Ok(x) => x,
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;
    ///     }
    /// };
    ///
    /// // iterate through the trajectory
    /// // `raw_frame` is either mutable reference to the `System` in the current state
    /// // or `ReadXtcError` in case reading of the frame has failed
    /// for raw_frame in trajectory {
    ///     match raw_frame {
    ///         // calculate center of the system in the current frame and print it
    ///         Ok(frame) => println!("{:?}", frame.group_get_center("all")),
    ///         Err(e) => {
    ///             eprintln!("{}", e);
    ///             return;
    ///         }
    ///     }
    /// }
    /// ```
    /// Much more concise way using the `?` operator.
    /// ```no_run
    /// use groan_rs::prelude::*;
    /// use groan_rs::errors::ReadXdrError;
    ///
    /// fn example_fn() -> Result<(), ReadXdrError> {
    ///     // load system from file
    ///     let mut system = System::from_file("system.gro").unwrap();
    ///
    ///     // loop through the trajectory
    ///     for raw_frame in system.xtc_iter("trajectory.xtc")? {
    ///         let frame = raw_frame?;
    ///         println!("{:?}", frame.group_get_center("all"));
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    ///
    /// You can also iterate over just a part of the trajectory.
    /// Here, only frames in the time range 10-100 ns will be read.
    /// **Do not use the `skip` method for skipping over the initial frames
    /// of the trajectory as that is very inefficient.**
    /// ```no_run
    /// use groan_rs::prelude::*;
    /// use groan_rs::errors::ReadXdrError;
    ///
    /// fn example_fn() -> Result<(), ReadXdrError> {
    ///     let mut system = System::from_file("system.gro").unwrap();
    ///
    ///     for raw_frame in system
    ///         .xtc_iter("trajectory.xtc")?
    ///         .with_range(10_000.0, 100_000.0)?
    ///     {
    ///         let frame = raw_frame?;
    ///         println!("{:?}", frame.group_get_center("all"));
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    ///
    /// ## Warning
    /// - Only orthogonal simulation boxes are currently supported!
    ///
    /// ## Notes
    /// - The function checks whether the number of atoms in the system corresponds to the number of atoms in the xtc file.
    /// - The `System` structure is modified while iterating through the xtc file.
    pub fn xtc_iter(&mut self, filename: impl AsRef<Path>) -> Result<XtcReader, ReadXdrError> {
        XtcReader::new(self, filename)
    }
}

/**************************/
/*       WRITING XTC      */
/**************************/

/// Structure for writing xtc files.
/// Each `XtcWriter` instance is tightly coupled with a corresponding `System` structure.
/// If you make updates to the `System` structure, such as during iteration with `System::xtc_iter()`,
/// and subsequently write an XTC frame using `XtcWriter::write_frame()`, the modifications
/// made to the `System` will be reflected in the written frame.
///
/// `XtcWriter` implements the `XdrWriter` trait.
pub struct XtcWriter {
    system: *const System,
    xtc: XdrFile,
}

impl XdrWriter for XtcWriter {
    /// Open a new xtc file for writing.
    ///
    /// ## Returns
    /// An instance of `XtcWriter` structure or `WriteXdrError` in case the file can't be created.
    ///
    /// ## Example
    /// Create a new xtc file for writing and associate a system with it.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let system = System::from_file("system.gro").unwrap();
    ///
    /// let mut writer = XtcWriter::new(&system, "output.xtc").unwrap();
    /// ```
    fn new(system: &System, filename: impl AsRef<Path>) -> Result<XtcWriter, WriteXdrError> {
        // create the xtc file and save a handle to it
        let xtc = match XdrFile::open_xdr(filename.as_ref(), OpenMode::Write) {
            Ok(x) => x,
            Err(XdrError::FileNotFound(x)) => return Err(WriteXdrError::CouldNotCreate(x)),
            Err(XdrError::InvalidPath(x)) => return Err(WriteXdrError::InvalidPath(x)),
        };

        Ok(XtcWriter { system, xtc })
    }

    /// Write the current state of the system into an open xtc file.
    ///
    /// ## Returns
    /// - `Ok` if the frame has been successfully written. Otherwise `WriteXdrError`.
    ///
    /// ## Example
    /// Reading and writing an xtc file.
    /// ```no_run
    /// use groan_rs::prelude::*;
    /// use std::error::Error;
    ///
    /// fn example_fn() -> Result<(), Box<dyn Error>> {
    ///     // load system from file
    ///     let mut system = System::from_file("system.gro")?;
    ///
    ///     // create an xtc file for writing and associate it with the system
    ///     let mut writer = XtcWriter::new(&system, "output.xtc")?;
    ///
    ///     // loop through the trajectory
    ///     for raw_frame in system.xtc_iter("trajectory.xtc")? {
    ///         // check for errors
    ///         let _ = raw_frame?;
    ///
    ///         // write the current frame into `output.xtc`
    ///         writer.write_frame()?;
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    ///
    /// ## Notes
    /// - Precision for writing the xtc file is taken from the `System` structure.
    fn write_frame(&mut self) -> Result<(), WriteXdrError> {
        unsafe {
            let n_atoms = (*self.system).get_n_atoms();

            // prepare coordinate matrix
            let mut coordinates = vec![[0.0, 0.0, 0.0]; n_atoms];
            for (i, atom) in (*self.system).atoms_iter().enumerate() {
                let pos = atom.get_position();
                coordinates[i] = [pos.x, pos.y, pos.z];
            }

            // write the xtc frame
            let return_code = xdrfile::write_xtc(
                self.xtc.handle,
                n_atoms as c_int,
                (*self.system).get_simulation_step() as i32,
                (*self.system).get_simulation_time(),
                &mut xdrfile::simbox2matrix((*self.system).get_box_as_ref()),
                coordinates.as_mut_ptr(),
                (*self.system).get_precision() as f32,
            );

            if return_code != 0 {
                return Err(WriteXdrError::CouldNotWrite);
            }
        }

        Ok(())
    }
}

/// Structure for writing groups of atoms into xtc files.
/// Each `XtcGroupWriter` is tightly coupled with a corresponding `Group` from `System` structure.
/// If you make updates to the `System` structure, such as during the iteration with `System::xtc_iter()`,
/// and subsequently write an XTC frame using `XtcGroupWriter::write_frame()`, the modifications
/// made to the `System` will be reflected in the written frame.
///
/// Note that the purpose of the `XtcGroupWriter` is writing valid xtc files with consistent number of atoms.
/// Therefore, the `XtcGroupWriter` always works with the original provided group of atoms.
/// If you change the meaning of this group after constructing `XtcGroupWriter` (by overwriting the group),
/// `XtcGroupWriter` will still use the original group of atoms.
/// If you completely remove the original group, `XtcGroupWriter` will still maintain a working copy of it.
///
/// `XtcGroupWriter` implements the `XdrGroupWriter` trait.
pub struct XtcGroupWriter {
    system: *const System,
    xtc: XdrFile,
    /// This is a deep copy of the group from the system. `XtcGroupWriter` must always work, even if the user removes the group from the `System` or overwrites it.
    group: Group,
}

impl XdrGroupWriter for XtcGroupWriter {
    /// Open a new xtc file for writing and associate a specific group from a specific system with it.
    ///
    /// ## Returns
    /// An instance of `XtcGroupWriter` structure or `WriteXdrError` in case the file can't be created
    /// or the group does not exist.
    ///
    /// ## Example
    /// Create a new xtc file for writing and associate a group with it.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("My Group", "resid 1-4").unwrap();
    ///
    /// let mut writer = XtcGroupWriter::new(&system, "My Group", "output.xtc").unwrap();
    /// ```
    fn new(
        system: &System,
        group_name: &str,
        filename: impl AsRef<Path>,
    ) -> Result<XtcGroupWriter, WriteXdrError> {
        // get copy of the group
        let group = match system.get_groups_as_ref().get(group_name) {
            None => return Err(WriteXdrError::GroupNotFound(group_name.to_owned())),
            Some(g) => g.clone(),
        };

        // create the xtc file and save a handle to it
        let xtc = match XdrFile::open_xdr(filename.as_ref(), OpenMode::Write) {
            Ok(x) => x,
            Err(XdrError::FileNotFound(x)) => return Err(WriteXdrError::CouldNotCreate(x)),
            Err(XdrError::InvalidPath(x)) => return Err(WriteXdrError::InvalidPath(x)),
        };

        Ok(XtcGroupWriter { system, xtc, group })
    }

    /// Write the current state of the group into an open xtc file.
    ///
    /// ## Returns
    /// - `Ok` if the frame has been successfully written. Otherwise `WriteXdrError`.
    ///
    /// ## Example
    /// Reading an xtc file and writing only the atoms corresponding to an ndx group `Protein` into the output xtc file.
    /// ```no_run
    /// use groan_rs::prelude::*;
    /// use std::error::Error;
    ///
    /// fn example_fn() -> Result<(), Box<dyn Error>> {
    ///     // load system from file
    ///     let mut system = System::from_file("system.gro")?;
    ///     system.read_ndx("index.ndx")?;
    ///
    ///     // create an xtc file for writing and associate it with the group `Protein`
    ///     let mut writer = XtcGroupWriter::new(&system, "Protein", "output.xtc")?;
    ///
    ///     // loop through the trajectory
    ///     for raw_frame in system.xtc_iter("trajectory.xtc")? {
    ///         // check for errors
    ///         let _ = raw_frame?;
    ///
    ///         // write the current state of the group `Protein` into `output.xtc`
    ///         writer.write_frame()?;
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    ///
    /// ## Notes
    /// - Precision for writing the xtc file is taken from the `System` structure.
    fn write_frame(&mut self) -> Result<(), WriteXdrError> {
        unsafe {
            let n_atoms = self.group.get_n_atoms();

            // create an iterator over the atoms of the group
            let iterator = AtomIterator::new(
                (*self.system).get_atoms_as_ref(),
                self.group.get_atom_ranges(),
                (*self.system).get_box_as_ref(),
            );

            // prepare coordinate matrix
            let mut coordinates = vec![[0.0, 0.0, 0.0]; n_atoms];
            for (i, atom) in iterator.enumerate() {
                let pos = atom.get_position();
                coordinates[i] = [pos.x, pos.y, pos.z];
            }

            // write the xtc frame
            let return_code = xdrfile::write_xtc(
                self.xtc.handle,
                n_atoms as c_int,
                (*self.system).get_simulation_step() as i32,
                (*self.system).get_simulation_time(),
                &mut xdrfile::simbox2matrix((*self.system).get_box_as_ref()),
                coordinates.as_mut_ptr(),
                (*self.system).get_precision() as f32,
            );

            if return_code != 0 {
                return Err(WriteXdrError::CouldNotWrite);
            }
        }

        Ok(())
    }
}

/******************************/
/*     PRIVATE FUNCTIONS      */
/******************************/

impl XdrFile {
    /// Check that the number of atoms in an unopened xtc file matches the expected number.
    fn check_xtc(filename: impl AsRef<Path>, n_atoms: usize) -> Result<bool, ReadXdrError> {
        unsafe {
            let c_path = match xdrfile::path2cstring(filename.as_ref()) {
                Ok(x) => x,
                Err(_) => return Err(ReadXdrError::InvalidPath(Box::from(filename.as_ref()))),
            };

            let mut xtc_atoms: c_int = 0;

            if xdrfile::read_xtc_natoms(c_path.as_ptr(), &mut xtc_atoms) != 0 {
                // reading the file failed
                return Err(ReadXdrError::FileNotFound(Box::from(filename.as_ref())));
            }

            // if reading was successful
            if xtc_atoms == n_atoms as c_int {
                Ok(true)
            } else {
                Ok(false)
            }
        }
    }
}

/// Structure containing data read from an xtc file.
#[derive(Debug)]
struct XtcFrameData {
    step: c_int,
    time: c_float,
    boxvector: [[c_float; 3usize]; 3usize],
    precision: c_float,
    coordinates: Vec<[f32; 3]>,
}

impl XdrFrameData for XtcFrameData {
    /// Read `XtcFrameData` from an xtc frame.
    ///
    /// ## Returns
    /// - `None` if the file has been read completely
    /// - `Some(XtcFrameData)` if the frame has been successfully read.
    /// - `Some(ReadXdrError)` if the frame could not be read.
    fn from_frame(xdrfile: &mut XdrFile, n_atoms: usize) -> Option<Result<Self, ReadXdrError>> {
        let mut step: c_int = 0;
        let mut time: c_float = 0.0;
        let mut boxvector: [[c_float; 3usize]; 3usize] = [[0.0; 3]; 3];
        let mut precision: c_float = 0.0;
        let mut coordinates = vec![[0.0, 0.0, 0.0]; n_atoms];

        unsafe {
            let return_code = xdrfile::read_xtc(
                xdrfile.handle,
                n_atoms as c_int,
                &mut step,
                &mut time,
                &mut boxvector,
                coordinates.as_mut_ptr(),
                &mut precision,
            );

            match return_code {
                // reading successful and there is more to read
                0 => Some(Ok(XtcFrameData {
                    step,
                    time,
                    boxvector,
                    precision,
                    coordinates,
                })),
                // file is completely read
                11 => None,
                // error occured
                _ => Some(Err(ReadXdrError::FrameNotFound)),
            }
        }
    }

    /// Update the `System` structure based on data from `XtcFrameData`.
    fn update_system(self, system: &mut System) {
        unsafe {
            for (i, atom) in system.get_atoms_as_ref_mut().iter_mut().enumerate() {
                atom.set_position(&Vector3D::from(*self.coordinates.get_unchecked(i)));
            }

            // update the system
            system.set_simulation_step(self.step as u64);
            system.set_simulation_time(self.time);
            system.set_box(xdrfile::matrix2simbox(self.boxvector));
            system.set_precision(self.precision as u64);
        }
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::assert_approx_eq;
    use std::fs::File;
    use tempfile::NamedTempFile;

    #[test]
    fn read_xtc() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        // first frame
        system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .next();

        assert_eq!(system.get_simulation_step(), 0);
        assert_eq!(system.get_precision(), 100);
        assert_approx_eq!(f32, system.get_simulation_time(), 0.0);
        assert_approx_eq!(f32, system.get_box_as_ref().x, 13.034535);
        assert_approx_eq!(f32, system.get_box_as_ref().y, 13.034535);
        assert_approx_eq!(f32, system.get_box_as_ref().z, 11.228164);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_approx_eq!(f32, atom1.get_position().x, 5.97);
        assert_approx_eq!(f32, atom1.get_position().y, 7.03);
        assert_approx_eq!(f32, atom1.get_position().z, 7.69);

        assert_approx_eq!(f32, atom1.get_velocity().x, -0.0683);
        assert_approx_eq!(f32, atom1.get_velocity().y, 0.1133);
        assert_approx_eq!(f32, atom1.get_velocity().z, 0.0005);

        assert_approx_eq!(f32, atom1.get_force().x, 0.0);
        assert_approx_eq!(f32, atom1.get_force().y, 0.0);
        assert_approx_eq!(f32, atom1.get_force().z, 0.0);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_approx_eq!(f32, atom2.get_position().x, 7.06);
        assert_approx_eq!(f32, atom2.get_position().y, 0.42);
        assert_approx_eq!(f32, atom2.get_position().z, 9.38);

        assert_approx_eq!(f32, atom2.get_velocity().x, 0.0712);
        assert_approx_eq!(f32, atom2.get_velocity().y, 0.2294);
        assert_approx_eq!(f32, atom2.get_velocity().z, -0.1673);

        assert_approx_eq!(f32, atom2.get_force().x, 0.0);
        assert_approx_eq!(f32, atom2.get_force().y, 0.0);
        assert_approx_eq!(f32, atom2.get_force().z, 0.0);

        // last frame
        system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .last();

        assert_eq!(system.get_simulation_step(), 50000);
        assert_eq!(system.get_precision(), 100);
        assert_approx_eq!(f32, system.get_simulation_time(), 1000.0);
        assert_approx_eq!(f32, system.get_box_as_ref().x, 13.02659);
        assert_approx_eq!(f32, system.get_box_as_ref().y, 13.02659);
        assert_approx_eq!(f32, system.get_box_as_ref().z, 11.250414);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_approx_eq!(f32, atom1.get_position().x, 5.67);
        assert_approx_eq!(f32, atom1.get_position().y, 6.31);
        assert_approx_eq!(f32, atom1.get_position().z, 6.96);

        assert_approx_eq!(f32, atom1.get_velocity().x, -0.0683);
        assert_approx_eq!(f32, atom1.get_velocity().y, 0.1133);
        assert_approx_eq!(f32, atom1.get_velocity().z, 0.0005);

        assert_approx_eq!(f32, atom1.get_force().x, 0.0);
        assert_approx_eq!(f32, atom1.get_force().y, 0.0);
        assert_approx_eq!(f32, atom1.get_force().z, 0.0);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_approx_eq!(f32, atom2.get_position().x, 8.36);
        assert_approx_eq!(f32, atom2.get_position().y, 1.18);
        assert_approx_eq!(f32, atom2.get_position().z, 0.71);

        assert_approx_eq!(f32, atom2.get_velocity().x, 0.0712);
        assert_approx_eq!(f32, atom2.get_velocity().y, 0.2294);
        assert_approx_eq!(f32, atom2.get_velocity().z, -0.1673);

        assert_approx_eq!(f32, atom2.get_force().x, 0.0);
        assert_approx_eq!(f32, atom2.get_force().y, 0.0);
        assert_approx_eq!(f32, atom2.get_force().z, 0.0);
    }

    #[test]
    fn read_xtc_iter() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        for (i, raw) in system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .enumerate()
        {
            let frame = raw.unwrap();

            assert_eq!(frame.get_simulation_time() as usize, i * 100);
        }
    }

    #[test]
    fn read_xtc_range_full() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        for (raw_frame1, raw_frame2) in system1
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .zip(
                system2
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(0.0, 10000.0)
                    .unwrap(),
            )
        {
            let frame1 = raw_frame1.unwrap();
            let frame2 = raw_frame2.unwrap();

            for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                assert_eq!(atom1.get_residue_number(), atom2.get_residue_number());
                assert_eq!(atom1.get_residue_name(), atom2.get_residue_name());
                assert_eq!(atom1.get_atom_number(), atom2.get_atom_number());
                assert_eq!(atom1.get_atom_name(), atom2.get_atom_name());
                assert_eq!(atom1.get_chain(), atom2.get_chain());

                assert_approx_eq!(f32, atom1.get_position().x, atom2.get_position().x);
                assert_approx_eq!(f32, atom1.get_position().y, atom2.get_position().y);
                assert_approx_eq!(f32, atom1.get_position().z, atom2.get_position().z);

                assert_approx_eq!(f32, atom1.get_velocity().x, atom2.get_velocity().x);
                assert_approx_eq!(f32, atom1.get_velocity().y, atom2.get_velocity().y);
                assert_approx_eq!(f32, atom1.get_velocity().z, atom2.get_velocity().z);

                assert_approx_eq!(f32, atom1.get_force().x, atom2.get_force().x);
                assert_approx_eq!(f32, atom1.get_force().y, atom2.get_force().y);
                assert_approx_eq!(f32, atom1.get_force().z, atom2.get_force().z);
            }
        }
    }

    #[test]
    fn read_xtc_range() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        let mut i = 0;

        for (raw_frame1, raw_frame2) in system1
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .skip(3)
            .take(6)
            .zip(
                system2
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(300.0, 800.0)
                    .unwrap(),
            )
        {
            i += 1;

            let frame1 = raw_frame1.unwrap();
            let frame2 = raw_frame2.unwrap();

            for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                assert_eq!(atom1.get_residue_number(), atom2.get_residue_number());
                assert_eq!(atom1.get_residue_name(), atom2.get_residue_name());
                assert_eq!(atom1.get_atom_number(), atom2.get_atom_number());
                assert_eq!(atom1.get_atom_name(), atom2.get_atom_name());
                assert_eq!(atom1.get_chain(), atom2.get_chain());

                assert_approx_eq!(f32, atom1.get_position().x, atom2.get_position().x);
                assert_approx_eq!(f32, atom1.get_position().y, atom2.get_position().y);
                assert_approx_eq!(f32, atom1.get_position().z, atom2.get_position().z);

                assert_approx_eq!(f32, atom1.get_velocity().x, atom2.get_velocity().x);
                assert_approx_eq!(f32, atom1.get_velocity().y, atom2.get_velocity().y);
                assert_approx_eq!(f32, atom1.get_velocity().z, atom2.get_velocity().z);

                assert_approx_eq!(f32, atom1.get_force().x, atom2.get_force().x);
                assert_approx_eq!(f32, atom1.get_force().y, atom2.get_force().y);
                assert_approx_eq!(f32, atom1.get_force().z, atom2.get_force().z);
            }
        }

        assert_eq!(i, 6);
    }

    #[test]
    fn read_xtc_range_full_iter() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        for (i, raw) in system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .with_range(0.0, 1000.0)
            .unwrap()
            .enumerate()
        {
            let frame = raw.unwrap();

            assert_eq!(frame.get_simulation_time() as usize, i * 100);
        }
    }

    #[test]
    fn read_xtc_range_iter() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        for (i, raw) in system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .with_range(300.0, 800.0)
            .unwrap()
            .enumerate()
        {
            let frame = raw.unwrap();

            assert_eq!(frame.get_simulation_time() as usize, (i + 3) * 100);
        }
    }

    #[test]
    fn read_xtc_range_negative() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .with_range(-300.0, 800.0)
        {
            Ok(_) => panic!("Iterator should not have been constructed."),
            Err(ReadXdrError::TimeRangeNegative(_)) => (),
            Err(e) => panic!("Incorrect error type {} returned", e),
        }

        match system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .with_range(300.0, -800.0)
        {
            Ok(_) => panic!("Iterator should not have been constructed."),
            Err(ReadXdrError::TimeRangeNegative(_)) => (),
            Err(e) => panic!("Incorrect error type {} returned", e),
        }
    }

    #[test]
    fn read_xtc_range_end_start() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .with_range(800.0, 300.0)
        {
            Ok(_) => panic!("Iterator should not have been constructed."),
            Err(ReadXdrError::InvalidTimeRange(_, _)) => (),
            Err(e) => panic!("Incorrect error type {} returned", e),
        }
    }

    #[test]
    fn read_xtc_range_start_not_found() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .with_range(12000.0, 20000.0)
        {
            Ok(_) => panic!("Iterator should not have been constructed."),
            Err(ReadXdrError::StartNotFound(_)) => (),
            Err(e) => panic!("Incorrect error type {} returned", e),
        }
    }

    #[test]
    fn read_xtc_unmatching() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();

        match system.xtc_iter("test_files/short_trajectory.xtc") {
            Err(ReadXdrError::AtomsNumberMismatch(_)) => (),
            _ => panic!("XTC file should not be valid."),
        }
    }

    #[test]
    fn read_xtc_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.xtc_iter("test_files/nonexistent.xtc") {
            Err(ReadXdrError::FileNotFound(_)) => (),
            _ => panic!("XTC file should not exist."),
        }
    }

    #[test]
    fn write_xtc() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        let mut writer = XtcWriter::new(&system, path_to_output).unwrap();

        for _ in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            writer.write_frame().unwrap();
        }

        // we must close the file, otherwise metadata do not get updated
        drop(writer);

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_invalid_path() {
        let system = System::from_file("test_files/example.gro").unwrap();

        match XtcWriter::new(&system, "test_files/nonexistent/output.xtc") {
            Err(WriteXdrError::CouldNotCreate(_)) => (),
            _ => panic!("Output XTC file should not have been created."),
        }
    }

    #[test]
    fn write_group_xtc() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        let mut writer = XtcGroupWriter::new(&system, "Protein", path_to_output).unwrap();

        for _ in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            writer.write_frame().unwrap();
        }

        // we must close the file, otherwise metadata do not get updated
        drop(writer);

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory_protein.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group_xtc_all() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        let mut writer = XtcGroupWriter::new(&system, "all", path_to_output).unwrap();

        for _ in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            writer.write_frame().unwrap();
        }

        // we must close the file, otherwise metadata do not get updated
        drop(writer);

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group_xtc_replace() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        let mut writer = XtcGroupWriter::new(&system, "Protein", path_to_output).unwrap();

        // replace the protein group with something else; this should not change the output of the XtcGroupWriter
        if let Ok(_) = system.group_create("Protein", "serial 1") {
            panic!("Function should return warning but it did not.");
        }

        for _ in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            writer.write_frame().unwrap();
        }

        // we must close the file, otherwise metadata do not get updated
        drop(writer);

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory_protein.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group_xtc_remove() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        let mut writer = XtcGroupWriter::new(&system, "Protein", path_to_output).unwrap();

        // remove the protein group from the system; this should not change the output of the XtcGroupWriter
        unsafe {
            let val = system.get_groups_as_ref_mut().remove("Protein").unwrap();
            assert_eq!(val.get_atom_ranges(), writer.group.get_atom_ranges());
            assert!(!system.group_exists("Protein"));
        }

        for _ in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            writer.write_frame().unwrap();
        }

        // we must close the file, otherwise metadata do not get updated
        drop(writer);

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory_protein.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group_xtc_nonexistent() {
        let system = System::from_file("test_files/example.gro").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        match XtcGroupWriter::new(&system, "Protein", path_to_output) {
            Err(WriteXdrError::GroupNotFound(g)) => assert_eq!(g, "Protein"),
            _ => panic!("Output XTC file should not have been created."),
        }
    }
}
