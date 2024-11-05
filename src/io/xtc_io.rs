// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of functions for reading and writing xtc files.

use std::marker::PhantomData;
use std::os::raw::{c_float, c_int};
use std::path::Path;

use crate::errors::{ReadTrajError, TrajError, WriteTrajError};
use crate::io::traj_cat::TrajConcatenator;
use crate::io::traj_read::{
    FrameData, FrameDataTime, TrajRangeRead, TrajRead, TrajReadOpen, TrajStepRead, TrajStepTimeRead,
};
use crate::io::xdrfile::{self, CXdrFile, OpenMode, XdrFile};
use crate::prelude::{AtomIterator, TrajReader};
use crate::structures::group::Group;
use crate::structures::{simbox::SimBox, vector3d::Vector3D};
use crate::system::System;

use super::traj_write::{PrivateTrajWrite, TrajWrite};

/**************************/
/*       READING XTC      */
/**************************/

/// Structure containing data read from an xtc file.
#[derive(Debug)]
pub struct XtcFrameData {
    step: c_int,
    time: c_float,
    simbox: SimBox,
    precision: c_float,
    coordinates: Vec<[f32; 3]>,
}

impl FrameData for XtcFrameData {
    type TrajFile = CXdrFile;

    /// Read `XtcFrameData` from an xtc frame.
    ///
    /// ## Returns
    /// - `None` if the file has been read completely
    /// - `Some(XtcFrameData)` if the frame has been successfully read.
    /// - `Some(ReadTrajError)` if the frame could not be read.
    fn from_frame(
        xdrfile: &mut Self::TrajFile,
        system: &System,
    ) -> Option<Result<Self, ReadTrajError>> {
        let mut step: c_int = 0;
        let mut time: c_float = 0.0;
        let mut boxvector: [[c_float; 3usize]; 3usize] = [[0.0; 3]; 3];
        let mut precision: c_float = 0.0;
        let mut coordinates = vec![[0.0, 0.0, 0.0]; system.get_n_atoms()];

        unsafe {
            let return_code = xdrfile::read_xtc(
                xdrfile,
                system.get_n_atoms() as c_int,
                &mut step,
                &mut time,
                &mut boxvector,
                coordinates.as_mut_ptr(),
                &mut precision,
            );

            match return_code {
                // reading successful and there is more to read
                0 => {
                    let simbox = match xdrfile::matrix2simbox(boxvector) {
                        Ok(x) => x,
                        Err(e) => return Some(Err(e)),
                    };

                    Some(Ok(XtcFrameData {
                        step,
                        time,
                        simbox,
                        precision,
                        coordinates,
                    }))
                }
                // file is completely read
                11 => None,
                // error occured
                _ => Some(Err(ReadTrajError::FrameNotFound)),
            }
        }
    }

    /// Update the `System` structure based on data from `XtcFrameData`.
    fn update_system(self, system: &mut System) {
        unsafe {
            for (i, atom) in system.get_atoms_mut().iter_mut().enumerate() {
                atom.set_position(Vector3D::from(*self.coordinates.get_unchecked(i)));
                atom.reset_velocity();
                atom.reset_force();
            }

            // update the system
            system.set_simulation_step(self.step as u64);
            system.set_simulation_time(self.time);
            system.set_box(self.simbox);
            system.set_precision(self.precision as u64);
        }
    }
}

impl FrameDataTime for XtcFrameData {
    fn get_time(&self) -> f32 {
        self.time
    }
}

/// Iterator over an xtc file.
pub struct XtcReader<'a> {
    system: *mut System,
    xtc: XdrFile,
    phantom: PhantomData<&'a mut System>,
}

impl<'a> TrajRead<'a> for XtcReader<'a> {
    type FrameData = XtcFrameData;

    fn get_system(&mut self) -> *mut System {
        self.system
    }

    fn get_file_handle(&mut self) -> &mut CXdrFile {
        unsafe { self.xtc.handle.as_mut().unwrap() }
    }
}

impl<'a> TrajReadOpen<'a> for XtcReader<'a> {
    /// Create an iterator over an xtc file.
    fn new(system: &'a mut System, filename: impl AsRef<Path>) -> Result<XtcReader, ReadTrajError> {
        let n_atoms = system.get_n_atoms();

        // sanity check the number of atoms
        match XdrFile::check_xtc(filename.as_ref(), n_atoms) {
            Err(e) => return Err(e),
            Ok(false) => {
                return Err(ReadTrajError::AtomsNumberMismatch(Box::from(
                    filename.as_ref(),
                )))
            }
            Ok(true) => (),
        };

        // open the xtc file and save the handle to it
        let xtc = match XdrFile::open_xdr(filename.as_ref(), OpenMode::Read) {
            Ok(x) => x,
            Err(TrajError::FileNotFound(x)) => return Err(ReadTrajError::FileNotFound(x)),
            Err(TrajError::InvalidPath(x)) => return Err(ReadTrajError::InvalidPath(x)),
        };

        let xtc_reader = XtcReader {
            system,
            xtc,
            phantom: PhantomData,
        };

        Ok(xtc_reader)
    }
}

impl<'a> TrajRangeRead<'a> for XtcReader<'a> {
    fn jump_to_start(&mut self, start_time: f32) -> Result<(), ReadTrajError> {
        unsafe {
            if xdrfile::xtc_jump_to_start(self.get_file_handle(), start_time as c_float) != 0 {
                Err(ReadTrajError::StartNotFound(start_time.to_string()))
            } else {
                Ok(())
            }
        }
    }
}

impl<'a> TrajStepRead<'a> for XtcReader<'a> {
    fn skip_frame(&mut self) -> Result<bool, ReadTrajError> {
        unsafe {
            match xdrfile::xtc_skip_frame(self.get_file_handle()) {
                0 => Ok(true),
                1 => Err(ReadTrajError::SkipFailed),
                2 => Ok(false),
                number => panic!(
                    "FATAL GROAN ERROR | XtcReader::skip_frame | `xdrfile::xtc_skip_frame` returned an unsupported number '{}'",
                    number
                ),
            }
        }
    }
}

impl<'a> TrajStepTimeRead<'a> for XtcReader<'a> {
    fn skip_frame_time(&mut self) -> Result<Option<f32>, ReadTrajError> {
        unsafe {
            let mut time: c_float = 0.0;

            match xdrfile::xtc_skip_frame_with_time(self.get_file_handle(), &mut time as *mut c_float) {
                0 => Ok(Some(time)),
                1 => Err(ReadTrajError::SkipFailed),
                2 => Ok(None),
                number => panic!(
                    "FATAL GROAN ERROR | XtcReader::skip_frame_time | `xdrfile::xtc_skip_frame_with_time` returned an unsupported number '{}'",
                    number
                ),
            }
        }
    }
}

/// ## Methods for reading xtc files.
impl System {
    /// Create an `XtcReader` structure which is an iterator over an xtc file.
    ///
    /// ## Returns
    /// `TrajReader<XtcReader>` if the xtc file exists and matches the structure file.
    /// Else returns `ReadTrajError`.
    ///
    /// ## Examples
    /// Iterating through an xtc trajectory and calculating
    /// and printing the current center of geometry of the system.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
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
    /// # use groan_rs::prelude::*;
    /// # use groan_rs::errors::ReadTrajError;
    /// #
    /// fn example_fn() -> Result<(), ReadTrajError> {
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
    /// The `with_range` method is very efficient for the xtc files
    /// as the particle coordinates from the xtc frames outside
    /// of the range will not be read.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use groan_rs::errors::ReadTrajError;
    /// #
    /// fn example_fn() -> Result<(), ReadTrajError> {
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
    /// Furthermore, you can efficiently skip over some frames of the trajectory.
    /// Here, only every 10th frame of the trajectory will be read.
    /// The `with_step` method is very efficient for the xtc files as the
    /// particle coordinates from the skipped over frames will not be read.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use groan_rs::errors::ReadTrajError;
    /// #
    /// fn example_fn() -> Result<(), ReadTrajError> {
    ///     let mut system = System::from_file("system.gro").unwrap();
    ///
    ///     for raw_frame in system
    ///         .xtc_iter("trajectory.xtc")?
    ///         .with_step(10)?
    ///     {
    ///         let frame = raw_frame?;
    ///         println!("{:?}", frame.group_get_center("all"));
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    ///
    /// Finally, you can combine `with_range` and `with_step` to
    /// iterate only through a specific time range of the trajectory
    /// while reading only every `step`th frame.
    /// Here, only every other frame will be read and the iteration
    /// will start at time 10 ns and end at time 100 ns.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use groan_rs::errors::ReadTrajError;
    /// #
    /// fn example_fn() -> Result<(), ReadTrajError> {
    ///     let mut system = System::from_file("system.gro").unwrap();
    ///
    ///     for raw_frame in system
    ///         .xtc_iter("trajectory.xtc")?
    ///         .with_range(10_000.0, 100_000.0)?
    ///         .with_step(2)?
    ///     {
    ///         let frame = raw_frame?;
    ///         println!("{:?}", frame.group_get_center("all"));
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    ///
    /// ## Notes
    /// - The function checks whether the number of atoms in the system corresponds to the number of atoms in the xtc file.
    /// - The `System` structure is modified while iterating through the xtc file.
    /// - The `velocity` and `force` information is set to `None` for all atoms as it is not available in the xtc file.
    pub fn xtc_iter(
        &mut self,
        filename: impl AsRef<Path>,
    ) -> Result<TrajReader<XtcReader>, ReadTrajError> {
        Ok(TrajReader::wrap_traj(XtcReader::new(self, filename)?))
    }

    /// Iterate through multiple xtc files.
    /// Any duplicate frames at the boundaries of the trajectories are skipped.
    ///
    /// ## Example
    /// Iterate through multiple xtc files and calculate and print the current center of geometry of the system.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use groan_rs::errors::ReadTrajError;
    /// #
    /// fn example_fn() -> Result<(), ReadTrajError> {
    ///     let mut system = System::from_file("system.gro").unwrap();
    ///     let trajectories = vec!["md0001.xtc", "md0002.xtc", "md0003.xtc"];     
    ///
    ///     for raw_frame in system.xtc_cat_iter(&trajectories).unwrap() {
    ///         let frame = raw_frame?;
    ///         println!("{:?}", frame.group_get_center("all"));
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    ///
    /// You can also iterate through only part of the trajectory, only read every Nth trajectory frame,
    /// and print progress of the trajectory reading. For more information, see [`System::xtc_iter`](`System::xtc_iter`).
    ///
    /// For more information about trajectory concatenation, see [`System::traj_cat_iter`](`System::traj_cat_iter`).
    pub fn xtc_cat_iter<'a>(
        &'a mut self,
        filenames: &[impl AsRef<Path>],
    ) -> Result<TrajReader<'a, TrajConcatenator<'a, XtcReader>>, ReadTrajError> {
        self.traj_cat_iter::<XtcReader>(filenames)
    }
}

/**************************/
/*       WRITING XTC      */
/**************************/

impl System {
    #[inline(always)]
    pub fn xtc_writer_init(&mut self, filename: impl AsRef<Path>) -> Result<(), WriteTrajError> {
        self.traj_writer_init::<XtcWriter>(filename)
    }

    #[inline(always)]
    pub fn xtc_group_writer_init(
        &mut self,
        filename: impl AsRef<Path>,
        group: &str,
    ) -> Result<(), WriteTrajError> {
        self.traj_group_writer_init::<XtcWriter>(filename, group)
    }
}

pub struct XtcWriter {
    xtc: XdrFile,
    // deep copy of the group from `System`
    group: Group,
}

impl TrajWrite for XtcWriter {}

impl PrivateTrajWrite for XtcWriter {
    /// Open a new xtc file for writing.
    fn new(
        system: &System,
        filename: impl AsRef<Path>,
        group: Option<&str>,
    ) -> Result<Self, WriteTrajError>
    where
        Self: Sized,
    {
        // create the xtc file and save a handle to it
        let xtc = match XdrFile::open_xdr(filename.as_ref(), OpenMode::Write) {
            Ok(x) => x,
            Err(TrajError::FileNotFound(x)) => return Err(WriteTrajError::CouldNotCreate(x)),
            Err(TrajError::InvalidPath(x)) => return Err(WriteTrajError::InvalidPath(x)),
        };

        // get the requested group from the system or use `all`
        let group = match group {
            Some(x) => system
                .get_groups()
                .get(x)
                .ok_or_else(|| WriteTrajError::GroupNotFound(x.to_owned()))?
                .clone(),
            None => system
                .get_groups()
                .get("all")
                .expect("FATAL GROAN ERROR | XtcWriter::new | Group `all` should exist.")
                .clone(),
        };

        Ok(Self { xtc, group })
    }

    /// Write the current state of the system into an open xtc file.
    fn write_frame(&mut self, system: &System) -> Result<(), WriteTrajError> {
        let n_atoms = self.group.get_n_atoms();

        // prepare coordinate matrix
        let iterator =
            AtomIterator::new(system.get_atoms(), self.group.get_atoms(), system.get_box());
        let mut coordinates = vec![[0.0, 0.0, 0.0]; n_atoms];
        for atom in iterator {
            if let Some(pos) = atom.get_position() {
                coordinates[atom.get_index()] = [pos.x, pos.y, pos.z];
            }
        }

        // write the xtc frame
        let return_code = unsafe {
            xdrfile::write_xtc(
                self.xtc.handle,
                n_atoms as c_int,
                system.get_simulation_step() as i32,
                system.get_simulation_time(),
                &mut xdrfile::simbox2matrix(system.get_box()),
                coordinates.as_mut_ptr(),
                system.get_precision() as f32,
            )
        };

        if return_code != 0 {
            Err(WriteTrajError::CouldNotWrite)
        } else {
            Ok(())
        }
    }
}

/******************************/
/*     PRIVATE FUNCTIONS      */
/******************************/

impl XdrFile {
    /// Check that the number of atoms in an unopened xtc file matches the expected number.
    fn check_xtc(filename: impl AsRef<Path>, n_atoms: usize) -> Result<bool, ReadTrajError> {
        unsafe {
            let c_path = match xdrfile::path2cstring(filename.as_ref()) {
                Ok(x) => x,
                Err(_) => return Err(ReadTrajError::InvalidPath(Box::from(filename.as_ref()))),
            };

            let mut xtc_atoms: c_int = 0;

            if xdrfile::read_xtc_natoms(c_path.as_ptr(), &mut xtc_atoms) != 0 {
                // reading the file failed
                return Err(ReadTrajError::FileNotFound(Box::from(filename.as_ref())));
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

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::assert_approx_eq;
    use std::fs::File;
    use tempfile::NamedTempFile;

    use crate::test_utilities::utilities::{compare_atoms, compare_box};

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
        let simbox = system.get_box().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.034535);
        assert_approx_eq!(f32, simbox.y, 13.034535);
        assert_approx_eq!(f32, simbox.z, 11.228164);

        let atom1 = &system.get_atoms()[0];
        let atom2 = &system.get_atoms()[16843];

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_approx_eq!(f32, atom1.get_position().unwrap().x, 5.97);
        assert_approx_eq!(f32, atom1.get_position().unwrap().y, 7.03);
        assert_approx_eq!(f32, atom1.get_position().unwrap().z, 7.69);

        assert_eq!(atom1.get_velocity(), None);

        assert_eq!(atom1.get_force(), None);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_approx_eq!(f32, atom2.get_position().unwrap().x, 7.06);
        assert_approx_eq!(f32, atom2.get_position().unwrap().y, 0.42);
        assert_approx_eq!(f32, atom2.get_position().unwrap().z, 9.38);

        assert_eq!(atom2.get_velocity(), None);

        assert_eq!(atom2.get_force(), None);

        // last frame
        system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .last();

        assert_eq!(system.get_simulation_step(), 50000);
        assert_eq!(system.get_precision(), 100);
        assert_approx_eq!(f32, system.get_simulation_time(), 1000.0);
        let simbox = system.get_box().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.02659);
        assert_approx_eq!(f32, simbox.y, 13.02659);
        assert_approx_eq!(f32, simbox.z, 11.250414);

        let atom1 = &system.get_atoms()[0];
        let atom2 = &system.get_atoms()[16843];

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_approx_eq!(f32, atom1.get_position().unwrap().x, 5.67);
        assert_approx_eq!(f32, atom1.get_position().unwrap().y, 6.31);
        assert_approx_eq!(f32, atom1.get_position().unwrap().z, 6.96);

        assert_eq!(atom1.get_velocity(), None);

        assert_eq!(atom1.get_force(), None);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_approx_eq!(f32, atom2.get_position().unwrap().x, 8.36);
        assert_approx_eq!(f32, atom2.get_position().unwrap().y, 1.18);
        assert_approx_eq!(f32, atom2.get_position().unwrap().z, 0.71);

        assert_eq!(atom2.get_velocity(), None);

        assert_eq!(atom2.get_force(), None);
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
                compare_atoms(atom1, atom2);
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
                compare_atoms(atom1, atom2);
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
            Err(ReadTrajError::TimeRangeNegative(_)) => (),
            Err(e) => panic!("Incorrect error type {} returned", e),
        }

        match system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .with_range(300.0, -800.0)
        {
            Ok(_) => panic!("Iterator should not have been constructed."),
            Err(ReadTrajError::TimeRangeNegative(_)) => (),
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
            Err(ReadTrajError::InvalidTimeRange(_, _)) => (),
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
            Err(ReadTrajError::StartNotFound(_)) => (),
            Err(e) => panic!("Incorrect error type {} returned", e),
        }
    }

    #[test]
    fn read_xtc_step_1() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        for (raw1, raw2) in system1
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .with_step(1)
            .unwrap()
            .zip(system2.xtc_iter("test_files/short_trajectory.xtc").unwrap())
        {
            let frame1 = raw1.unwrap();
            let frame2 = raw2.unwrap();

            for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                compare_atoms(atom1, atom2);
            }
        }
    }

    #[test]
    fn read_xtc_step_3() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        for (raw1, raw2) in system1
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .with_step(3)
            .unwrap()
            .zip(
                system2
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .step_by(3),
            )
        {
            let frame1 = raw1.unwrap();
            let frame2 = raw2.unwrap();

            for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                compare_atoms(atom1, atom2);
            }
        }
    }

    #[test]
    fn read_xtc_step_23() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        for (raw1, raw2) in system1
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .with_step(23)
            .unwrap()
            .zip(
                system2
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .step_by(23),
            )
        {
            let frame1 = raw1.unwrap();
            let frame2 = raw2.unwrap();

            for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                compare_atoms(atom1, atom2);
            }
        }
    }

    #[test]
    fn read_xtc_step_0() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .with_step(0)
        {
            Ok(_) => panic!("Should have failed."),
            Err(ReadTrajError::InvalidStep(s)) => assert_eq!(s, 0),
            Err(e) => panic!("Incorrect error type {} returned.", e),
        }
    }

    #[test]
    fn read_xtc_range_step() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        let mut i = 0;

        for (raw_frame1, raw_frame2) in system1
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .skip(3)
            .step_by(2)
            .zip(
                system2
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(300.0, 800.0)
                    .unwrap()
                    .with_step(2)
                    .unwrap(),
            )
        {
            i += 1;

            let frame1 = raw_frame1.unwrap();
            let frame2 = raw_frame2.unwrap();

            for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                compare_atoms(atom1, atom2);
            }
        }

        assert_eq!(i, 3);
    }

    /// Tests that the order of `with_step` and `with_range` does not matter.
    #[test]
    fn read_xtc_range_step_step_range() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        let mut i = 0;

        for (raw_frame1, raw_frame2) in system1
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .with_step(2)
            .unwrap()
            .with_range(300.0, 800.0)
            .unwrap()
            .zip(
                system2
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(300.0, 800.0)
                    .unwrap()
                    .with_step(2)
                    .unwrap(),
            )
        {
            i += 1;

            let frame1 = raw_frame1.unwrap();
            let frame2 = raw_frame2.unwrap();

            for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                compare_atoms(atom1, atom2);
            }
        }

        assert_eq!(i, 3);
    }

    #[test]
    fn read_xtc_unmatching() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();

        match system.xtc_iter("test_files/short_trajectory.xtc") {
            Err(ReadTrajError::AtomsNumberMismatch(_)) => (),
            _ => panic!("XTC file should not be valid."),
        }
    }

    #[test]
    fn read_xtc_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.xtc_iter("test_files/nonexistent.xtc") {
            Err(ReadTrajError::FileNotFound(_)) => (),
            _ => panic!("XTC file should not exist."),
        }
    }

    #[test]
    fn read_xtc_triclinic() {
        let mut system = System::from_file("test_files/triclinic.gro").unwrap();

        // second frame
        let frame = system
            .xtc_iter("test_files/triclinic_trajectory.xtc")
            .unwrap()
            .nth(1)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 5000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 100.0);

        let simbox = frame.get_box().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 5.2868834);
        assert_approx_eq!(f32, simbox.v2y, 4.7799735);
        assert_approx_eq!(f32, simbox.v3z, 2.2256064);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 0.8428372);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, 1.0159061);
        assert_approx_eq!(f32, simbox.v3y, -1.6872015);

        // last frame
        let frame = system
            .xtc_iter("test_files/triclinic_trajectory.xtc")
            .unwrap()
            .nth(10)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 50000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 1000.0);

        let simbox = frame.get_box().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 5.2712817);
        assert_approx_eq!(f32, simbox.v2y, 4.7658677);
        assert_approx_eq!(f32, simbox.v3z, 2.1743093);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 0.84035);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, 1.0129081);
        assert_approx_eq!(f32, simbox.v3y, -1.6822226);
    }

    #[test]
    fn read_xtc_octahedron() {
        let mut system = System::from_file("test_files/octahedron.gro").unwrap();

        // second frame
        let frame = system
            .xtc_iter("test_files/octahedron_trajectory.xtc")
            .unwrap()
            .nth(1)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 5000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 100.0);

        let simbox = frame.get_box().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 6.266603);
        assert_approx_eq!(f32, simbox.v2y, 5.908211);
        assert_approx_eq!(f32, simbox.v3z, 5.1106043);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 2.0888677);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, -2.0888677);
        assert_approx_eq!(f32, simbox.v3y, 2.9541006);

        // last frame
        let frame = system
            .xtc_iter("test_files/octahedron_trajectory.xtc")
            .unwrap()
            .nth(10)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 50000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 1000.0);

        let simbox = frame.get_box().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 6.2004085);
        assert_approx_eq!(f32, simbox.v2y, 5.8458023);
        assert_approx_eq!(f32, simbox.v3z, 5.0840497);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 2.066803);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, -2.066803);
        assert_approx_eq!(f32, simbox.v3y, 2.9228961);
    }

    #[test]
    fn read_xtc_dodecahedron() {
        let mut system = System::from_file("test_files/dodecahedron.gro").unwrap();

        // second frame
        let frame = system
            .xtc_iter("test_files/dodecahedron_trajectory.xtc")
            .unwrap()
            .nth(1)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 5000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 100.0);

        let simbox = frame.get_box().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 6.260709);
        assert_approx_eq!(f32, simbox.v2y, 6.260709);
        assert_approx_eq!(f32, simbox.v3z, 4.4316807);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 0.0000000);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, 3.1303544);
        assert_approx_eq!(f32, simbox.v3y, 3.1303544);

        // last frame
        let frame = system
            .xtc_iter("test_files/dodecahedron_trajectory.xtc")
            .unwrap()
            .nth(10)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 50000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 1000.0);

        let simbox = frame.get_box().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 6.2197995);
        assert_approx_eq!(f32, simbox.v2y, 6.2197995);
        assert_approx_eq!(f32, simbox.v3z, 4.4066653);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 0.0000000);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, 3.1098998);
        assert_approx_eq!(f32, simbox.v3y, 3.1098998);
    }

    #[test]
    fn cat_xtc() {
        let mut system_single = System::from_file("test_files/example.gro").unwrap();
        let mut system_cat = System::from_file("test_files/example.gro").unwrap();

        let traj_single = system_single
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap();
        let traj_cat = system_cat
            .xtc_cat_iter(&[
                "test_files/split/traj1.xtc",
                "test_files/split/traj2.xtc",
                "test_files/split/traj3.xtc",
                "test_files/split/traj4.xtc",
                "test_files/split/traj5.xtc",
                "test_files/split/traj6.xtc",
            ])
            .unwrap();

        for (frame_single, frame_cat) in traj_single.zip(traj_cat) {
            let frame_single = frame_single.unwrap();
            let frame_cat = frame_cat.unwrap();

            assert_approx_eq!(
                f32,
                frame_single.get_simulation_time(),
                frame_cat.get_simulation_time()
            );
            assert_eq!(
                frame_single.get_simulation_step(),
                frame_cat.get_simulation_step()
            );

            compare_box(
                frame_single.get_box().unwrap(),
                frame_cat.get_box().unwrap(),
            );

            for (atom_single, atom_cat) in frame_single.atoms_iter().zip(frame_cat.atoms_iter()) {
                compare_atoms(atom_single, atom_cat);
            }
        }
    }

    #[test]
    fn write_xtc() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        system.xtc_writer_init(path_to_output).unwrap();

        for frame in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            let frame = frame.unwrap();

            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_invalid_path() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.xtc_writer_init("test_files/nonexistent/output.xtc") {
            Err(WriteTrajError::CouldNotCreate(_)) => (),
            _ => panic!("Output XTC file should not have been created."),
        }
    }

    #[test]
    fn write_group_xtc() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        system
            .xtc_group_writer_init(path_to_output, "Protein")
            .unwrap();

        for frame in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            let frame = frame.unwrap();
            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory_protein.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group_xtc_all() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        system.xtc_group_writer_init(path_to_output, "all").unwrap();

        for frame in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            let frame = frame.unwrap();

            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group_xtc_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        match system.xtc_group_writer_init(path_to_output, "Protein") {
            Err(WriteTrajError::GroupNotFound(g)) => assert_eq!(g, "Protein"),
            _ => panic!("Output XTC file should not have been created."),
        }
    }

    #[test]
    fn write_group_xtc_replace() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        system
            .xtc_group_writer_init(path_to_output, "Protein")
            .unwrap();

        // replace the protein group with something else; this should not change the output of the trajectory writing
        if system.group_create("Protein", "serial 1").is_ok() {
            panic!("Function should return warning but it did not.");
        }

        for frame in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            let frame = frame.unwrap();
            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

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

        system
            .xtc_group_writer_init(path_to_output, "Protein")
            .unwrap();

        // remove the `Protein` group; this should not change the output of the trajectory writing
        system.group_remove("Protein").unwrap();

        for frame in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            let frame = frame.unwrap();
            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory_protein.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_xtc_triclinic() {
        let mut system = System::from_file("test_files/triclinic.gro").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        system.xtc_writer_init(path_to_output).unwrap();

        for frame in system
            .xtc_iter("test_files/triclinic_trajectory.xtc")
            .unwrap()
        {
            let frame = frame.unwrap();

            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/triclinic_trajectory.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_xtc_octahedron() {
        let mut system = System::from_file("test_files/octahedron.gro").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        system.xtc_writer_init(path_to_output).unwrap();

        for frame in system
            .xtc_iter("test_files/octahedron_trajectory.xtc")
            .unwrap()
        {
            let frame = frame.unwrap();

            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/octahedron_trajectory.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_xtc_dodecahedron() {
        let mut system = System::from_file("test_files/dodecahedron.gro").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        system.xtc_writer_init(path_to_output).unwrap();

        for frame in system
            .xtc_iter("test_files/dodecahedron_trajectory.xtc")
            .unwrap()
        {
            let frame = frame.unwrap();

            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/dodecahedron_trajectory.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
}
