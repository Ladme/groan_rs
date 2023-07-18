// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of functions for reading and writing xtc files.

use std::ffi::{CString, NulError};
use std::marker::PhantomData;
use std::os::raw::{c_float, c_int};
use std::path::Path;

use crate::atom::Atom;
use crate::c_xdrfile;
use crate::errors::{ReadXtcError, WriteXtcError, XtcError};
use crate::simbox::SimBox;
use crate::system::System;
use crate::vector3d::Vector3D;

/**************************/
/*       READING XTC      */
/**************************/

/// Iterator over an xtc file.
pub struct XtcReader<'a> {
    system: *mut System,
    xtc: Xdrfile,
    phantom: PhantomData<&'a mut System>,
}

impl<'a> XtcReader<'a> {
    /// Create an iterator over an xtc file.
    fn new(system: &'a mut System, filename: impl AsRef<Path>) -> Result<Self, ReadXtcError> {
        let n_atoms = system.get_n_atoms();

        // sanity check the number of atoms
        match Xdrfile::check_xtc(filename.as_ref(), n_atoms) {
            Err(e) => return Err(e),
            Ok(false) => {
                return Err(ReadXtcError::AtomsNumberMismatch(Box::from(
                    filename.as_ref(),
                )))
            }
            Ok(true) => (),
        };

        // open the xtc file and save the handle to it
        let xtc = match Xdrfile::open_xtc(filename.as_ref(), OpenMode::Read) {
            Ok(x) => x,
            Err(XtcError::FileNotFound(x)) => return Err(ReadXtcError::FileNotFound(x)),
            Err(XtcError::InvalidPath(x)) => return Err(ReadXtcError::InvalidPath(x)),
        };

        Ok(XtcReader {
            system,
            xtc,
            phantom: PhantomData,
        })
    }
}

impl<'a> Iterator for XtcReader<'a> {
    type Item = Result<&'a mut System, ReadXtcError>;

    /// Read next frame in an xtc file.
    ///
    /// ## Returns
    /// `None` in case the file has been fully read.
    /// `Some(&mut System)` in case the frame has been read successfully.
    /// `Some(ReadXtcError)` in case an error occured while reading.
    fn next(&mut self) -> Option<Self::Item> {
        unsafe {
            let n_atoms = (*self.system).get_n_atoms();

            let atoms = (*self.system).get_atoms_as_ref_mut() as *mut Vec<Atom>;

            let mut step: c_int = 0;
            let mut time: c_float = 0.0;
            let mut boxvector: [[c_float; 3usize]; 3usize] = [[0.0; 3]; 3];
            let mut precision: c_float = 0.0;
            let mut coordinates = vec![[0.0, 0.0, 0.0]; n_atoms as usize];

            // read xtc frame
            let return_code = c_xdrfile::read_xtc(
                self.xtc.handle,
                n_atoms as c_int,
                &mut step,
                &mut time,
                &mut boxvector,
                coordinates.as_mut_ptr(),
                &mut precision,
            );

            match return_code {
                // reading successful and there is more to read
                0 => (),
                // file is completely read
                11 => return None,
                // error occured
                _ => {
                    return Some(Err(ReadXtcError::FrameNotFound));
                }
            }

            for (i, atom) in (*atoms).iter_mut().enumerate() {
                atom.set_position(&Vector3D::from(*coordinates.get_unchecked(i)));
            }

            // update the system
            (*self.system).set_simulation_step(step as u64);
            (*self.system).set_simulation_time(time);
            (*self.system).set_box(matrix2simbox(boxvector));
            (*self.system).set_precision(precision as u64);

            Some(Ok(&mut *self.system))
        }
    }
}

impl System {
    /// Create an `XtcReader` structure which is an iterator over an xtc file.
    ///
    /// ## Returns
    /// `XtcReader` if the xtc file exists and matches the gro file.
    /// Else returns `ReadXtcError`.
    ///
    /// ## Example
    /// Iterating through an xtc trajectory and calculating
    /// and printing the current center of geometry of the system.
    /// ```no_run
    /// use groan_rs::System;
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
    /// use groan_rs::System;
    /// use groan_rs::errors::ReadXtcError;
    ///
    /// fn example_fn() -> Result<(), ReadXtcError> {
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
    /// ## Warning
    /// Only orthogonal simulation boxes are currently supported!
    ///
    /// ## Notes
    /// - The function checks whether the number of atoms in the system corresponds to the number of atoms in the xtc file.
    /// - The `System` structure is modified while iterating through the xtc file.
    pub fn xtc_iter(&mut self, filename: impl AsRef<Path>) -> Result<XtcReader, ReadXtcError> {
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
pub struct XtcWriter {
    system: *const System,
    xtc: Xdrfile,
}

impl XtcWriter {
    /// Open a new xtc file for writing.
    ///
    /// ## Returns
    /// An instance of `XtcWriter` structure or `WriteXtcError` in case the file can't be created.
    ///
    /// ## Example
    /// Create a new xtc file for writing and associate a system with it.
    /// ```no_run
    /// use groan_rs::{System, XtcWriter};
    ///
    /// let system = System::from_file("system.gro").unwrap();
    ///
    /// let mut writer = XtcWriter::new(&system, "output.xtc").unwrap();
    /// ```
    pub fn new(system: &System, filename: impl AsRef<Path>) -> Result<XtcWriter, WriteXtcError> {
        // create the xtc file and save the handle to it
        let xtc = match Xdrfile::open_xtc(filename.as_ref(), OpenMode::Write) {
            Ok(x) => x,
            Err(XtcError::FileNotFound(x)) => return Err(WriteXtcError::CouldNotCreate(x)),
            Err(XtcError::InvalidPath(x)) => return Err(WriteXtcError::InvalidPath(x)),
        };

        Ok(XtcWriter { system, xtc })
    }

    /// Write the current state of the system into an open xtc file.
    ///
    /// ## Example
    /// Reading and writing an xtc file.
    /// ```no_run
    /// use groan_rs::{System, XtcWriter};
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
    ///         writer.write_frame(None)?;
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    pub fn write_frame(&mut self, precision: Option<u64>) -> Result<(), WriteXtcError> {
        unsafe {
            let n_atoms = (*self.system).get_n_atoms();
            let real_precision = precision.unwrap_or((*self.system).get_precision());

            // prepare coordinate matrix
            let mut coordinates = vec![[0.0, 0.0, 0.0]; n_atoms as usize];
            for (i, atom) in (*self.system).atoms_iter().enumerate() {
                let pos = atom.get_position();
                coordinates[i] = [pos.x, pos.y, pos.z];
            }

            // write the xtc frame
            let return_code = c_xdrfile::write_xtc(
                self.xtc.handle,
                n_atoms as c_int,
                (*self.system).get_simulation_step() as i32,
                (*self.system).get_simulation_time(),
                &mut simbox2matrix((*self.system).get_box_as_ref()),
                coordinates.as_mut_ptr(),
                real_precision as f32,
            );

            if return_code != 0 {
                return Err(WriteXtcError::CouldNotWrite);
            }
        }

        Ok(())
    }
}

/******************************/
/*     PRIVATE FUNCTIONS      */
/******************************/

#[derive(Debug, PartialEq, Clone, Copy)]
enum OpenMode {
    Read,
    Write,
}

#[derive(Debug)]
pub struct Xdrfile {
    handle: *mut c_xdrfile::XDRFILE,
}

impl Xdrfile {
    /// Check that the number of atoms in an unopened xtc file matches the expected number.
    fn check_xtc(filename: impl AsRef<Path>, n_atoms: usize) -> Result<bool, ReadXtcError> {
        unsafe {
            let c_path = match path2cstring(filename.as_ref()) {
                Ok(x) => x,
                Err(_) => return Err(ReadXtcError::InvalidPath(Box::from(filename.as_ref()))),
            };

            let mut xtc_atoms: c_int = 0;

            if c_xdrfile::read_xtc_natoms(c_path.as_ptr(), &mut xtc_atoms) != 0 {
                // reading the file failed
                return Err(ReadXtcError::FileNotFound(Box::from(filename.as_ref())));
            }

            // if reading was successful
            if xtc_atoms == n_atoms as c_int {
                Ok(true)
            } else {
                Ok(false)
            }
        }
    }

    /// Open an xtc returning a handle to the file.
    fn open_xtc(filename: impl AsRef<Path>, mode: OpenMode) -> Result<Self, XtcError> {
        unsafe {
            let c_path = match path2cstring(filename.as_ref()) {
                Ok(x) => x,
                Err(_) => return Err(XtcError::InvalidPath(Box::from(filename.as_ref()))),
            };

            let handle = c_xdrfile::xdrfile_open(c_path.as_ptr(), mode2cstring(mode).as_ptr());

            if !handle.is_null() {
                Ok(Xdrfile { handle })
            } else {
                Err(XtcError::FileNotFound(Box::from(filename.as_ref())))
            }
        }
    }
}

impl Drop for Xdrfile {
    /// Close the file as Xdrfile gets dropped.
    fn drop(&mut self) {
        unsafe {
            c_xdrfile::xdrfile_close(self.handle);
        }
    }
}

/// Convert Rust path to null-terminated C string.
fn path2cstring(path: impl AsRef<Path>) -> Result<CString, NulError> {
    CString::new(
        path.as_ref()
            .to_str()
            .expect("Internal Error. Could not convert path to string.")
            .as_bytes(),
    )
}

/// Convert Xdrfile OpenMode to C string.
fn mode2cstring(mode: OpenMode) -> CString {
    match mode {
        OpenMode::Read => std::ffi::CStr::from_bytes_with_nul(b"r\0")
            .unwrap()
            .to_owned(),
        OpenMode::Write => std::ffi::CStr::from_bytes_with_nul(b"w\0")
            .unwrap()
            .to_owned(),
    }
}

/// Convert C box matrix from an xtc file to SimBox.
///
/// ## Warning
/// Currently only works with orthogonal simulation boxes.
fn matrix2simbox(matrix: [[c_float; 3usize]; 3usize]) -> SimBox {
    [matrix[0][0], matrix[1][1], matrix[2][2]].into()
}

/// Convert SimBox to C box matrix for an xtc file.
///
/// ## Warning
/// Currently only works with orthogonal simulation boxes.
fn simbox2matrix(simbox: &SimBox) -> [[c_float; 3usize]; 3usize] {
    [
        [simbox.x, 0.0, 0.0],
        [0.0, simbox.y, 0.0],
        [0.0, 0.0, simbox.z],
    ]
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
    fn read_xtc_unmatching() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();

        match system.xtc_iter("test_files/short_trajectory.xtc") {
            Err(ReadXtcError::AtomsNumberMismatch(_)) => (),
            _ => panic!("XTC file should not be valid."),
        }
    }

    #[test]
    fn read_xtc_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.xtc_iter("test_files/nonexistent.xtc") {
            Err(ReadXtcError::FileNotFound(_)) => (),
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
            writer.write_frame(None).unwrap();
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
            Err(WriteXtcError::CouldNotCreate(_)) => (),
            _ => panic!("Output XTC file should not have been created."),
        }
    }
}
