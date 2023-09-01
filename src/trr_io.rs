// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of functions for reading and writing trr files.

use std::marker::PhantomData;
use std::os::raw::{c_float, c_int};
use std::path::Path;

use crate::atom::Atom;
use crate::errors::{ReadXdrError, WriteXdrError, XdrError};
use crate::system::System;
use crate::vector3d::Vector3D;
use crate::xdrfile::{self, OpenMode, XdrFile, XdrReader, XdrWriter};

/**************************/
/*      READING TRR       */
/**************************/

/// Iterator over a trr file.
pub struct TrrReader<'a> {
    system: *mut System,
    trr: XdrFile,
    phantom: PhantomData<&'a mut System>,
}

impl<'a> XdrReader<'a> for TrrReader<'a> {
    /// Create an iterator over a trr file.
    fn new(system: &'a mut System, filename: impl AsRef<Path>) -> Result<Self, ReadXdrError> {
        let n_atoms = system.get_n_atoms();

        // sanity check the number of atoms
        match XdrFile::check_trr(filename.as_ref(), n_atoms) {
            Err(e) => return Err(e),
            Ok(false) => {
                return Err(ReadXdrError::AtomsNumberMismatch(Box::from(
                    filename.as_ref(),
                )))
            }
            Ok(true) => (),
        };

        // open the trr file and save the handle to it
        let trr = match XdrFile::open_xdr(filename.as_ref(), OpenMode::Read) {
            Ok(x) => x,
            Err(XdrError::FileNotFound(x)) => return Err(ReadXdrError::FileNotFound(x)),
            Err(XdrError::InvalidPath(x)) => return Err(ReadXdrError::InvalidPath(x)),
        };

        Ok(TrrReader {
            system,
            trr,
            phantom: PhantomData,
        })
    }
}

impl<'a> Iterator for TrrReader<'a> {
    type Item = Result<&'a mut System, ReadXdrError>;

    /// Read next frame in a trr file.
    ///
    /// ## Returns
    /// `None` in case the file has been fully read.
    /// `Some(&mut System)` in case the frame has been read successfully.
    /// `Some(ReadXdrError)` in case an error occured while reading.
    fn next(&mut self) -> Option<Self::Item> {
        unsafe {
            let n_atoms = (*self.system).get_n_atoms();

            let atoms = (*self.system).get_atoms_as_ref_mut() as *mut Vec<Atom>;

            let mut step: c_int = 0;
            let mut time: c_float = 0.0;
            let mut boxvector: [[c_float; 3usize]; 3usize] = [[0.0; 3]; 3];
            let mut lambda: c_float = 0.0;
            let mut coordinates = vec![[0.0, 0.0, 0.0]; n_atoms as usize];
            let mut velocities = vec![[0.0, 0.0, 0.0]; n_atoms as usize];
            let mut forces = vec![[0.0, 0.0, 0.0]; n_atoms as usize];

            // read xtc frame
            let return_code = xdrfile::read_trr(
                self.trr.handle,
                n_atoms as c_int,
                &mut step,
                &mut time,
                &mut lambda,
                &mut boxvector,
                coordinates.as_mut_ptr(),
                velocities.as_mut_ptr(),
                forces.as_mut_ptr(),
            );

            match return_code {
                // reading successful and there is more to read
                0 => (),
                // file is completely read
                4 | 11 => return None,
                // error occured
                _ => return Some(Err(ReadXdrError::FrameNotFound)),
            }

            for (i, atom) in (*atoms).iter_mut().enumerate() {
                atom.set_position(&Vector3D::from(*coordinates.get_unchecked(i)));
                atom.set_velocity(&Vector3D::from(*velocities.get_unchecked(i)));
                atom.set_force(&Vector3D::from(*forces.get_unchecked(i)));
            }

            // update the system
            (*self.system).set_simulation_step(step as u64);
            (*self.system).set_simulation_time(time);
            (*self.system).set_box(xdrfile::matrix2simbox(boxvector));
            (*self.system).set_lambda(lambda);

            Some(Ok(&mut *self.system))
        }
    }
}

impl System {
    /// Create a `TrrReader` structure which is an iterator over a trr file.
    ///
    /// ## Returns
    /// `TrrReader` if the trr file exists and matches the structure file.
    /// Else returns `ReadXdrError`.
    ///
    /// ## Example
    /// Iterating through a trr trajectory and calculating
    /// and printing the current center of geometry of the system.
    ///
    /// ```no_run
    /// use groan_rs::prelude::*;
    /// use groan_rs::errors::ReadXdrError;
    ///
    /// fn example_fn() -> Result<(), ReadXdrError> {
    ///     // load system from file
    ///     let mut system = System::from_file("system.gro").unwrap();
    ///
    ///     // loop through the trajectory
    ///     for raw_frame in system.trr_iter("trajectory.trr")? {
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
    /// - The function checks whether the number of atoms in the system corresponds to the number of atoms in the trr file.
    /// - The `System` structure is modified while iterating through the trr file.
    /// - The trr file does not need to have positions, velocities, and forces provided for each frame.
    /// In case any of the properties is missing, it is set to 0 for all particles.
    pub fn trr_iter(&mut self, filename: impl AsRef<Path>) -> Result<TrrReader, ReadXdrError> {
        TrrReader::new(self, filename)
    }
}

/**************************/
/*       WRITING TRR      */
/**************************/

/// Structure for writing trr files.
/// Each `TrrWriter` instance is tightly coupled with a corresponding `System` structure.
/// If you make updates to the `System` structure, such as during iteration with `System::trr_iter()`,
/// and subsequently write a TRR frame using `TrrWriter::write_frame()`, the modifications
/// made to the `System` will be reflected in the written frame.
///
/// `TrrWriter` implements the `XdrWriter` trait.
pub struct TrrWriter {
    system: *const System,
    trr: XdrFile,
}

impl XdrWriter for TrrWriter {
    /// Open a new trr file for writing.
    ///
    /// ## Returns
    /// An instance of `TrrWriter` structure or `WriteTrrError` in case the file can't be created.
    ///
    /// ## Example
    /// Create a new trr file for writing and associate a system with it.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let system = System::from_file("system.gro").unwrap();
    ///
    /// let mut writer = TrrWriter::new(&system, "output.trr").unwrap();
    /// ```
    fn new(system: &System, filename: impl AsRef<Path>) -> Result<TrrWriter, WriteXdrError> {
        // create the trr file and save the handle to it
        let trr = match XdrFile::open_xdr(filename.as_ref(), OpenMode::Write) {
            Ok(x) => x,
            Err(XdrError::FileNotFound(x)) => return Err(WriteXdrError::CouldNotCreate(x)),
            Err(XdrError::InvalidPath(x)) => return Err(WriteXdrError::InvalidPath(x)),
        };

        Ok(TrrWriter { system, trr })
    }

    /// Write the current state of the system into an open trr file.
    ///
    /// ## Example
    /// Reading and writing a trr file.
    /// ```no_run
    /// use groan_rs::prelude::*;
    /// use std::error::Error;
    ///
    /// fn example_fn() -> Result<(), Box<dyn Error>> {
    ///     // load system from file
    ///     let mut system = System::from_file("system.gro")?;
    ///
    ///     // create a trr file for writing and associate it with the system
    ///     let mut writer = TrrWriter::new(&system, "output.trr")?;
    ///
    ///     // loop through the trajectory
    ///     for raw_frame in system.trr_iter("trajectory.trr")? {
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
    /// - While Gromacs supports writing trr files containing only some of the particle properties (e.g. just velocities),
    /// `groan_rs` will always write full trr file with all the properties, even if they are zeroed.
    fn write_frame(&mut self) -> Result<(), WriteXdrError> {
        unsafe {
            let n_atoms = (*self.system).get_n_atoms();

            // prepare coordinate, velocity and forces matrix
            let mut coordinates = vec![[0.0, 0.0, 0.0]; n_atoms as usize];
            let mut velocities = vec![[0.0, 0.0, 0.0]; n_atoms as usize];
            let mut forces = vec![[0.0, 0.0, 0.0]; n_atoms as usize];

            for (i, atom) in (*self.system).atoms_iter().enumerate() {
                let pos = atom.get_position();
                coordinates[i] = [pos.x, pos.y, pos.z];

                let vel = atom.get_velocity();
                velocities[i] = [vel.x, vel.y, vel.z];

                let force = atom.get_force();
                forces[i] = [force.x, force.y, force.z];
            }

            // write the xtc frame
            let return_code = xdrfile::write_trr(
                self.trr.handle,
                n_atoms as c_int,
                (*self.system).get_simulation_step() as i32,
                (*self.system).get_simulation_time(),
                (*self.system).get_lambda(),
                &mut xdrfile::simbox2matrix((*self.system).get_box_as_ref()),
                coordinates.as_mut_ptr(),
                velocities.as_mut_ptr(),
                forces.as_mut_ptr(),
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
    /// Check that the number of atoms in an unopened trr file matches the expected number.
    fn check_trr(filename: impl AsRef<Path>, n_atoms: usize) -> Result<bool, ReadXdrError> {
        unsafe {
            let c_path = match xdrfile::path2cstring(filename.as_ref()) {
                Ok(x) => x,
                Err(_) => return Err(ReadXdrError::InvalidPath(Box::from(filename.as_ref()))),
            };

            let mut trr_atoms: c_int = 0;

            if xdrfile::read_trr_natoms(c_path.as_ptr(), &mut trr_atoms) != 0 {
                // reading the file failed
                return Err(ReadXdrError::FileNotFound(Box::from(filename.as_ref())));
            }

            // if reading was successful
            if trr_atoms == n_atoms as c_int {
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
    use tempfile::NamedTempFile;

    #[test]
    fn read_trr() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        // first frame
        system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .next();

        assert_eq!(system.get_simulation_step(), 0);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 0.0);
        assert_approx_eq!(f32, system.get_box_as_ref().x, 13.01331);
        assert_approx_eq!(f32, system.get_box_as_ref().y, 13.01331);
        assert_approx_eq!(f32, system.get_box_as_ref().z, 11.25347);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_approx_eq!(f32, atom1.get_position().x, 9.497);
        assert_approx_eq!(f32, atom1.get_position().y, 1.989);
        assert_approx_eq!(f32, atom1.get_position().z, 7.498);

        assert_approx_eq!(f32, atom1.get_velocity().x, -0.0683);
        assert_approx_eq!(f32, atom1.get_velocity().y, 0.1133);
        assert_approx_eq!(f32, atom1.get_velocity().z, 0.0005);

        assert_approx_eq!(f32, atom1.get_force().x, -6.2916107);
        assert_approx_eq!(f32, atom1.get_force().y, -276.57983);
        assert_approx_eq!(f32, atom1.get_force().z, -306.23727);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_approx_eq!(f32, atom2.get_position().x, 8.829);
        assert_approx_eq!(f32, atom2.get_position().y, 11.186);
        assert_approx_eq!(f32, atom2.get_position().z, 2.075);

        assert_approx_eq!(f32, atom2.get_velocity().x, 0.0712);
        assert_approx_eq!(f32, atom2.get_velocity().y, 0.2294);
        assert_approx_eq!(f32, atom2.get_velocity().z, -0.1673);

        assert_approx_eq!(f32, atom2.get_force().x, -21.009035);
        assert_approx_eq!(f32, atom2.get_force().y, -6.7285156);
        assert_approx_eq!(f32, atom2.get_force().z, -68.827545);

        // second frame
        system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .nth(1);

        assert_eq!(system.get_simulation_step(), 6000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 120.0);
        assert_approx_eq!(f32, system.get_box_as_ref().x, 13.024242);
        assert_approx_eq!(f32, system.get_box_as_ref().y, 13.024242);
        assert_approx_eq!(f32, system.get_box_as_ref().z, 11.242146);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_approx_eq!(f32, atom1.get_position().x, 0.0);
        assert_approx_eq!(f32, atom1.get_position().y, 0.0);
        assert_approx_eq!(f32, atom1.get_position().z, 0.0);

        assert_approx_eq!(f32, atom1.get_velocity().x, 0.22166125);
        assert_approx_eq!(f32, atom1.get_velocity().y, 0.22522248);
        assert_approx_eq!(f32, atom1.get_velocity().z, 0.19859326);

        assert_approx_eq!(f32, atom1.get_force().x, 0.0);
        assert_approx_eq!(f32, atom1.get_force().y, 0.0);
        assert_approx_eq!(f32, atom1.get_force().z, 0.0);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_approx_eq!(f32, atom2.get_position().x, 0.0);
        assert_approx_eq!(f32, atom2.get_position().y, 0.0);
        assert_approx_eq!(f32, atom2.get_position().z, 0.0);

        assert_approx_eq!(f32, atom2.get_velocity().x, 0.22474734);
        assert_approx_eq!(f32, atom2.get_velocity().y, -0.1732943);
        assert_approx_eq!(f32, atom2.get_velocity().z, -0.1461453);

        assert_approx_eq!(f32, atom2.get_force().x, 0.0);
        assert_approx_eq!(f32, atom2.get_force().y, 0.0);
        assert_approx_eq!(f32, atom2.get_force().z, 0.0);

        // third frame
        system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .nth(2);

        assert_eq!(system.get_simulation_step(), 8000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 160.0);
        assert_approx_eq!(f32, system.get_box_as_ref().x, 13.076236);
        assert_approx_eq!(f32, system.get_box_as_ref().y, 13.076236);
        assert_approx_eq!(f32, system.get_box_as_ref().z, 11.13604);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_approx_eq!(f32, atom1.get_position().x, 0.0);
        assert_approx_eq!(f32, atom1.get_position().y, 0.0);
        assert_approx_eq!(f32, atom1.get_position().z, 0.0);

        assert_approx_eq!(f32, atom1.get_velocity().x, 0.0);
        assert_approx_eq!(f32, atom1.get_velocity().y, 0.0);
        assert_approx_eq!(f32, atom1.get_velocity().z, 0.0);

        assert_approx_eq!(f32, atom1.get_force().x, -167.09401);
        assert_approx_eq!(f32, atom1.get_force().y, -214.71092);
        assert_approx_eq!(f32, atom1.get_force().z, -78.804085);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_approx_eq!(f32, atom2.get_position().x, 0.0);
        assert_approx_eq!(f32, atom2.get_position().y, 0.0);
        assert_approx_eq!(f32, atom2.get_position().z, 0.0);

        assert_approx_eq!(f32, atom2.get_velocity().x, 0.0);
        assert_approx_eq!(f32, atom2.get_velocity().y, 0.0);
        assert_approx_eq!(f32, atom2.get_velocity().z, 0.0);

        assert_approx_eq!(f32, atom2.get_force().x, 230.31451);
        assert_approx_eq!(f32, atom2.get_force().y, -0.87537766);
        assert_approx_eq!(f32, atom2.get_force().z, 72.7905);

        // fourth frame
        system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .nth(3);

        assert_eq!(system.get_simulation_step(), 12000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 240.0);
        assert_approx_eq!(f32, system.get_box_as_ref().x, 13.083817);
        assert_approx_eq!(f32, system.get_box_as_ref().y, 13.083817);
        assert_approx_eq!(f32, system.get_box_as_ref().z, 11.159238);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_approx_eq!(f32, atom1.get_position().x, 9.498894);
        assert_approx_eq!(f32, atom1.get_position().y, 1.8789341);
        assert_approx_eq!(f32, atom1.get_position().z, 7.577659);

        assert_approx_eq!(f32, atom1.get_velocity().x, 0.0472764);
        assert_approx_eq!(f32, atom1.get_velocity().y, 0.003011168);
        assert_approx_eq!(f32, atom1.get_velocity().z, 0.10009501);

        assert_approx_eq!(f32, atom1.get_force().x, 0.0);
        assert_approx_eq!(f32, atom1.get_force().y, 0.0);
        assert_approx_eq!(f32, atom1.get_force().z, 0.0);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_approx_eq!(f32, atom2.get_position().x, 8.397229);
        assert_approx_eq!(f32, atom2.get_position().y, 10.933028);
        assert_approx_eq!(f32, atom2.get_position().z, 2.1274538);

        assert_approx_eq!(f32, atom2.get_velocity().x, 0.39095137);
        assert_approx_eq!(f32, atom2.get_velocity().y, -0.6620998);
        assert_approx_eq!(f32, atom2.get_velocity().z, -0.33029458);

        assert_approx_eq!(f32, atom2.get_force().x, 0.0);
        assert_approx_eq!(f32, atom2.get_force().y, 0.0);
        assert_approx_eq!(f32, atom2.get_force().z, 0.0);

        // last frame
        system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .last();

        assert_eq!(system.get_simulation_step(), 32000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 640.0);
        assert_approx_eq!(f32, system.get_box_as_ref().x, 12.965868);
        assert_approx_eq!(f32, system.get_box_as_ref().y, 12.965868);
        assert_approx_eq!(f32, system.get_box_as_ref().z, 11.348931);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_approx_eq!(f32, atom1.get_position().x, 0.0);
        assert_approx_eq!(f32, atom1.get_position().y, 0.0);
        assert_approx_eq!(f32, atom1.get_position().z, 0.0);

        assert_approx_eq!(f32, atom1.get_velocity().x, 0.0);
        assert_approx_eq!(f32, atom1.get_velocity().y, 0.0);
        assert_approx_eq!(f32, atom1.get_velocity().z, 0.0);

        assert_approx_eq!(f32, atom1.get_force().x, 133.31625);
        assert_approx_eq!(f32, atom1.get_force().y, 66.783325);
        assert_approx_eq!(f32, atom1.get_force().z, 181.96724);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_approx_eq!(f32, atom2.get_position().x, 0.0);
        assert_approx_eq!(f32, atom2.get_position().y, 0.0);
        assert_approx_eq!(f32, atom2.get_position().z, 0.0);

        assert_approx_eq!(f32, atom2.get_velocity().x, 0.0);
        assert_approx_eq!(f32, atom2.get_velocity().y, 0.0);
        assert_approx_eq!(f32, atom2.get_velocity().z, 0.0);

        assert_approx_eq!(f32, atom2.get_force().x, -4.2425976);
        assert_approx_eq!(f32, atom2.get_force().y, 182.99162);
        assert_approx_eq!(f32, atom2.get_force().z, -12.333496);
    }

    #[test]
    fn read_trr_unmatching() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();

        match system.trr_iter("test_files/short_trajectory.trr") {
            Err(ReadXdrError::AtomsNumberMismatch(_)) => (),
            _ => panic!("TRR file should not be valid."),
        }
    }

    #[test]
    fn read_trr_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.trr_iter("test_files/nonexistent.trr") {
            Err(ReadXdrError::FileNotFound(_)) => (),
            _ => panic!("TRR file should not exist."),
        }
    }

    #[test]
    fn write_trr() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        let mut writer = TrrWriter::new(&system, path_to_output).unwrap();

        for _ in system.trr_iter("test_files/short_trajectory.trr").unwrap() {
            writer.write_frame().unwrap();
        }

        // we must close the file, otherwise metadata do not get updated
        drop(writer);

        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        for (raw_frame1, raw_frame2) in system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .zip(system2.trr_iter(path_to_output).unwrap())
        {
            let frame1 = raw_frame1.unwrap();
            let frame2 = raw_frame2.unwrap();

            assert_eq!(frame1.get_simulation_step(), frame2.get_simulation_step());
            assert_approx_eq!(
                f32,
                frame1.get_simulation_time(),
                frame2.get_simulation_time()
            );

            let box1 = frame1.get_box_as_ref();
            let box2 = frame2.get_box_as_ref();
            assert_approx_eq!(f32, box1.x, box2.x);
            assert_approx_eq!(f32, box1.y, box2.y);
            assert_approx_eq!(f32, box1.z, box2.z);

            for (a1, a2) in frame1
                .get_atoms_as_ref()
                .iter()
                .zip(frame2.get_atoms_as_ref().iter())
            {
                assert_eq!(a1.get_atom_name(), a2.get_atom_name());
                assert_eq!(a1.get_atom_number(), a2.get_atom_number());

                let pos1 = a1.get_position();
                let pos2 = a2.get_position();
                assert_approx_eq!(f32, pos1.x, pos2.x);
                assert_approx_eq!(f32, pos1.y, pos2.y);
                assert_approx_eq!(f32, pos1.z, pos2.z);

                let vel1 = a1.get_velocity();
                let vel2 = a2.get_velocity();
                assert_approx_eq!(f32, vel1.x, vel2.x);
                assert_approx_eq!(f32, vel1.y, vel2.y);
                assert_approx_eq!(f32, vel1.z, vel2.z);

                let force1 = a1.get_force();
                let force2 = a2.get_force();
                assert_approx_eq!(f32, force1.x, force2.x);
                assert_approx_eq!(f32, force1.y, force2.y);
                assert_approx_eq!(f32, force1.z, force2.z);
            }
        }
    }

    #[test]
    fn write_invalid_path() {
        let system = System::from_file("test_files/example.gro").unwrap();

        match TrrWriter::new(&system, "test_files/nonexistent/output.trr") {
            Err(WriteXdrError::CouldNotCreate(_)) => (),
            _ => panic!("Output TRR file should not have been created."),
        }
    }
}
