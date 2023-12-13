// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of functions for reading and writing trr files.

use std::marker::PhantomData;
use std::os::raw::{c_float, c_int};
use std::path::Path;

use crate::errors::{ReadTrajError, TrajError, WriteTrajError};
use crate::io::traj_io::{
    FrameData, FrameDataTime, TrajGroupWrite, TrajRangeRead, TrajRead, TrajReader, TrajStepRead,
    TrajWrite,
};
use crate::io::xdrfile::{self, CXdrFile, OpenMode, XdrFile};
use crate::structures::iterators::AtomIterator;
use crate::structures::{group::Group, simbox::SimBox, vector3d::Vector3D};
use crate::system::general::System;

/**************************/
/*      READING TRR       */
/**************************/

/// Structure containing data read from a trr file.
#[derive(Debug)]
pub struct TrrFrameData {
    step: c_int,
    time: c_float,
    simbox: SimBox,
    lambda: c_float,
    coordinates: Vec<[f32; 3]>,
    velocities: Vec<[f32; 3]>,
    forces: Vec<[f32; 3]>,
}

impl FrameData for TrrFrameData {
    type TrajFile = CXdrFile;

    /// Read `TrrFrameData` from a trr frame.
    ///
    /// ## Returns
    /// - `None` if the file has been read completely
    /// - `Some(TrrFrameData)` if the frame has been successfully read.
    /// - `Some(ReadTrajError)` if the frame could not be read.
    fn from_frame(
        xdrfile: &mut Self::TrajFile,
        system: &System,
    ) -> Option<Result<Self, ReadTrajError>> {
        let mut step: c_int = 0;
        let mut time: c_float = 0.0;
        let mut boxvector: [[c_float; 3usize]; 3usize] = [[0.0; 3]; 3];
        let mut lambda: c_float = 0.0;
        let mut coordinates = vec![[0.0, 0.0, 0.0]; system.get_n_atoms()];
        let mut velocities = vec![[0.0, 0.0, 0.0]; system.get_n_atoms()];
        let mut forces = vec![[0.0, 0.0, 0.0]; system.get_n_atoms()];

        unsafe {
            let return_code = xdrfile::read_trr(
                xdrfile,
                system.get_n_atoms() as c_int,
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
                0 => {
                    let simbox = match xdrfile::matrix2simbox(boxvector) {
                        Ok(x) => x,
                        Err(e) => return Some(Err(e)),
                    };

                    Some(Ok(TrrFrameData {
                        step,
                        time,
                        simbox,
                        lambda,
                        coordinates,
                        velocities,
                        forces,
                    }))
                }
                // file is completely read
                4 | 11 => None,
                // error occured
                _ => Some(Err(ReadTrajError::FrameNotFound)),
            }
        }
    }

    /// Update the `System` structure based on data from `TrrFrameData`.
    fn update_system(self, system: &mut System) {
        unsafe {
            for (i, atom) in system.get_atoms_as_ref_mut().iter_mut().enumerate() {
                let pos = Vector3D::from(*self.coordinates.get_unchecked(i));
                let vel = Vector3D::from(*self.velocities.get_unchecked(i));
                let force = Vector3D::from(*self.forces.get_unchecked(i));

                if pos.is_zero() {
                    atom.reset_position();
                } else {
                    atom.set_position(pos);
                }

                if vel.is_zero() {
                    atom.reset_velocity();
                } else {
                    atom.set_velocity(vel);
                }

                if force.is_zero() {
                    atom.reset_force();
                } else {
                    atom.set_force(force);
                }
            }

            // update the system
            system.set_simulation_step(self.step as u64);
            system.set_simulation_time(self.time);
            system.set_box(self.simbox);
            system.set_lambda(self.lambda);
        }
    }
}

impl FrameDataTime for TrrFrameData {
    fn get_time(&self) -> f32 {
        self.time
    }
}

/// Iterator over a trr file.
pub struct TrrReader<'a> {
    system: *mut System,
    trr: XdrFile,
    phantom: PhantomData<&'a mut System>,
}

impl<'a> TrajRead<'a> for TrrReader<'a> {
    type FrameData = TrrFrameData;

    /// Create an iterator over a trr file.
    fn new(system: &'a mut System, filename: impl AsRef<Path>) -> Result<TrrReader, ReadTrajError> {
        let n_atoms = system.get_n_atoms();

        // sanity check the number of atoms
        match XdrFile::check_trr(filename.as_ref(), n_atoms) {
            Err(e) => return Err(e),
            Ok(false) => {
                return Err(ReadTrajError::AtomsNumberMismatch(Box::from(
                    filename.as_ref(),
                )))
            }
            Ok(true) => (),
        };

        // open the trr file and save the handle to it
        let trr = match XdrFile::open_xdr(filename.as_ref(), OpenMode::Read) {
            Ok(x) => x,
            Err(TrajError::FileNotFound(x)) => return Err(ReadTrajError::FileNotFound(x)),
            Err(TrajError::InvalidPath(x)) => return Err(ReadTrajError::InvalidPath(x)),
        };

        let trr_reader = TrrReader {
            system,
            trr,
            phantom: PhantomData,
        };

        Ok(trr_reader)
    }

    fn get_system(&mut self) -> *mut System {
        self.system
    }

    fn get_file_handle(&mut self) -> &mut CXdrFile {
        unsafe { self.trr.handle.as_mut().unwrap() }
    }
}

impl<'a> TrajRangeRead<'a> for TrrReader<'a> {
    fn jump_to_start(&mut self, start_time: f32) -> Result<(), ReadTrajError> {
        unsafe {
            if xdrfile::trr_jump_to_start(self.get_file_handle(), start_time as c_float) != 0 {
                Err(ReadTrajError::StartNotFound(start_time.to_string()))
            } else {
                Ok(())
            }
        }
    }
}

impl<'a> TrajStepRead<'a> for TrrReader<'a> {
    fn skip_frame(&mut self) -> Result<bool, ReadTrajError> {
        unsafe {
            match xdrfile::trr_skip_frame(self.get_file_handle()) {
                0 => Ok(true),
                1 => Err(ReadTrajError::SkipFailed),
                2 => Ok(false),
                number => panic!(
                    "FATAL GROAN ERROR | TrrReader::skip_frame | `xdrfile::trr_skip_frame` returned an unsupported number '{}'.",
                    number
                ),
            }
        }
    }
}

/// ## Methods for reading trr files.
impl System {
    /// Create a `TrrReader` structure which is an iterator over a trr file.
    ///
    /// ## Returns
    /// `TrajReader<TrrReader>` if the trr file exists and matches the structure file.
    /// Else returns `ReadTrajError`.
    ///
    /// ## Examples
    /// Iterating through a trr trajectory and calculating
    /// and printing the current center of geometry of the system.
    ///
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use groan_rs::errors::ReadTrajError;
    /// #
    /// fn example_fn() -> Result<(), ReadTrajError> {
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
    /// You can also iterate over just a part of the trajectory.
    /// Here, only frames in the time range 10-100 ns will be read.
    /// The `with_range` method is very efficient for the trr files
    /// as the particle properties from the trr frames outside of
    /// the range will not be read.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use groan_rs::errors::ReadTrajError;
    /// #
    /// fn example_fn() -> Result<(), ReadTrajError> {
    ///     let mut system = System::from_file("system.gro").unwrap();
    ///
    ///     for raw_frame in system
    ///         .trr_iter("trajectory.trr")?
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
    /// The `with_step` method is very efficient for the trr files as the
    /// particle properties from the skipped over frames will not be read.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use groan_rs::errors::ReadTrajError;
    /// #
    /// fn example_fn() -> Result<(), ReadTrajError> {
    ///     let mut system = System::from_file("system.gro").unwrap();
    ///
    ///     for raw_frame in system
    ///         .trr_iter("trajectory.trr")?
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
    ///         .trr_iter("trajectory.trr")?
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
    /// ## Warning
    /// Only orthogonal simulation boxes are currently supported!
    ///
    /// ## Notes
    /// - The function checks whether the number of atoms in the system corresponds to the number of atoms in the trr file.
    /// - The `System` structure is modified while iterating through the trr file.
    /// - The trr file does not need to have positions, velocities, and forces provided for each frame.
    /// In case any of these properties is missing, it is set to `None` for all atoms.
    pub fn trr_iter(
        &mut self,
        filename: impl AsRef<Path>,
    ) -> Result<TrajReader<TrrReader>, ReadTrajError> {
        Ok(TrajReader::wrap_traj(TrrReader::new(self, filename)?))
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
/// `TrrWriter` implements the `TrajWrite` trait.
pub struct TrrWriter {
    system: *const System,
    trr: XdrFile,
}

impl TrajWrite for TrrWriter {
    /// Open a new trr file for writing.
    ///
    /// ## Returns
    /// An instance of `TrrWriter` structure or `WriteTrajError` in case the file can't be created.
    ///
    /// ## Example
    /// Create a new trr file for writing and associate a system with it.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let system = System::from_file("system.gro").unwrap();
    /// let mut writer = TrrWriter::new(&system, "output.trr").unwrap();
    /// ```
    fn new(system: &System, filename: impl AsRef<Path>) -> Result<TrrWriter, WriteTrajError> {
        // create the trr file and save a handle to it
        let trr = match XdrFile::open_xdr(filename.as_ref(), OpenMode::Write) {
            Ok(x) => x,
            Err(TrajError::FileNotFound(x)) => return Err(WriteTrajError::CouldNotCreate(x)),
            Err(TrajError::InvalidPath(x)) => return Err(WriteTrajError::InvalidPath(x)),
        };

        Ok(TrrWriter { system, trr })
    }

    /// Write the current state of the system into an open trr file.
    ///
    /// ## Returns
    /// - `Ok` if the frame has been successfully written. Otherwise `WriteTrajError`.
    ///
    /// ## Example
    /// Reading and writing a trr file.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use std::error::Error;
    /// #
    /// fn example_fn() -> Result<(), Box<dyn Error + Send + Sync>> {
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
    fn write_frame(&mut self) -> Result<(), WriteTrajError> {
        unsafe {
            let n_atoms = (*self.system).get_n_atoms();

            // prepare coordinate, velocity and forces matrix
            let mut coordinates = vec![[0.0, 0.0, 0.0]; n_atoms];
            let mut velocities = vec![[0.0, 0.0, 0.0]; n_atoms];
            let mut forces = vec![[0.0, 0.0, 0.0]; n_atoms];

            for (i, atom) in (*self.system).atoms_iter().enumerate() {
                if let Some(pos) = atom.get_position() {
                    coordinates[i] = [pos.x, pos.y, pos.z];
                }

                if let Some(vel) = atom.get_velocity() {
                    velocities[i] = [vel.x, vel.y, vel.z];
                }

                if let Some(force) = atom.get_force() {
                    forces[i] = [force.x, force.y, force.z]
                }
            }

            // write the trr frame
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
                return Err(WriteTrajError::CouldNotWrite);
            }
        }

        Ok(())
    }
}

/// Structure for writing groups of atoms into trr files.
/// Each `TrrGroupWriter` is tightly coupled with a corresponding `Group` from `System` structure.
/// If you make updates to the `System` structure, such as during the iteration with `System::trr_iter()`,
/// and subsequently write a TRR frame using `TrrGroupWriter::write_frame()`, the modifications
/// made to the `System` will be reflected in the written frame.
///
/// Note that the purpose of the `TrrGroupWriter` is writing valid trr files with consistent number of atoms.
/// Therefore, the `TrrGroupWriter` always works with the original provided group of atoms.
/// If you change the meaning of this group after constructing `TrrGroupWriter` (by overwriting the group),
/// `TrrGroupWriter` will still use the original group of atoms.
/// If you completely remove the original group, `TrrGroupWriter` will still maintain a working copy of it.
///
/// `TrrGroupWriter` implements the `TrajGroupWrite` trait.
pub struct TrrGroupWriter {
    system: *const System,
    trr: XdrFile,
    /// This is a deep copy of the group from the system. `TrrGroupWriter` must always work, even if the user removes the group from the `System` or overwrites it.
    group: Group,
}

impl TrajGroupWrite for TrrGroupWriter {
    /// Open a new trr file for writing and associate a specific group from a specific system with it.
    ///
    /// ## Returns
    /// An instance of `TrrGroupWriter` structure or `WriteTrajError` in case the file can't be created
    /// or the group does not exist.
    ///
    /// ## Example
    /// Create a new trr file for writing and associate a group with it.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("My Group", "resid 1-4").unwrap();
    ///
    /// let mut writer = TrrGroupWriter::new(&system, "My Group", "output.trr").unwrap();
    /// ```
    fn new(
        system: &System,
        group_name: &str,
        filename: impl AsRef<Path>,
    ) -> Result<TrrGroupWriter, WriteTrajError> {
        // get copy of the group
        let group = match system.get_groups_as_ref().get(group_name) {
            None => return Err(WriteTrajError::GroupNotFound(group_name.to_owned())),
            Some(g) => g.clone(),
        };

        // create the trr file and save a handle to it
        let trr = match XdrFile::open_xdr(filename.as_ref(), OpenMode::Write) {
            Ok(x) => x,
            Err(TrajError::FileNotFound(x)) => return Err(WriteTrajError::CouldNotCreate(x)),
            Err(TrajError::InvalidPath(x)) => return Err(WriteTrajError::InvalidPath(x)),
        };

        Ok(TrrGroupWriter { system, trr, group })
    }

    /// Write the current state of the group into an open trr file.
    ///
    /// ## Returns
    /// - `Ok` if the frame has been successfully written. Otherwise `WriteTrajError`.
    ///
    /// ## Example
    /// Reading a trr file and writing only the atoms corresponding to an ndx group `Protein` into the output trr file.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use std::error::Error;
    /// #
    /// fn example_fn() -> Result<(), Box<dyn Error + Send + Sync>> {
    ///     // load system from file
    ///     let mut system = System::from_file("system.gro")?;
    ///     system.read_ndx("index.ndx")?;
    ///
    ///     // create a trr file for writing and associate it with the group `Protein`
    ///     let mut writer = TrrGroupWriter::new(&system, "Protein", "output.trr")?;
    ///
    ///     // loop through the trajectory
    ///     for raw_frame in system.trr_iter("trajectory.trr")? {
    ///         // check for errors
    ///         let _ = raw_frame?;
    ///
    ///         // write the current state of the group `Protein` into `output.trr`
    ///         writer.write_frame()?;
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    /// ## Notes
    /// - While Gromacs supports writing trr files containing only some of the particle properties (e.g. just velocities),
    /// `groan_rs` will always write full trr file with all the properties, even if they are zeroed.
    fn write_frame(&mut self) -> Result<(), WriteTrajError> {
        unsafe {
            let n_atoms = self.group.get_n_atoms();

            // create an iterator over the atoms of the group
            let iterator = AtomIterator::new(
                (*self.system).get_atoms_as_ref(),
                self.group.get_atoms(),
                (*self.system).get_box_as_ref(),
            );

            // prepare coordinate, velocity and forces matrix
            let mut coordinates = vec![[0.0, 0.0, 0.0]; n_atoms];
            let mut velocities = vec![[0.0, 0.0, 0.0]; n_atoms];
            let mut forces = vec![[0.0, 0.0, 0.0]; n_atoms];

            for (i, atom) in iterator.enumerate() {
                if let Some(pos) = atom.get_position() {
                    coordinates[i] = [pos.x, pos.y, pos.z];
                }

                if let Some(vel) = atom.get_velocity() {
                    velocities[i] = [vel.x, vel.y, vel.z];
                }

                if let Some(force) = atom.get_force() {
                    forces[i] = [force.x, force.y, force.z]
                }
            }

            // write the trr frame
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
                return Err(WriteTrajError::CouldNotWrite);
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
    fn check_trr(filename: impl AsRef<Path>, n_atoms: usize) -> Result<bool, ReadTrajError> {
        unsafe {
            let c_path = match xdrfile::path2cstring(filename.as_ref()) {
                Ok(x) => x,
                Err(_) => return Err(ReadTrajError::InvalidPath(Box::from(filename.as_ref()))),
            };

            let mut trr_atoms: c_int = 0;

            if xdrfile::read_trr_natoms(c_path.as_ptr(), &mut trr_atoms) != 0 {
                // reading the file failed
                return Err(ReadTrajError::FileNotFound(Box::from(filename.as_ref())));
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
    use crate::test_utilities::utilities::compare_atoms;
    use float_cmp::assert_approx_eq;
    use std::fs::File;
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
        let simbox = system.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.01331);
        assert_approx_eq!(f32, simbox.y, 13.01331);
        assert_approx_eq!(f32, simbox.z, 11.25347);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_approx_eq!(f32, atom1.get_position().unwrap().x, 9.497);
        assert_approx_eq!(f32, atom1.get_position().unwrap().y, 1.989);
        assert_approx_eq!(f32, atom1.get_position().unwrap().z, 7.498);

        assert_approx_eq!(f32, atom1.get_velocity().unwrap().x, -0.0683);
        assert_approx_eq!(f32, atom1.get_velocity().unwrap().y, 0.1133);
        assert_approx_eq!(f32, atom1.get_velocity().unwrap().z, 0.0005);

        assert_approx_eq!(f32, atom1.get_force().unwrap().x, -6.2916107);
        assert_approx_eq!(f32, atom1.get_force().unwrap().y, -276.57983);
        assert_approx_eq!(f32, atom1.get_force().unwrap().z, -306.23727);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_approx_eq!(f32, atom2.get_position().unwrap().x, 8.829);
        assert_approx_eq!(f32, atom2.get_position().unwrap().y, 11.186);
        assert_approx_eq!(f32, atom2.get_position().unwrap().z, 2.075);

        assert_approx_eq!(f32, atom2.get_velocity().unwrap().x, 0.0712);
        assert_approx_eq!(f32, atom2.get_velocity().unwrap().y, 0.2294);
        assert_approx_eq!(f32, atom2.get_velocity().unwrap().z, -0.1673);

        assert_approx_eq!(f32, atom2.get_force().unwrap().x, -21.009035);
        assert_approx_eq!(f32, atom2.get_force().unwrap().y, -6.7285156);
        assert_approx_eq!(f32, atom2.get_force().unwrap().z, -68.827545);

        // second frame
        system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .nth(1);

        assert_eq!(system.get_simulation_step(), 6000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 120.0);
        let simbox = system.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.024242);
        assert_approx_eq!(f32, simbox.y, 13.024242);
        assert_approx_eq!(f32, simbox.z, 11.242146);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_eq!(atom1.get_position(), None);

        assert_approx_eq!(f32, atom1.get_velocity().unwrap().x, 0.22166125);
        assert_approx_eq!(f32, atom1.get_velocity().unwrap().y, 0.22522248);
        assert_approx_eq!(f32, atom1.get_velocity().unwrap().z, 0.19859326);

        assert_eq!(atom1.get_force(), None);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_eq!(atom2.get_position(), None);

        assert_approx_eq!(f32, atom2.get_velocity().unwrap().x, 0.22474734);
        assert_approx_eq!(f32, atom2.get_velocity().unwrap().y, -0.1732943);
        assert_approx_eq!(f32, atom2.get_velocity().unwrap().z, -0.1461453);

        assert_eq!(atom2.get_force(), None);

        // third frame
        system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .nth(2);

        assert_eq!(system.get_simulation_step(), 8000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 160.0);
        let simbox = system.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.076236);
        assert_approx_eq!(f32, simbox.y, 13.076236);
        assert_approx_eq!(f32, simbox.z, 11.13604);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_eq!(atom1.get_position(), None);

        assert_eq!(atom1.get_velocity(), None);

        assert_approx_eq!(f32, atom1.get_force().unwrap().x, -167.09401);
        assert_approx_eq!(f32, atom1.get_force().unwrap().y, -214.71092);
        assert_approx_eq!(f32, atom1.get_force().unwrap().z, -78.804085);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_eq!(atom2.get_position(), None);

        assert_eq!(atom2.get_velocity(), None);

        assert_approx_eq!(f32, atom2.get_force().unwrap().x, 230.31451);
        assert_approx_eq!(f32, atom2.get_force().unwrap().y, -0.87537766);
        assert_approx_eq!(f32, atom2.get_force().unwrap().z, 72.7905);

        // fourth frame
        system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .nth(3);

        assert_eq!(system.get_simulation_step(), 12000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 240.0);
        let simbox = system.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.083817);
        assert_approx_eq!(f32, simbox.y, 13.083817);
        assert_approx_eq!(f32, simbox.z, 11.159238);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_approx_eq!(f32, atom1.get_position().unwrap().x, 9.498894);
        assert_approx_eq!(f32, atom1.get_position().unwrap().y, 1.8789341);
        assert_approx_eq!(f32, atom1.get_position().unwrap().z, 7.577659);

        assert_approx_eq!(f32, atom1.get_velocity().unwrap().x, 0.0472764);
        assert_approx_eq!(f32, atom1.get_velocity().unwrap().y, 0.003011168);
        assert_approx_eq!(f32, atom1.get_velocity().unwrap().z, 0.10009501);

        assert_eq!(atom1.get_force(), None);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_approx_eq!(f32, atom2.get_position().unwrap().x, 8.397229);
        assert_approx_eq!(f32, atom2.get_position().unwrap().y, 10.933028);
        assert_approx_eq!(f32, atom2.get_position().unwrap().z, 2.1274538);

        assert_approx_eq!(f32, atom2.get_velocity().unwrap().x, 0.39095137);
        assert_approx_eq!(f32, atom2.get_velocity().unwrap().y, -0.6620998);
        assert_approx_eq!(f32, atom2.get_velocity().unwrap().z, -0.33029458);

        assert_eq!(atom2.get_force(), None);

        // last frame
        system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .last();

        assert_eq!(system.get_simulation_step(), 32000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 640.0);
        let simbox = system.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.x, 12.965868);
        assert_approx_eq!(f32, simbox.y, 12.965868);
        assert_approx_eq!(f32, simbox.z, 11.348931);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_eq!(atom1.get_position(), None);

        assert_eq!(atom1.get_velocity(), None);

        assert_approx_eq!(f32, atom1.get_force().unwrap().x, 133.31625);
        assert_approx_eq!(f32, atom1.get_force().unwrap().y, 66.783325);
        assert_approx_eq!(f32, atom1.get_force().unwrap().z, 181.96724);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_eq!(atom2.get_position(), None);

        assert_eq!(atom2.get_velocity(), None);

        assert_approx_eq!(f32, atom2.get_force().unwrap().x, -4.2425976);
        assert_approx_eq!(f32, atom2.get_force().unwrap().y, 182.99162);
        assert_approx_eq!(f32, atom2.get_force().unwrap().z, -12.333496);
    }

    #[test]
    fn read_trr_double_precision() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        // first frame
        system
            .trr_iter("test_files/short_trajectory_double.trr")
            .unwrap()
            .next();

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(system.get_simulation_step(), 0);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 0.0);
        let simbox = system.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.01331);
        assert_approx_eq!(f32, simbox.y, 13.01331);
        assert_approx_eq!(f32, simbox.z, 11.25347);

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_approx_eq!(f32, atom1.get_position().unwrap().x, 9.497161);
        assert_approx_eq!(f32, atom1.get_position().unwrap().y, 1.9891102);
        assert_approx_eq!(f32, atom1.get_position().unwrap().z, 7.497941);

        assert_approx_eq!(f32, atom1.get_velocity().unwrap().x, -0.06389237);
        assert_approx_eq!(f32, atom1.get_velocity().unwrap().y, 0.054320477);
        assert_approx_eq!(f32, atom1.get_velocity().unwrap().z, 0.008154817);

        assert_approx_eq!(f32, atom1.get_force().unwrap().x, -6.330056);
        assert_approx_eq!(f32, atom1.get_force().unwrap().y, -278.8763);
        assert_approx_eq!(f32, atom1.get_force().unwrap().z, -305.94952);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_approx_eq!(f32, atom2.get_position().unwrap().x, 8.829);
        assert_approx_eq!(f32, atom2.get_position().unwrap().y, 11.186);
        assert_approx_eq!(f32, atom2.get_position().unwrap().z, 2.075);

        assert_approx_eq!(f32, atom2.get_velocity().unwrap().x, 0.16692738);
        assert_approx_eq!(f32, atom2.get_velocity().unwrap().y, 0.1674121);
        assert_approx_eq!(f32, atom2.get_velocity().unwrap().z, -0.27088445);

        assert_approx_eq!(f32, atom2.get_force().unwrap().x, -21.007483);
        assert_approx_eq!(f32, atom2.get_force().unwrap().y, -6.727664);
        assert_approx_eq!(f32, atom2.get_force().unwrap().z, -68.82874);

        // second frame
        system
            .trr_iter("test_files/short_trajectory_double.trr")
            .unwrap()
            .nth(1);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(system.get_simulation_step(), 6000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 120.0);
        let simbox = system.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.00128);
        assert_approx_eq!(f32, simbox.y, 13.00128);
        assert_approx_eq!(f32, simbox.z, 11.277555);

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_eq!(atom1.get_position(), None);

        assert_approx_eq!(f32, atom1.get_velocity().unwrap().x, 0.14590994);
        assert_approx_eq!(f32, atom1.get_velocity().unwrap().y, 0.02281682);
        assert_approx_eq!(f32, atom1.get_velocity().unwrap().z, 0.118289664);

        assert_eq!(atom1.get_force(), None);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_eq!(atom2.get_position(), None);

        assert_approx_eq!(f32, atom2.get_velocity().unwrap().x, -0.10735397);
        assert_approx_eq!(f32, atom2.get_velocity().unwrap().y, -0.024522306);
        assert_approx_eq!(f32, atom2.get_velocity().unwrap().z, 0.37654695);

        assert_eq!(atom2.get_force(), None);

        // third frame
        system
            .trr_iter("test_files/short_trajectory_double.trr")
            .unwrap()
            .nth(2);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(system.get_simulation_step(), 8000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 160.0);
        let simbox = system.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.117418);
        assert_approx_eq!(f32, simbox.y, 13.117418);
        assert_approx_eq!(f32, simbox.z, 11.052416);

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_eq!(atom1.get_position(), None);

        assert_eq!(atom1.get_velocity(), None);

        assert_approx_eq!(f32, atom1.get_force().unwrap().x, -1.2237711);
        assert_approx_eq!(f32, atom1.get_force().unwrap().y, -132.20737);
        assert_approx_eq!(f32, atom1.get_force().unwrap().z, 83.95251);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_eq!(atom2.get_position(), None);

        assert_eq!(atom2.get_velocity(), None);

        assert_approx_eq!(f32, atom2.get_force().unwrap().x, 0.82113415);
        assert_approx_eq!(f32, atom2.get_force().unwrap().y, -31.931189);
        assert_approx_eq!(f32, atom2.get_force().unwrap().z, -17.756308);

        // fourth frame
        system
            .trr_iter("test_files/short_trajectory_double.trr")
            .unwrap()
            .nth(3);

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(system.get_simulation_step(), 12000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 240.0);
        let simbox = system.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.020565);
        assert_approx_eq!(f32, simbox.y, 13.020565);
        assert_approx_eq!(f32, simbox.z, 11.256957);

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_approx_eq!(f32, atom1.get_position().unwrap().x, 8.99423);
        assert_approx_eq!(f32, atom1.get_position().unwrap().y, 1.8623593);
        assert_approx_eq!(f32, atom1.get_position().unwrap().z, 7.1457996);

        assert_approx_eq!(f32, atom1.get_velocity().unwrap().x, -0.17861652);
        assert_approx_eq!(f32, atom1.get_velocity().unwrap().y, 0.092642576);
        assert_approx_eq!(f32, atom1.get_velocity().unwrap().z, 0.0057291547);

        assert_eq!(atom1.get_force(), None);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_approx_eq!(f32, atom2.get_position().unwrap().x, 9.4999075);
        assert_approx_eq!(f32, atom2.get_position().unwrap().y, 11.252099);
        assert_approx_eq!(f32, atom2.get_position().unwrap().z, 1.6288123);

        assert_approx_eq!(f32, atom2.get_velocity().unwrap().x, -0.072873585);
        assert_approx_eq!(f32, atom2.get_velocity().unwrap().y, -0.28077266);
        assert_approx_eq!(f32, atom2.get_velocity().unwrap().z, -0.06289119);

        assert_eq!(atom2.get_force(), None);

        // last frame
        system
            .trr_iter("test_files/short_trajectory_double.trr")
            .unwrap()
            .last();

        let atom1 = &system.get_atoms_as_ref()[0];
        let atom2 = &system.get_atoms_as_ref()[16843];

        assert_eq!(system.get_simulation_step(), 32000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 640.0);
        let simbox = system.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.068602);
        assert_approx_eq!(f32, simbox.y, 13.068602);
        assert_approx_eq!(f32, simbox.z, 11.147263);

        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_atom_name(), "BB");

        assert_eq!(atom1.get_position(), None);

        assert_eq!(atom1.get_velocity(), None);

        assert_approx_eq!(f32, atom1.get_force().unwrap().x, -253.34071);
        assert_approx_eq!(f32, atom1.get_force().unwrap().y, -54.76411);
        assert_approx_eq!(f32, atom1.get_force().unwrap().z, 167.09177);

        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_atom_name(), "CL");

        assert_eq!(atom2.get_position(), None);

        assert_eq!(atom2.get_velocity(), None);

        assert_approx_eq!(f32, atom2.get_force().unwrap().x, -13.865962);
        assert_approx_eq!(f32, atom2.get_force().unwrap().y, -36.480534);
        assert_approx_eq!(f32, atom2.get_force().unwrap().z, -88.47915);
    }

    #[test]
    fn read_trr_iter() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let times = [0, 120, 160, 240, 320, 360, 480, 600, 640];

        for (i, raw) in system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .enumerate()
        {
            let frame = raw.unwrap();

            assert_eq!(frame.get_simulation_time() as usize, times[i]);
        }
    }

    #[test]
    fn read_trr_iter_double_precision() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let times = [0, 120, 160, 240, 320, 360, 480, 600, 640];

        for (i, raw) in system
            .trr_iter("test_files/short_trajectory_double.trr")
            .unwrap()
            .enumerate()
        {
            let frame = raw.unwrap();

            assert_eq!(frame.get_simulation_time() as usize, times[i]);
        }
    }

    #[test]
    fn read_trr_range_full() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        for (raw_frame1, raw_frame2) in system1
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .zip(
                system2
                    .trr_iter("test_files/short_trajectory.trr")
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
    fn read_trr_range() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        let mut i = 0;

        for (raw_frame1, raw_frame2) in system1
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .skip(3)
            .zip(
                system2
                    .trr_iter("test_files/short_trajectory.trr")
                    .unwrap()
                    .with_range(200.0, 600.0)
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

        assert_eq!(i, 5);
    }

    #[test]
    fn read_trr_range_double_precision() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        let mut i = 0;

        for (raw_frame1, raw_frame2) in system1
            .trr_iter("test_files/short_trajectory_double.trr")
            .unwrap()
            .skip(3)
            .zip(
                system2
                    .trr_iter("test_files/short_trajectory_double.trr")
                    .unwrap()
                    .with_range(200.0, 600.0)
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

        assert_eq!(i, 5);
    }

    #[test]
    fn read_trr_range_full_iter() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let times = [0, 120, 160, 240, 320, 360, 480, 600, 640];

        for (i, raw) in system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .with_range(0.0, 1000.0)
            .unwrap()
            .enumerate()
        {
            let frame = raw.unwrap();

            assert_eq!(frame.get_simulation_time() as usize, times[i]);
        }
    }

    #[test]
    fn read_trr_range_iter() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let times = [240, 320, 360, 480, 600];

        for (i, raw) in system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .with_range(200.0, 600.0)
            .unwrap()
            .enumerate()
        {
            let frame = raw.unwrap();

            assert_eq!(frame.get_simulation_time() as usize, times[i]);
        }
    }

    #[test]
    fn read_trr_range_negative() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .with_range(-300.0, 800.0)
        {
            Ok(_) => panic!("Iterator should not have been constructed."),
            Err(ReadTrajError::TimeRangeNegative(_)) => (),
            Err(e) => panic!("Incorrect error type {} returned", e),
        }

        match system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .with_range(300.0, -800.0)
        {
            Ok(_) => panic!("Iterator should not have been constructed."),
            Err(ReadTrajError::TimeRangeNegative(_)) => (),
            Err(e) => panic!("Incorrect error type {} returned", e),
        }
    }

    #[test]
    fn read_trr_range_end_start() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .with_range(800.0, 300.0)
        {
            Ok(_) => panic!("Iterator should not have been constructed."),
            Err(ReadTrajError::InvalidTimeRange(_, _)) => (),
            Err(e) => panic!("Incorrect error type {} returned", e),
        }
    }

    #[test]
    fn read_trr_range_start_not_found() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .with_range(12000.0, 20000.0)
        {
            Ok(_) => panic!("Iterator should not have been constructed."),
            Err(ReadTrajError::StartNotFound(_)) => (),
            Err(e) => panic!("Incorrect error type {} returned", e),
        }
    }

    #[test]
    fn read_trr_step_1() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        for (raw1, raw2) in system1
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .with_step(1)
            .unwrap()
            .zip(system2.trr_iter("test_files/short_trajectory.trr").unwrap())
        {
            let frame1 = raw1.unwrap();
            let frame2 = raw2.unwrap();

            for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                compare_atoms(atom1, atom2);
            }
        }
    }

    #[test]
    fn read_trr_step_3() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        for (raw1, raw2) in system1
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .with_step(3)
            .unwrap()
            .zip(
                system2
                    .trr_iter("test_files/short_trajectory.trr")
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
    fn read_trr_step_3_double_precision() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        for (raw1, raw2) in system1
            .trr_iter("test_files/short_trajectory_double.trr")
            .unwrap()
            .with_step(3)
            .unwrap()
            .zip(
                system2
                    .trr_iter("test_files/short_trajectory_double.trr")
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
    fn read_trr_step_23() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        for (raw1, raw2) in system1
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .with_step(23)
            .unwrap()
            .zip(
                system2
                    .trr_iter("test_files/short_trajectory.trr")
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
    fn read_trr_step_0() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .with_step(0)
        {
            Ok(_) => panic!("Should have failed."),
            Err(ReadTrajError::InvalidStep(s)) => assert_eq!(s, 0),
            Err(e) => panic!("Incorrect error type {} returned.", e),
        }
    }

    #[test]
    fn read_trr_range_step() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        let mut i = 0;

        for (raw_frame1, raw_frame2) in system1
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .skip(3)
            .step_by(2)
            .zip(
                system2
                    .trr_iter("test_files/short_trajectory.trr")
                    .unwrap()
                    .with_range(200.0, 600.0)
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
    fn read_trr_range_step_step_range() {
        let mut system1 = System::from_file("test_files/example.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        let mut i = 0;

        for (raw_frame1, raw_frame2) in system1
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .with_step(2)
            .unwrap()
            .with_range(200.0, 600.0)
            .unwrap()
            .zip(
                system2
                    .trr_iter("test_files/short_trajectory.trr")
                    .unwrap()
                    .with_range(200.0, 600.0)
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
    fn read_trr_unmatching() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();

        match system.trr_iter("test_files/short_trajectory.trr") {
            Err(ReadTrajError::AtomsNumberMismatch(_)) => (),
            _ => panic!("TRR file should not be valid."),
        }
    }

    #[test]
    fn read_trr_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.trr_iter("test_files/nonexistent.trr") {
            Err(ReadTrajError::FileNotFound(_)) => (),
            _ => panic!("TRR file should not exist."),
        }
    }

    fn check_trr(
        system: &mut System,
        path_to_output: impl AsRef<Path>,
        check_filename: impl AsRef<Path>,
        check_trr_iter: impl AsRef<Path>,
    ) {
        let mut system2 = System::from_file(check_filename).unwrap();

        for (raw_frame1, raw_frame2) in system
            .trr_iter(check_trr_iter)
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

            let box1 = frame1.get_box_as_ref().unwrap();
            let box2 = frame2.get_box_as_ref().unwrap();
            assert_approx_eq!(f32, box1.x, box2.x);
            assert_approx_eq!(f32, box1.y, box2.y);
            assert_approx_eq!(f32, box1.z, box2.z);

            for (a1, a2) in frame1
                .get_atoms_as_ref()
                .iter()
                .zip(frame2.get_atoms_as_ref().iter())
            {
                compare_atoms(a1, a2);
            }
        }
    }

    #[test]
    fn read_trr_triclinic() {
        let mut system = System::from_file("test_files/triclinic.gro").unwrap();

        // second frame
        let frame = system
            .trr_iter("test_files/triclinic_trajectory.trr")
            .unwrap()
            .nth(1)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 6000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 120.0);

        let simbox = frame.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 5.2841544);
        assert_approx_eq!(f32, simbox.v2y, 4.7775064);
        assert_approx_eq!(f32, simbox.v3z, 2.2274022);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 0.8424022);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, 1.0153817);
        assert_approx_eq!(f32, simbox.v3y, -1.6863307);

        // last frame
        let frame = system
            .trr_iter("test_files/triclinic_trajectory.trr")
            .unwrap()
            .nth(12)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 48000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 960.0);

        let simbox = frame.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 5.2607775);
        assert_approx_eq!(f32, simbox.v2y, 4.7563710);
        assert_approx_eq!(f32, simbox.v3z, 2.1813555);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 0.8386754);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, 1.0108896);
        assert_approx_eq!(f32, simbox.v3y, -1.6788703);
    }

    #[test]
    fn read_trr_octahedron() {
        let mut system = System::from_file("test_files/octahedron.gro").unwrap();

        // second frame
        let frame = system
            .trr_iter("test_files/octahedron_trajectory.trr")
            .unwrap()
            .nth(1)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 6000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 120.0);

        let simbox = frame.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 6.2648335);
        assert_approx_eq!(f32, simbox.v2y, 5.9065430);
        assert_approx_eq!(f32, simbox.v3z, 5.1107280);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 2.0882778);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, -2.0882778);
        assert_approx_eq!(f32, simbox.v3y, 2.9532664);

        // last frame
        let frame = system
            .trr_iter("test_files/octahedron_trajectory.trr")
            .unwrap()
            .nth(12)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 48000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 960.0);

        let simbox = frame.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 6.2079663);
        assert_approx_eq!(f32, simbox.v2y, 5.8529277);
        assert_approx_eq!(f32, simbox.v3z, 5.0791510);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 2.0693220);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, -2.0693220);
        assert_approx_eq!(f32, simbox.v3y, 2.9264590);
    }

    #[test]
    fn read_trr_dodecahedron() {
        let mut system = System::from_file("test_files/dodecahedron.gro").unwrap();

        // second frame
        let frame = system
            .trr_iter("test_files/dodecahedron_trajectory.trr")
            .unwrap()
            .nth(1)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 6000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 120.0);

        let simbox = frame.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 6.2602970);
        assert_approx_eq!(f32, simbox.v2y, 6.2602970);
        assert_approx_eq!(f32, simbox.v3z, 4.4314090);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 0.0000000);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, 3.1301484);
        assert_approx_eq!(f32, simbox.v3y, 3.1301484);

        // last frame
        let frame = system
            .trr_iter("test_files/dodecahedron_trajectory.trr")
            .unwrap()
            .nth(12)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 48000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 960.0);

        let simbox = frame.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 6.2228966);
        assert_approx_eq!(f32, simbox.v2y, 6.2228966);
        assert_approx_eq!(f32, simbox.v3z, 4.4067430);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 0.0000000);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, 3.1114483);
        assert_approx_eq!(f32, simbox.v3y, 3.1114483);
    }

    #[test]
    fn read_trr_triclinic_double_precision() {
        let mut system = System::from_file("test_files/triclinic.gro").unwrap();

        // second frame
        let frame = system
            .trr_iter("test_files/triclinic_trajectory_double_precision.trr")
            .unwrap()
            .nth(1)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 6000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 120.0);

        let simbox = frame.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 5.2975800);
        assert_approx_eq!(f32, simbox.v2y, 4.7896442);
        assert_approx_eq!(f32, simbox.v3z, 2.1716056);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 0.8445424);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, 1.0179615);
        assert_approx_eq!(f32, simbox.v3y, -1.6906150);

        // last frame
        let frame = system
            .trr_iter("test_files/triclinic_trajectory_double_precision.trr")
            .unwrap()
            .nth(12)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 48000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 960.0);

        let simbox = frame.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 5.1767554);
        assert_approx_eq!(f32, simbox.v2y, 4.6804047);
        assert_approx_eq!(f32, simbox.v3z, 2.0012338);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 0.8252806);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, 0.9947443);
        assert_approx_eq!(f32, simbox.v3y, -1.6520565);
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

        check_trr(
            &mut system,
            path_to_output,
            "test_files/example.gro",
            "test_files/short_trajectory.trr",
        );
    }

    #[test]
    fn write_invalid_path() {
        let system = System::from_file("test_files/example.gro").unwrap();

        match TrrWriter::new(&system, "test_files/nonexistent/output.trr") {
            Err(WriteTrajError::CouldNotCreate(_)) => (),
            _ => panic!("Output TRR file should not have been created."),
        }
    }

    #[test]
    fn write_group_trr() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        let mut writer = TrrGroupWriter::new(&system, "Protein", path_to_output).unwrap();

        for _ in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            writer.write_frame().unwrap();
        }

        // we must close the file, otherwise metadata do not get updated
        drop(writer);

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory_protein.trr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group_trr_all() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        let mut writer = TrrGroupWriter::new(&system, "all", path_to_output).unwrap();

        for _ in system.trr_iter("test_files/short_trajectory.trr").unwrap() {
            writer.write_frame().unwrap();
        }

        // we must close the file, otherwise metadata do not get updated
        drop(writer);

        check_trr(
            &mut system,
            path_to_output,
            "test_files/example.gro",
            "test_files/short_trajectory.trr",
        );
    }

    #[test]
    fn write_group_trr_replace() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        let mut writer = TrrGroupWriter::new(&system, "Protein", path_to_output).unwrap();

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
        let mut expected = File::open("test_files/short_trajectory_protein.trr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group_trr_remove() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        let mut writer = TrrGroupWriter::new(&system, "Protein", path_to_output).unwrap();

        // remove the protein group from the system; this should not change the output of the XtcGroupWriter
        unsafe {
            let val = system.get_groups_as_ref_mut().remove("Protein").unwrap();
            assert_eq!(val.get_atoms(), writer.group.get_atoms());
            assert!(!system.group_exists("Protein"));
        }

        for _ in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            writer.write_frame().unwrap();
        }

        // we must close the file, otherwise metadata do not get updated
        drop(writer);

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory_protein.trr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group_trr_nonexistent() {
        let system = System::from_file("test_files/example.gro").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        match TrrGroupWriter::new(&system, "Protein", path_to_output) {
            Err(WriteTrajError::GroupNotFound(g)) => assert_eq!(g, "Protein"),
            _ => panic!("Output XTC file should not have been created."),
        }
    }

    #[test]
    fn write_trr_triclinic() {
        let mut system = System::from_file("test_files/triclinic.gro").unwrap();

        let trr_output = NamedTempFile::new().unwrap();
        let path_to_output = trr_output.path();

        let mut writer = TrrWriter::new(&system, path_to_output).unwrap();

        for frame in system
            .trr_iter("test_files/triclinic_trajectory.trr")
            .unwrap()
        {
            let _ = frame.unwrap();

            writer.write_frame().unwrap();
        }

        drop(writer);

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/triclinic_trajectory_full.trr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_trr_octahedron() {
        let mut system = System::from_file("test_files/octahedron.gro").unwrap();

        let trr_output = NamedTempFile::new().unwrap();
        let path_to_output = trr_output.path();

        let mut writer = TrrWriter::new(&system, path_to_output).unwrap();

        for frame in system
            .trr_iter("test_files/octahedron_trajectory.trr")
            .unwrap()
        {
            let _ = frame.unwrap();

            writer.write_frame().unwrap();
        }

        drop(writer);

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/octahedron_trajectory_full.trr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_trr_dodecahedron() {
        let mut system = System::from_file("test_files/dodecahedron.gro").unwrap();

        let trr_output = NamedTempFile::new().unwrap();
        let path_to_output = trr_output.path();

        let mut writer = TrrWriter::new(&system, path_to_output).unwrap();

        for frame in system
            .trr_iter("test_files/dodecahedron_trajectory.trr")
            .unwrap()
        {
            let _ = frame.unwrap();

            writer.write_frame().unwrap();
        }

        drop(writer);

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/dodecahedron_trajectory_full.trr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
}
