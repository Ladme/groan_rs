// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of functions for reading and writing trr files.

use std::marker::PhantomData;
use std::os::raw::{c_float, c_int};
use std::path::Path;

use crate::errors::{ReadTrajError, TrajError, WriteTrajError};
use crate::io::traj_cat::TrajConcatenator;
use crate::io::traj_read::{
    FrameData, FrameDataTime, TrajRangeRead, TrajRead, TrajReadOpen, TrajReader, TrajStepRead,
    TrajStepTimeRead,
};
use crate::io::xdrfile::{self, OpenMode, XdrFile};
use crate::prelude::AtomIterator;
use crate::structures::group::Group;
use crate::structures::{simbox::SimBox, vector3d::Vector3D};
use crate::system::System;

use super::traj_write::{PrivateTrajWrite, TrajWrite};

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
    type TrajFile = XdrFile;

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
                xdrfile.handle,
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
            for (i, atom) in system.get_atoms_mut().iter_mut().enumerate() {
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

    fn get_system(&mut self) -> *mut System {
        self.system
    }

    fn get_file_handle(&mut self) -> &mut XdrFile {
        &mut self.trr
    }
}

impl<'a> TrajReadOpen<'a> for TrrReader<'a> {
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
}

impl<'a> TrajRangeRead<'a> for TrrReader<'a> {
    fn jump_to_start(&mut self, start_time: f32) -> Result<(), ReadTrajError> {
        unsafe {
            if xdrfile::trr_jump_to_start(self.get_file_handle().handle, start_time as c_float) != 0
            {
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
            match xdrfile::trr_skip_frame(self.get_file_handle().handle) {
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

impl<'a> TrajStepTimeRead<'a> for TrrReader<'a> {
    fn skip_frame_time(&mut self) -> Result<Option<f32>, ReadTrajError> {
        unsafe {
            let mut time: c_float = 0.0;

            match xdrfile::trr_skip_frame_with_time(self.get_file_handle().handle, &mut time as *mut c_float) {
                0 => Ok(Some(time)),
                1 => Err(ReadTrajError::SkipFailed),
                2 => Ok(None),
                number => panic!(
                    "FATAL GROAN ERROR | TrrReader::skip_frame_time | `xdrfile::trr_skip_frame_with_time` returned an unsupported number '{}'",
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
    /// ## Notes
    /// - The function checks whether the number of atoms in the system corresponds to the number of atoms in the trr file.
    /// - The `System` structure is modified while iterating through the trr file.
    /// - The trr file does not need to have positions, velocities, and forces provided for each frame.
    ///   In case any of these properties is missing, it is set to `None` for all atoms.
    pub fn trr_iter(
        &mut self,
        filename: impl AsRef<Path>,
    ) -> Result<TrajReader<TrrReader>, ReadTrajError> {
        Ok(TrajReader::wrap_traj(TrrReader::new(self, filename)?))
    }

    /// Iterate through multiple trr files.
    /// Any duplicate frames at the boundaries of the trajectories are skipped.
    ///
    /// ## Example
    /// Iterate through multiple trr files and calculate and print the current center of geometry of the system.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use groan_rs::errors::ReadTrajError;
    /// #
    /// fn example_fn() -> Result<(), ReadTrajError> {
    ///     let mut system = System::from_file("system.gro").unwrap();
    ///     let trajectories = vec!["md0001.trr", "md0002.trr", "md0003.trr"];     
    ///
    ///     for raw_frame in system.trr_cat_iter(&trajectories).unwrap() {
    ///         let frame = raw_frame?;
    ///         println!("{:?}", frame.group_get_center("all"));
    ///     }
    ///
    ///     Ok(())
    /// }
    /// ```
    ///
    /// You can also iterate through only part of the trajectory, only read every Nth trajectory frame,
    /// and print progress of the trajectory reading. For more information, see [`System::trr_iter`](`System::trr_iter`).
    ///
    /// For more information about trajectory concatenation, see [`System::traj_cat_iter`](`System::traj_cat_iter`).
    pub fn trr_cat_iter<'a>(
        &'a mut self,
        filenames: &[impl AsRef<Path>],
    ) -> Result<TrajReader<'a, TrajConcatenator<'a, TrrReader>>, ReadTrajError> {
        self.traj_cat_iter::<TrrReader>(filenames)
    }
}

/**************************/
/*       WRITING TRR      */
/**************************/

impl System {
    /// Initializes a TRR trajectory writer and associates it with `System`.
    ///
    /// This is a convenience method for [`System::traj_writer_init`] with `TrrWriter`, writing in TRR format.
    #[inline(always)]
    pub fn trr_writer_init(&mut self, filename: impl AsRef<Path>) -> Result<(), WriteTrajError> {
        self.traj_writer_init::<TrrWriter>(filename)
    }

    /// Initializes a TRR trajectory writer for a specific group of atoms within `System`.
    ///
    /// This is a convenience method for [`System::traj_group_writer_init`] with `TrrWriter`, writing in TRR format.
    #[inline(always)]
    pub fn trr_group_writer_init(
        &mut self,
        filename: impl AsRef<Path>,
        group: &str,
    ) -> Result<(), WriteTrajError> {
        self.traj_group_writer_init::<TrrWriter>(filename, group)
    }
}

pub struct TrrWriter {
    trr: XdrFile,
    // deep copy of the group from `System`
    group: Group,
}

impl TrajWrite for TrrWriter {}

impl PrivateTrajWrite for TrrWriter {
    fn new(
        system: &System,
        filename: impl AsRef<Path>,
        group: Option<&str>,
    ) -> Result<Self, WriteTrajError>
    where
        Self: Sized,
    {
        // get the requested group from the system or use `all`
        // this has to be done before opening the trr file
        let group = match group {
            Some(x) => system
                .get_groups()
                .get(x)
                .ok_or_else(|| WriteTrajError::GroupNotFound(x.to_owned()))?
                .clone(),
            None => system
                .get_groups()
                .get("all")
                .expect("FATAL GROAN ERROR | TrrWriter::new | Group `all` should exist.")
                .clone(),
        };

        // create the trr file and save a handle to it
        let trr = match XdrFile::open_xdr(filename.as_ref(), OpenMode::Write) {
            Ok(x) => x,
            Err(TrajError::FileNotFound(x)) => return Err(WriteTrajError::CouldNotCreate(x)),
            Err(TrajError::InvalidPath(x)) => return Err(WriteTrajError::InvalidPath(x)),
        };

        Ok(Self { trr, group })
    }

    fn write_frame(&mut self, system: &System) -> Result<(), WriteTrajError> {
        let n_atoms = self.group.get_n_atoms();

        // prepare coordinate, velocity and forces matrix
        let mut coordinates = vec![[0.0, 0.0, 0.0]; n_atoms];
        let mut velocities = vec![[0.0, 0.0, 0.0]; n_atoms];
        let mut forces = vec![[0.0, 0.0, 0.0]; n_atoms];

        let iterator =
            AtomIterator::new(system.get_atoms(), self.group.get_atoms(), system.get_box());

        for atom in iterator {
            if let Some(pos) = atom.get_position() {
                coordinates[atom.get_index()] = [pos.x, pos.y, pos.z];
            }

            if let Some(vel) = atom.get_velocity() {
                velocities[atom.get_index()] = [vel.x, vel.y, vel.z];
            }

            if let Some(force) = atom.get_force() {
                forces[atom.get_index()] = [force.x, force.y, force.z]
            }
        }

        // write the trr frame
        let return_code = unsafe {
            xdrfile::write_trr(
                self.trr.handle,
                n_atoms as c_int,
                system.get_simulation_step() as i32,
                system.get_simulation_time(),
                system.get_lambda(),
                &mut xdrfile::simbox2matrix(system.get_box()),
                coordinates.as_mut_ptr(),
                velocities.as_mut_ptr(),
                forces.as_mut_ptr(),
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
    use crate::test_utilities::utilities::{compare_atoms, compare_box};
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
        let simbox = system.get_box().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.01331);
        assert_approx_eq!(f32, simbox.y, 13.01331);
        assert_approx_eq!(f32, simbox.z, 11.25347);

        let atom1 = &system.get_atoms()[0];
        let atom2 = &system.get_atoms()[16843];

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
        let simbox = system.get_box().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.024242);
        assert_approx_eq!(f32, simbox.y, 13.024242);
        assert_approx_eq!(f32, simbox.z, 11.242146);

        let atom1 = &system.get_atoms()[0];
        let atom2 = &system.get_atoms()[16843];

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
        let simbox = system.get_box().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.076236);
        assert_approx_eq!(f32, simbox.y, 13.076236);
        assert_approx_eq!(f32, simbox.z, 11.13604);

        let atom1 = &system.get_atoms()[0];
        let atom2 = &system.get_atoms()[16843];

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
        let simbox = system.get_box().unwrap();
        assert_approx_eq!(f32, simbox.x, 13.083817);
        assert_approx_eq!(f32, simbox.y, 13.083817);
        assert_approx_eq!(f32, simbox.z, 11.159238);

        let atom1 = &system.get_atoms()[0];
        let atom2 = &system.get_atoms()[16843];

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
        let simbox = system.get_box().unwrap();
        assert_approx_eq!(f32, simbox.x, 12.965868);
        assert_approx_eq!(f32, simbox.y, 12.965868);
        assert_approx_eq!(f32, simbox.z, 11.348931);

        let atom1 = &system.get_atoms()[0];
        let atom2 = &system.get_atoms()[16843];

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

        let atom1 = &system.get_atoms()[0];
        let atom2 = &system.get_atoms()[16843];

        assert_eq!(system.get_simulation_step(), 0);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 0.0);
        let simbox = system.get_box().unwrap();
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

        let atom1 = &system.get_atoms()[0];
        let atom2 = &system.get_atoms()[16843];

        assert_eq!(system.get_simulation_step(), 6000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 120.0);
        let simbox = system.get_box().unwrap();
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

        let atom1 = &system.get_atoms()[0];
        let atom2 = &system.get_atoms()[16843];

        assert_eq!(system.get_simulation_step(), 8000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 160.0);
        let simbox = system.get_box().unwrap();
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

        let atom1 = &system.get_atoms()[0];
        let atom2 = &system.get_atoms()[16843];

        assert_eq!(system.get_simulation_step(), 12000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 240.0);
        let simbox = system.get_box().unwrap();
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

        let atom1 = &system.get_atoms()[0];
        let atom2 = &system.get_atoms()[16843];

        assert_eq!(system.get_simulation_step(), 32000);
        assert_eq!(system.get_lambda(), 0.0);
        assert_approx_eq!(f32, system.get_simulation_time(), 640.0);
        let simbox = system.get_box().unwrap();
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

            let box1 = frame1.get_box().unwrap();
            let box2 = frame2.get_box().unwrap();
            assert_approx_eq!(f32, box1.x, box2.x);
            assert_approx_eq!(f32, box1.y, box2.y);
            assert_approx_eq!(f32, box1.z, box2.z);

            for (a1, a2) in frame1.get_atoms().iter().zip(frame2.get_atoms().iter()) {
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

        let simbox = frame.get_box().unwrap();
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

        let simbox = frame.get_box().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 5.2607775);
        assert_approx_eq!(f32, simbox.v2y, 4.756371);
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

        let simbox = frame.get_box().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 6.2648335);
        assert_approx_eq!(f32, simbox.v2y, 5.906543);
        assert_approx_eq!(f32, simbox.v3z, 5.110728);
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

        let simbox = frame.get_box().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 6.2079663);
        assert_approx_eq!(f32, simbox.v2y, 5.8529277);
        assert_approx_eq!(f32, simbox.v3z, 5.079151);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 2.069322);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, -2.069322);
        assert_approx_eq!(f32, simbox.v3y, 2.926459);
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

        let simbox = frame.get_box().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 6.260297);
        assert_approx_eq!(f32, simbox.v2y, 6.260297);
        assert_approx_eq!(f32, simbox.v3z, 4.431409);
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

        let simbox = frame.get_box().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 6.2228966);
        assert_approx_eq!(f32, simbox.v2y, 6.2228966);
        assert_approx_eq!(f32, simbox.v3z, 4.406743);
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

        let simbox = frame.get_box().unwrap();
        assert_approx_eq!(f32, simbox.v1x, 5.29758);
        assert_approx_eq!(f32, simbox.v2y, 4.7896442);
        assert_approx_eq!(f32, simbox.v3z, 2.1716056);
        assert_approx_eq!(f32, simbox.v1y, 0.0000000);
        assert_approx_eq!(f32, simbox.v1z, 0.0000000);
        assert_approx_eq!(f32, simbox.v2x, 0.8445424);
        assert_approx_eq!(f32, simbox.v2z, 0.0000000);
        assert_approx_eq!(f32, simbox.v3x, 1.0179615);
        assert_approx_eq!(f32, simbox.v3y, -1.690615);

        // last frame
        let frame = system
            .trr_iter("test_files/triclinic_trajectory_double_precision.trr")
            .unwrap()
            .nth(12)
            .unwrap()
            .unwrap();
        assert_eq!(frame.get_simulation_step(), 48000);
        assert_approx_eq!(f32, frame.get_simulation_time(), 960.0);

        let simbox = frame.get_box().unwrap();
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
    fn cat_trr() {
        let mut system_single = System::from_file("test_files/example.gro").unwrap();
        let mut system_cat = System::from_file("test_files/example.gro").unwrap();

        let traj_single = system_single
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap();
        let traj_cat = system_cat
            .trr_cat_iter(&[
                "test_files/split/traj1.trr",
                "test_files/split/traj2.trr",
                "test_files/split/traj3.trr",
                "test_files/split/traj4.trr",
                "test_files/split/traj5.trr",
                "test_files/split/traj6.trr",
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
    fn write_trr() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        system.trr_writer_init(path_to_output).unwrap();

        for frame in system.trr_iter("test_files/short_trajectory.trr").unwrap() {
            let frame = frame.unwrap();

            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

        check_trr(
            &mut system,
            path_to_output,
            "test_files/example.gro",
            "test_files/short_trajectory.trr",
        );
    }

    #[test]
    fn write_invalid_path() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.trr_writer_init("test_files/nonexistent/output.trr") {
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

        system
            .trr_group_writer_init(path_to_output, "Protein")
            .unwrap();

        for frame in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            let frame = frame.unwrap();

            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory_protein.trr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group_trr_all() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        system.trr_group_writer_init(path_to_output, "all").unwrap();

        for frame in system.trr_iter("test_files/short_trajectory.trr").unwrap() {
            let frame = frame.unwrap();

            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

        check_trr(
            &mut system,
            path_to_output,
            "test_files/example.gro",
            "test_files/short_trajectory.trr",
        );
    }

    #[test]
    fn write_group_trr_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.trr_group_writer_init("will_not_be_created.trr", "Protein") {
            Err(WriteTrajError::GroupNotFound(g)) => {
                assert_eq!(g, "Protein");
                assert!(!Path::new("will_not_be_created.trr").exists());
            }
            _ => panic!("Output XTC file should not have been created."),
        }
    }

    #[test]
    fn write_group_trr_replace() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        system
            .trr_group_writer_init(path_to_output, "Protein")
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
        let mut expected = File::open("test_files/short_trajectory_protein.trr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group_trr_remove() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();

        system
            .trr_group_writer_init(path_to_output, "Protein")
            .unwrap();

        // remove the `Protein` group; this should not change the output of the trajectory writing
        system.group_remove("Protein").unwrap();

        for frame in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            let frame = frame.unwrap();
            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory_protein.trr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_trr_triclinic() {
        let mut system = System::from_file("test_files/triclinic.gro").unwrap();

        let trr_output = NamedTempFile::new().unwrap();
        let path_to_output = trr_output.path();

        system.trr_writer_init(path_to_output).unwrap();

        for frame in system
            .trr_iter("test_files/triclinic_trajectory.trr")
            .unwrap()
        {
            let frame = frame.unwrap();

            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/triclinic_trajectory_full.trr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_trr_octahedron() {
        let mut system = System::from_file("test_files/octahedron.gro").unwrap();

        let trr_output = NamedTempFile::new().unwrap();
        let path_to_output = trr_output.path();

        system.trr_writer_init(path_to_output).unwrap();

        for frame in system
            .trr_iter("test_files/octahedron_trajectory.trr")
            .unwrap()
        {
            let frame = frame.unwrap();

            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/octahedron_trajectory_full.trr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_trr_dodecahedron() {
        let mut system = System::from_file("test_files/dodecahedron.gro").unwrap();

        let trr_output = NamedTempFile::new().unwrap();
        let path_to_output = trr_output.path();

        system.trr_writer_init(path_to_output).unwrap();

        for frame in system
            .trr_iter("test_files/dodecahedron_trajectory.trr")
            .unwrap()
        {
            let frame = frame.unwrap();

            frame.traj_write_frame().unwrap();
        }

        system.traj_close();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/dodecahedron_trajectory_full.trr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
}
