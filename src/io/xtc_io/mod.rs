// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Implementation of functions for reading and writing XTC files.

#[cfg(feature = "molly")]
pub mod molly_xtc;
use std::path::Path;

#[cfg(feature = "molly")]
pub use molly_xtc::XtcReader;

#[cfg(not(feature = "molly"))]
pub mod xdrfile_xtc;
#[cfg(not(feature = "molly"))]
pub use xdrfile_xtc::XtcReader;

#[cfg(not(feature = "no-xdrfile"))]
pub use xtc_write::XtcWriter;

use crate::{
    errors::ReadTrajError,
    prelude::{TrajFullReadOpen, TrajReader},
    system::System,
};

use super::traj_cat::TrajConcatenator;

#[cfg(not(feature = "no-xdrfile"))]
use super::traj_write::PrivateTrajWrite;

/// ## Methods for reading xtc files.
impl System {
    /// Create an `XtcReader` structure which is an iterator over an xtc file.
    ///
    /// In case you are interested only in properties of a particular group of atoms,
    /// you might want to use [`System::group_xtc_iter`] instead.
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
    /// - Supports reading xtc files in the 2023 format (generally used for giant systems) **ONLY** when the `molly`
    ///   feature is enabled. **Reading of the 2023 format is not well tested. Be careful!**
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
    ) -> Result<TrajReader<'a, TrajConcatenator<'a, XtcReader<'a>>>, ReadTrajError> {
        self.traj_cat_iter::<XtcReader>(filenames)
    }
}

/**************************/
/*       WRITING XTC      */
/**************************/

#[cfg(not(feature = "no-xdrfile"))]
mod xtc_write {
    use std::{ffi::c_int, path::Path};

    use crate::{
        errors::{TrajError, WriteTrajError},
        io::xdrfile::{self, OpenMode, XdrFile},
        prelude::{AtomIterator, TrajWrite},
        structures::group::Group,
        system::System,
    };

    use super::PrivateTrajWrite;

    impl System {
        /// Initializes an XTC trajectory writer and associates it with `System`.
        ///
        /// This is a convenience method for [`System::traj_writer_init`] with `XtcWriter`, writing in XTC format.
        #[inline(always)]
        pub fn xtc_writer_init(
            &mut self,
            filename: impl AsRef<Path>,
        ) -> Result<(), WriteTrajError> {
            self.traj_writer_init::<XtcWriter>(filename)
        }

        /// Initializes an XTC trajectory writer for a specific group of atoms within `System`.
        ///
        /// This is a convenience method for [`System::traj_group_writer_init`] with `XtcWriter`, writing in XTC format.
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
            // get the requested group from the system or use `all`
            // this has to be done before opening the xtc file
            let group = match group {
                Some(x) => system
                    .get_groups()
                    .get(x)
                    .map_err(|_| WriteTrajError::GroupNotFound(x.to_owned()))?
                    .clone(),
                None => system
                    .get_groups()
                    .get("all")
                    .expect("FATAL GROAN ERROR | XtcWriter::new | Group `all` should exist.")
                    .clone(),
            };

            // create the xtc file and save a handle to it
            let xtc = match XdrFile::open_xdr(filename.as_ref(), OpenMode::Write) {
                Ok(x) => x,
                Err(TrajError::FileNotFound(x)) => return Err(WriteTrajError::CouldNotCreate(x)),
                Err(TrajError::InvalidPath(x)) => return Err(WriteTrajError::InvalidPath(x)),
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
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests_read {
    use super::*;
    use float_cmp::assert_approx_eq;

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
            println!(
                ">>>>>> {} {}",
                frame1.get_simulation_time(),
                frame2.get_simulation_time()
            );

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

    #[cfg(feature = "molly")]
    #[test]
    fn read_xtc_empty() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.xtc_iter("test_files/example_empty.gro") {
            Err(ReadTrajError::CouldNotReadMagic) => (),
            _ => panic!("Magic number reading should have failed."),
        }
    }

    #[cfg(not(feature = "molly"))]
    #[test]
    fn read_xtc_empty() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.xtc_iter("test_files/example_empty.gro") {
            Err(ReadTrajError::FileNotFound(_)) => (),
            _ => panic!("Should not be able to read the file."),
        }
    }

    #[cfg(feature = "molly")]
    #[test]
    fn read_xtc_not_xtc() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.xtc_iter("test_files/triclinic.gro") {
            Err(ReadTrajError::NotXtc(_)) => (),
            _ => panic!("This not-an-xtc-file should not have been interpreted as an xtc file."),
        }
    }

    #[cfg(not(feature = "molly"))]
    #[test]
    fn read_xtc_not_xtc() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.xtc_iter("test_files/triclinic.gro") {
            Err(ReadTrajError::FileNotFound(_)) => (),
            _ => panic!("This not-an-xtc-file should not have been interpreted as an xtc file."),
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
}

#[cfg(not(feature = "no-xdrfile"))]
#[cfg(test)]
mod tests_write {
    use std::fs::File;

    use tempfile::NamedTempFile;

    use crate::errors::WriteTrajError;

    use super::*;

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

        match system.xtc_group_writer_init("will_not_be_created.xtc", "Protein") {
            Err(WriteTrajError::GroupNotFound(g)) => {
                assert_eq!(g, "Protein");
                assert!(!Path::new("will_not_be_created.xtc").exists());
            }
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
