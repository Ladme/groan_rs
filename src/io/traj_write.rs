// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Traits and structures for writing trajectory files.

/**************************/
/*  WRITING TRAJECTORIES  */
/**************************/

use crate::files::FileType;
use crate::{errors::WriteTrajError, system::System};
use hashbrown::HashMap;
use std::fmt::Debug;
use std::path::Path;

#[cfg(feature = "parallel")]
pub(crate) use multi_threaded::SystemWriters;
#[cfg(not(feature = "parallel"))]
pub(crate) use single_threaded::SystemWriters;

use super::gro_io::GroWriter;
use super::trr_io::TrrWriter;
use super::xtc_io::XtcWriter;

/// ## Associating trajectory writers with System.
impl System {
    /// Initializes a file for trajectory writing and associates the resulting writer with `System`.
    ///
    /// ## Parameters
    /// - `filename`: The path to the output file that should be opened.
    /// - `Writer`: A type implementing `TrajWrite`, determining the output format and behavior.
    ///
    /// ## Returns
    /// - `Ok` if the trajectory writer is successfully created and associated with the `System`.
    /// - `WriteTrajError` if an error occurs, such as if a writer for the same file already exists.
    ///
    /// ## Examples
    /// To start writing trajectories, initialize a trajectory writer:
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.traj_writer_init::<XtcWriter>("output.xtc").unwrap();
    /// ```
    ///
    /// After initialization, an XTC trajectory writer is now associated with the system, targeting `output.xtc`.
    /// Additional trajectory writers can be attached to the same system:
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # let mut system = System::from_file("system.gro").unwrap();
    /// #
    /// // attach a TRR trajectory writer
    /// system.traj_writer_init::<TrrWriter>("output.trr").unwrap();
    /// // attach an XTC trajectory writer for a specific group
    /// system.group_create("Protein", "@protein").unwrap();
    /// system.traj_group_writer_init::<XtcWriter>("output_protein.xtc", "Protein").unwrap();
    /// ```
    ///
    /// Now three writers are associated with the system. To write the current system state into each writer:
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # let mut system = System::from_file("system.gro").unwrap();
    /// # system.traj_writer_init::<XtcWriter>("output.xtc").unwrap();
    /// #
    /// system.traj_write_frame().unwrap();
    /// ```
    ///
    /// This writes the current system state to `output.xtc` (XTC format) and `output.trr` (TRR format), and
    /// writes only the `Protein` group state to `output_protein.xtc`.
    ///
    /// To write the current system state to a specific writer, specify the writer’s filename:
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # let mut system = System::from_file("system.gro").unwrap();
    /// # system.traj_writer_init::<XtcWriter>("output.xtc").unwrap();
    /// #
    /// system.traj_write_frame_to_file("output.xtc").unwrap();
    /// ```
    /// This writes a frame only to `output.xtc` in XTC format.
    ///
    /// Trajectory files close automatically when the `System` goes out of scope. To manually close all writers:
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # let mut system = System::from_file("system.gro").unwrap();
    /// # system.traj_writer_init::<XtcWriter>("output.xtc").unwrap();
    /// #
    /// system.traj_close();
    /// ```
    ///
    /// To close a specific writer by its filename:
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # let mut system = System::from_file("system.gro").unwrap();
    /// # system.traj_writer_init::<XtcWriter>("output.xtc").unwrap();
    /// #
    /// system.traj_close_file("output.xtc");
    /// ```
    ///
    /// ## Notes
    /// - Multiple writers can be associated with the same system.
    /// - Each writer has a unique name based on the output file path. Attempting to add a writer for the same
    ///   file will result in an error.
    /// - See also [System::traj_writer_auto_init], [System::xtc_writer_init], [System::trr_writer_init], and [System::gro_writer_init].
    #[inline]
    pub fn traj_writer_init<Writer>(
        &mut self,
        filename: impl AsRef<Path>,
    ) -> Result<(), WriteTrajError>
    where
        Writer: TrajWrite + 'static,
    {
        let writers = self.get_writers_mut() as *mut SystemWriters;
        unsafe { &mut (*writers) }.try_insert::<Writer>(self, filename, None)
    }

    /// Initializes a trajectory writer for a specific group of atoms and associates it with `System`.
    ///
    /// This method functions similarly to [`System::traj_writer_init`], but restricts the writer to a specific
    /// group of atoms within the system.
    ///
    /// ## Parameters
    /// - `filename`: The path to the output file that should be opened.
    /// - `group`: The name of the atom group for which the trajectory writer will be initialized. This group
    ///   must already exist within the system.
    /// - `Writer`: A type implementing `TrajWrite`, determining the output format and behavior.
    ///
    /// ## Returns
    /// - `Ok` if the trajectory writer is successfully created for the specified group and associated
    ///   with the `System`.
    /// - `WriteTrajError` if an error occurs, such as if a writer for the same file already exists or
    ///    if the group does not exist.
    ///
    /// ## Notes
    /// - Once the writer is initialized, the group can be removed from the system. The trajectory writer retains
    ///   an independent copy of the group’s state, allowing it to function even if the original group is deleted.
    /// - For more details on general trajectory writer initialization and usage, see [`System::traj_writer_init`].
    #[inline]
    pub fn traj_group_writer_init<Writer>(
        &mut self,
        filename: impl AsRef<Path>,
        group: &str,
    ) -> Result<(), WriteTrajError>
    where
        Writer: TrajWrite + 'static,
    {
        let writers = self.get_writers_mut() as *mut SystemWriters;
        unsafe { &mut (*writers) }.try_insert::<Writer>(self, filename, Some(group))
    }

    /// Automatically initializes a trajectory writer based on the output file
    /// extension and associates it with `System`.
    ///
    /// This method determines the appropriate writer type (`XtcWriter`, `TrrWriter`, or `GroWriter`)
    /// based on the file extension and calls [`System::traj_writer_init`].
    /// Returns an error if the extension is unknown or unsupported.
    #[inline]
    pub fn traj_writer_auto_init(
        &mut self,
        filename: impl AsRef<Path>,
    ) -> Result<(), WriteTrajError> {
        match FileType::from_name(&filename) {
            FileType::XTC => self.traj_writer_init::<XtcWriter>(filename),
            FileType::TRR => self.traj_writer_init::<TrrWriter>(filename),
            FileType::GRO => self.traj_writer_init::<GroWriter>(filename),
            _ => Err(WriteTrajError::UnknownExtension(Box::from(
                filename.as_ref(),
            ))),
        }
    }

    /// Automatically initializes a trajectory writer for a specific group
    /// based on the output file extension and associates it with `System`.
    ///
    /// This method selects the appropriate writer type (`XtcWriter`, `TrrWriter`, or `GroWriter`)
    /// based on the file extension and calls [`System::traj_group_writer_init`].
    /// Returns an error if the extension is unknown or unsupported.
    #[inline]
    pub fn traj_group_writer_auto_init(
        &mut self,
        filename: impl AsRef<Path>,
        group: &str,
    ) -> Result<(), WriteTrajError> {
        match FileType::from_name(&filename) {
            FileType::XTC => self.traj_group_writer_init::<XtcWriter>(filename, group),
            FileType::TRR => self.traj_group_writer_init::<TrrWriter>(filename, group),
            FileType::GRO => self.traj_group_writer_init::<GroWriter>(filename, group),
            _ => Err(WriteTrajError::UnknownExtension(Box::from(
                filename.as_ref(),
            ))),
        }
    }

    /// Writes a single frame to a specified trajectory writer identified by its filename.
    ///
    /// This function writes the current state of the `System`
    /// to the writer associated with `name`.
    /// If no writer matches `name`, an error is returned.
    /// For examples, see [`System::traj_writer_init`].
    #[inline(always)]
    pub fn traj_write_frame_to_file(
        &mut self,
        name: impl AsRef<Path>,
    ) -> Result<(), WriteTrajError> {
        let system = self as *const System;
        let writer_name = crate::auxiliary::path2string(&name);

        self.get_writers_mut()
            .select_and_write(unsafe { &*system }, &writer_name)
    }

    /// Writes a single frame to all associated trajectory writers.
    ///
    /// This function writes the current state of the `System`
    /// to each trajectory writer attached to it.
    /// For examples, see [`System::traj_writer_init`].
    #[inline(always)]
    pub fn traj_write_frame(&mut self) -> Result<(), WriteTrajError> {
        let system = self as *const System;
        self.get_writers_mut().write_all(unsafe { &*system })
    }

    /// Closes a specific trajectory writer identified by its filename.
    ///
    /// This function closes the trajectory writer associated with `name`.
    /// If no writer matches `name`, an error is returned.
    /// For examples, see [`System::traj_writer_init`].
    #[inline(always)]
    pub fn traj_close_file(&mut self, name: impl AsRef<Path>) -> Result<(), WriteTrajError> {
        let writer_name = crate::auxiliary::path2string(&name);
        self.get_writers_mut().select_and_close(&writer_name)
    }

    /// Closes all trajectory writers associated with the `System`.
    ///
    /// This function closes every trajectory writer attached to the `System`,
    /// releasing any resources they hold.
    /// For examples, see [`System::traj_writer_init`].
    #[inline(always)]
    pub fn traj_close(&mut self) {
        self.get_writers_mut().close_all()
    }
}

#[cfg(not(feature = "parallel"))]
mod single_threaded {
    use std::{cell::RefCell, rc::Rc};

    use super::*;

    #[derive(Clone, Default)]
    pub(crate) struct SystemWriters(HashMap<String, Rc<RefCell<dyn TrajWrite>>>);

    impl SystemWriters {
        /// Select a writer associated with the system and write the current frame of the system into it.
        pub(crate) fn select_and_write(
            &self,
            system: &System,
            name: &str,
        ) -> Result<(), WriteTrajError> {
            let mut writer = match self.0.get(name) {
                Some(x) => x.borrow_mut(),
                None => return Err(WriteTrajError::WriterNotFound(name.to_owned())),
            };

            writer.write_frame(system)
        }

        /// Select a writer associated with the system and close it.
        /// Returns an error if the writer does not exist.
        pub(crate) fn select_and_close(&mut self, name: &str) -> Result<(), WriteTrajError> {
            self.0
                .remove(name)
                .ok_or_else(|| WriteTrajError::WriterNotFound(name.to_owned()))
                .map(|_| ())
        }

        /// Write the current frame of the system using all the writers associated with the system.
        pub(crate) fn write_all(&self, system: &System) -> Result<(), WriteTrajError> {
            for writer in self.0.values() {
                let mut writer = writer.borrow_mut();
                writer.write_frame(system)?;
            }

            Ok(())
        }

        /// Close all writers associated with this system.
        pub(crate) fn close_all(&mut self) {
            self.0.clear()
        }

        /// Open a new file and add the corresponding trajectory writer into `SystemWriters`.
        ///
        /// If the file was already opened using another writer that is already associated with the system,
        /// it is not reopened and no writer is added into `SystemWriters`. Instead, an error is returned.
        pub(super) fn try_insert<Writer>(
            &mut self,
            system: &System,
            filename: impl AsRef<Path>,
            group: Option<&str>,
        ) -> Result<(), WriteTrajError>
        where
            Writer: TrajWrite + 'static,
        {
            let name = crate::auxiliary::path2string(&filename);

            // check whether writer to the same file already exists
            // TODO: check actual paths instead of just names
            if self.0.contains_key(&name) {
                return Err(WriteTrajError::WriterAlreadyExists(name.to_owned()));
            }

            let writer = Writer::new(system, &filename, group)?;

            self.0.insert(name, Rc::new(RefCell::new(writer)));

            Ok(())
        }

        pub(crate) fn len(&self) -> usize {
            self.0.len()
        }
    }

    impl Debug for SystemWriters {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(
                f,
                "{} associated thread-unsafe trajectory writer(s)",
                self.0.len()
            )
        }
    }
}

#[cfg(feature = "parallel")]
mod multi_threaded {
    use super::*;
    use parking_lot::Mutex;
    use std::sync::Arc;

    #[derive(Clone, Default)]
    pub(crate) struct SystemWriters(HashMap<String, Arc<Mutex<dyn TrajWrite>>>);

    unsafe impl Send for SystemWriters {}
    unsafe impl Sync for SystemWriters {}

    impl SystemWriters {
        /// Select a writer associated with the system and write the current frame of the system into it.
        pub(crate) fn select_and_write(
            &self,
            system: &System,
            name: &str,
        ) -> Result<(), WriteTrajError> {
            let mut writer = match self.0.get(name) {
                Some(x) => x.lock(),
                None => return Err(WriteTrajError::WriterNotFound(name.to_owned())),
            };

            writer.write_frame(system)
        }

        /// Select a writer associated with the system and close it.
        /// Returns an error if the writer does not exist.
        pub(crate) fn select_and_close(&mut self, name: &str) -> Result<(), WriteTrajError> {
            self.0
                .remove(name)
                .ok_or_else(|| WriteTrajError::WriterNotFound(name.to_owned()))
                .map(|_| ())
        }

        /// Write the current frame of the system using all the writers associated with the system.
        pub(crate) fn write_all(&self, system: &System) -> Result<(), WriteTrajError> {
            for writer in self.0.values() {
                let mut writer = writer.lock();
                writer.write_frame(system)?;
            }

            Ok(())
        }

        /// Close all writers associated with this system.
        pub(crate) fn close_all(&mut self) {
            self.0.clear()
        }

        /// Open a new file and add the corresponding trajectory writer into `SystemWriters`.
        ///
        /// If the file was already opened using another writer that is already associated with the system,
        /// it is not reopened and no writer is added into `SystemWriters`. Instead, an error is returned.
        pub(super) fn try_insert<Writer>(
            &mut self,
            system: &System,
            filename: impl AsRef<Path>,
            group: Option<&str>,
        ) -> Result<(), WriteTrajError>
        where
            Writer: TrajWrite + 'static,
        {
            let name = crate::auxiliary::path2string(&filename);

            // check whether writer to the same file already exists
            // TODO: check actual paths instead of just names
            if self.0.contains_key(&name) {
                return Err(WriteTrajError::WriterAlreadyExists(name.to_owned()));
            }

            let writer = Writer::new(system, &filename, group)?;

            self.0.insert(name, Arc::new(Mutex::new(writer)));

            Ok(())
        }

        pub(crate) fn len(&self) -> usize {
            self.0.len()
        }
    }

    impl Debug for SystemWriters {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(
                f,
                "{} associated thread-safe trajectory writer(s)",
                self.0.len()
            )
        }
    }
}

/// Any structure implementing the `TrajWrite` trait can be used as a trajectory writer.
#[allow(private_bounds)]
pub trait TrajWrite: PrivateTrajWrite {}

/// Trait containing private methods implemented by all trajectory writers.
pub(super) trait PrivateTrajWrite {
    /// Create a new trajectory writer.
    /// The existence of the group must be checked before the file is created.
    fn new(
        system: &System,
        filename: impl AsRef<Path>,
        group: Option<&str>,
    ) -> Result<Self, WriteTrajError>
    where
        Self: Sized;

    /// Write the current state of the system into an open trajectory file.
    fn write_frame(&mut self, system: &System) -> Result<(), WriteTrajError>;
}

#[cfg(test)]
mod tests {
    use std::{fs::File, path::PathBuf};

    use tempfile::{Builder, NamedTempFile};

    use super::*;

    fn create_named_path(extension: &str) -> (NamedTempFile, PathBuf) {
        let output = Builder::new().suffix(extension).tempfile().unwrap();
        let path = output.path().to_owned();
        (output, path)
    }

    #[test]
    fn multiple_writers() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let (_xtc, p_xtc) = create_named_path(".xtc");
        let (_xtc_group, p_xtc_group) = create_named_path(".xtc");
        let (_gro_group, p_gro_group) = create_named_path(".gro");
        let (_xtc2, p_xtc2) = create_named_path(".xtc");
        let (_trr, p_trr) = create_named_path(".trr");

        system.xtc_writer_init(&p_xtc).unwrap();
        system
            .traj_group_writer_auto_init(&p_xtc_group, "Protein")
            .unwrap();
        system
            .gro_group_writer_init(&p_gro_group, "Protein")
            .unwrap();
        system.traj_writer_auto_init(&p_xtc2).unwrap();

        system.traj_writer_auto_init(&p_trr).unwrap();

        match system.trr_writer_init(&p_xtc) {
            Ok(_) => panic!("Function should have failed."),
            Err(WriteTrajError::WriterAlreadyExists(x)) => assert_eq!(&x, p_xtc.to_str().unwrap()),
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }

        assert_eq!(system.get_n_writers(), 5);

        system.traj_close_file(&p_trr).unwrap();

        assert_eq!(system.get_n_writers(), 4);

        for (i, frame) in system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .enumerate()
        {
            let frame = frame.unwrap();

            if i < 3 {
                frame.traj_write_frame().unwrap();
            } else {
                frame.traj_write_frame_to_file(&p_xtc).unwrap();
                frame.traj_write_frame_to_file(&p_xtc_group).unwrap();
                frame.traj_write_frame_to_file(&p_xtc2).unwrap();
            }

            match frame.traj_write_frame_to_file("nonexistent.trr") {
                Ok(_) => panic!("Function should have failed."),
                Err(WriteTrajError::WriterNotFound(x)) => assert_eq!(&x, "nonexistent.trr"),
                Err(e) => panic!("Unexpected error type `{}` returned.", e),
            }

            if i == 3 {
                frame.traj_close_file(&p_gro_group).unwrap();
                assert_eq!(frame.get_n_writers(), 3);
            }
        }

        system.traj_close();
        assert_eq!(system.get_n_writers(), 0);

        let mut result = File::open(p_xtc).unwrap();
        let mut expected = File::open("test_files/short_trajectory.xtc").unwrap();
        assert!(file_diff::diff_files(&mut result, &mut expected));

        let mut result = File::open(p_xtc2).unwrap();
        let mut expected = File::open("test_files/short_trajectory.xtc").unwrap();
        assert!(file_diff::diff_files(&mut result, &mut expected));

        let mut result = File::open(p_xtc_group).unwrap();
        let mut expected = File::open("test_files/short_trajectory_protein.xtc").unwrap();
        assert!(file_diff::diff_files(&mut result, &mut expected));

        let mut result = File::open(p_gro_group).unwrap();
        let mut expected = File::open("test_files/expected_protein_trajectory.gro").unwrap();
        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
}
