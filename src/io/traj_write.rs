// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

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

/// Associating trajectory writers with the system.
impl System {
    #[inline]
    pub fn traj_writer<Writer>(&mut self, filename: impl AsRef<Path>) -> Result<(), WriteTrajError>
    where
        Writer: TrajWrite + 'static,
    {
        self.get_writers_mut().try_insert::<Writer>(filename, None)
    }

    #[inline]
    pub fn traj_group_writer<Writer>(
        &mut self,
        filename: impl AsRef<Path>,
        group: &str,
    ) -> Result<(), WriteTrajError>
    where
        Writer: TrajWrite + 'static,
    {
        // check that the group exists
        if !self.group_exists(group) {
            return Err(WriteTrajError::GroupNotFound(group.to_owned()));
        }

        self.get_writers_mut()
            .try_insert::<Writer>(filename, Some(group))
    }

    #[inline(always)]
    pub fn traj_writer_auto(&mut self, filename: impl AsRef<Path>) -> Result<(), WriteTrajError> {
        match FileType::from_name(&filename) {
            FileType::XTC => self.traj_writer::<XtcWriter>(filename),
            FileType::TRR => self.traj_writer::<TrrWriter>(filename),
            FileType::GRO => self.traj_writer::<GroWriter>(filename),
            _ => Err(WriteTrajError::UnknownExtension(Box::from(
                filename.as_ref(),
            ))),
        }
    }

    #[inline(always)]
    pub fn traj_group_writer_auto(
        &mut self,
        filename: impl AsRef<Path>,
        group: &str,
    ) -> Result<(), WriteTrajError> {
        // check that the group exists
        if !self.group_exists(group) {
            return Err(WriteTrajError::GroupNotFound(group.to_owned()));
        }

        match FileType::from_name(&filename) {
            FileType::XTC => self.traj_group_writer::<XtcWriter>(filename, group),
            FileType::TRR => self.traj_group_writer::<TrrWriter>(filename, group),
            FileType::GRO => self.traj_group_writer::<GroWriter>(filename, group),
            _ => Err(WriteTrajError::UnknownExtension(Box::from(
                filename.as_ref(),
            ))),
        }
    }

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

    #[inline(always)]
    pub fn traj_write_frame(&mut self) -> Result<(), WriteTrajError> {
        let system = self as *const System;
        self.get_writers_mut().write_all(unsafe { &*system })
    }

    #[inline(always)]
    pub fn traj_close_file(&mut self, name: impl AsRef<Path>) -> Result<(), WriteTrajError> {
        let writer_name = crate::auxiliary::path2string(&name);
        self.get_writers_mut().select_and_close(&writer_name)
    }

    #[inline(always)]
    pub fn traj_close_all(&mut self) {
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

            let writer = Writer::new(&filename, group)?;

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

            let writer = Writer::new(&filename, group)?;

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

pub trait TrajWrite {
    /// Write the current state of the system into an open trajectory file.
    fn write_frame(&mut self, system: &System) -> Result<(), WriteTrajError>;

    /// Create a new trajectory writer.
    fn new(filename: impl AsRef<Path>, group: Option<&str>) -> Result<Self, WriteTrajError>
    where
        Self: Sized;
}

#[cfg(test)]
mod tests {
    use std::{fs::File, path::PathBuf};

    use tempfile::NamedTempFile;

    use super::*;

    fn create_named_path() -> (NamedTempFile, PathBuf) {
        let output = NamedTempFile::new().unwrap();
        let path = output.path().to_owned();
        (output, path)
    }

    #[test]
    fn multiple_writers() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let (_xtc, p_xtc) = create_named_path();
        let (_xtc_group, p_xtc_group) = create_named_path();
        let (_gro_group, p_gro_group) = create_named_path();
        let (_xtc2, p_xtc2) = create_named_path();

        system.xtc_writer(&p_xtc).unwrap();
        system.xtc_group_writer(&p_xtc_group, "Protein").unwrap();
        system.gro_group_writer(&p_gro_group, "Protein").unwrap();
        system.xtc_writer(&p_xtc2).unwrap();

        match system.trr_writer(&p_xtc) {
            Ok(_) => panic!("Function should have failed."),
            Err(WriteTrajError::WriterAlreadyExists(x)) => assert_eq!(&x, p_xtc.to_str().unwrap()),
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }

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

        system.traj_close_all();
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
