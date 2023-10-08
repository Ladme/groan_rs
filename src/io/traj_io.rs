// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Traits for reading generic trajectory files.

use crate::errors::{ReadTrajError, WriteTrajError};
use crate::io::xdrfile::CXdrFile;
use crate::system::general::System;
use std::marker::PhantomData;
use std::path::Path;

/**************************/
/*  READING TRAJECTORIES  */
/**************************/

/// Check that the specified times are valid. Returns Ok if valid, else ReadTrajError.
pub fn sanity_check_timerange(start_time: f32, end_time: f32) -> Result<(), ReadTrajError> {
    if start_time < 0.0 {
        return Err(ReadTrajError::TimeRangeNegative(start_time.to_string()));
    }

    if end_time < 0.0 {
        return Err(ReadTrajError::TimeRangeNegative(end_time.to_string()));
    }

    if start_time > end_time {
        return Err(ReadTrajError::InvalidTimeRange(
            start_time.to_string(),
            end_time.to_string(),
        ));
    }

    Ok(())
}

/*****************************/
/*  TrajRead and TrajReader  */
/*****************************/

/// Structure storing data from a single trajectory frame.
pub trait FrameData {
    /// Method specifying how a frame of the trajectory should be read and stored in the `FrameData` structure.
    unsafe fn from_frame(xdrfile: *mut CXdrFile, n_atoms: usize) -> Option<Result<Self, ReadTrajError>>
    where
        Self: Sized;

    /// Method specifying how the `System` structure should be updated based on the data in the `FrameData` structure.
    fn update_system(self, system: &mut System);
}

/// Any structure implementing `TrajRead` can be used to read a trajectory file.
pub trait TrajRead<'a> {
    type FrameData: FrameData;

    /// Open a trajectory file creating an iterator over it.
    ///
    /// ## Example
    /// Using `TrajRead::new` in a generic function.
    /// ```no_run
    /// use groan_rs::prelude::*;
    /// use std::path::Path;
    ///
    /// // this function can read any trajectory file implementing `TrajRead` trait
    /// fn read_trajectory<'a, Reader>(system: &'a mut System, file: impl AsRef<Path>)
    ///     where Reader: TrajRead<'a> + 'a
    /// {
    ///     // open the trajectory file for reading
    ///     let iterator = Reader::new(system, file).unwrap();
    ///
    ///     // read the trajectory file
    ///     for raw_frame in iterator {
    ///         let frame = raw_frame.unwrap();
    ///
    ///         // perform some operation with frame
    ///     }
    /// }
    ///
    /// // `read_trajectory` can be then called to read different trajectory files
    /// let mut system = System::from_file("system.gro").unwrap();
    /// read_trajectory::<XtcReader>(&mut system, "trajectory.xtc");
    /// read_trajectory::<TrrReader>(&mut system, "trajectory.trr");
    /// ```
    fn new(
        system: &'a mut System,
        filename: impl AsRef<Path>,
    ) -> Result<TrajReader<'a, Self>, ReadTrajError>
    where
        Self: Sized;

    /// Method specifying how to get a mutable pointer to `System` structure.
    fn get_system(&mut self) -> *mut System;

    /// Method specifying how to get a mutable handle to the file containing the trajectory.
    fn get_file_handle(&mut self) -> *mut CXdrFile;
}

/// Wrapper for any structure implementing `TrajRead` so the `Iterator` trait can be implemented for it.
pub struct TrajReader<'a, R: TrajRead<'a>> {
    pub traj_reader: R,
    _phantom: &'a PhantomData<R>,
}

impl<'a, R> TrajReader<'a, R>
where
    R: TrajRead<'a>,
{
    /// Pack trajectory reader implementing `TrajReader` into `TrajReader` structure.
    pub fn pack_traj(traj_reader: R) -> TrajReader<'a, R> {
        TrajReader {
            traj_reader,
            _phantom: &PhantomData,
        }
    }
}

/// Iterate the `TrajReader`.
impl<'a, R: TrajRead<'a>> Iterator for TrajReader<'a, R> {
    type Item = Result<&'a mut System, ReadTrajError>;

    fn next(&mut self) -> Option<Self::Item> {
        unsafe {
            let system = self.traj_reader.get_system();

            match R::FrameData::from_frame(
                self.traj_reader.get_file_handle(),
                (*system).get_n_atoms(),
            ) {
                None => None,
                Some(Err(e)) => Some(Err(e)),
                Some(Ok(data)) => {
                    data.update_system(&mut *system);
                    Some(Ok(&mut *system))
                }
            }
        }
    }
}

impl<'a, R> TrajReader<'a, R>
where
    R: TrajRangeRead<'a>,
{
    /// Convert `TrajReader` into `TrajRangeReader` structure.
    /// `start_time` and `end_time` must be specified in picoseconds.
    pub fn with_range(
        self,
        start_time: f32,
        end_time: f32,
    ) -> Result<TrajRangeReader<'a, R>, ReadTrajError> {
        sanity_check_timerange(start_time, end_time)?;

        let mut reader = TrajRangeReader {
            traj_reader: self.traj_reader,
            start_time,
            end_time,
            _phantom: &PhantomData,
        };

        // jump to the start of iteration
        reader.traj_reader.jump_to_start(start_time)?;

        Ok(reader)
    }
}

/***************************************/
/*  TrajRangeRead and TrajRangeReader  */
/****************************************/

/// Any structure implementing this trait can be used to construct `TrajRangeReader` structure from `TrajReader`.
pub trait TrajRangeRead<'a>: TrajRead<'a> {
    /// Method specifying how to navigate the trajectory file to reach the frame with `start_time`.
    /// In its simplest form, `skip` can be used, even though this is inefficient.
    fn jump_to_start(&mut self, start_time: f32) -> Result<(), ReadTrajError>;
}

/// Structure for partial reading of trajectory files using time ranges.
pub struct TrajRangeReader<'a, R: TrajRead<'a>> {
    pub traj_reader: R,
    pub start_time: f32,
    pub end_time: f32,
    pub _phantom: &'a PhantomData<R>,
}

/// Iterate the `TrajRangeReader`.
impl<'a, R> Iterator for TrajRangeReader<'a, R>
where
    R: TrajRead<'a>,
    R::FrameData: FrameData,
{
    type Item = Result<&'a mut System, ReadTrajError>;

    /// Read the next frame in the specified range of the trajectory.
    fn next(&mut self) -> Option<Self::Item> {
        unsafe {
            let system = self.traj_reader.get_system();

            match R::FrameData::from_frame(
                self.traj_reader.get_file_handle(),
                (*system).get_n_atoms(),
            ) {
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

/**************************/
/*  WRITING TRAJECTORIES  */
/**************************/

/// Any structure implementing `TrajWrite` can be used to write a trajectory file.
pub trait TrajWrite {
    /// Open a new trajectory file for writing.
    ///
    /// ## Example
    /// Using `TrajWrite::new` in a generic function.
    /// ```no_run
    /// use groan_rs::prelude::*;
    /// use std::path::Path;
    ///
    /// // this function can write any trajectory file implementing `TrajWrite` trait
    /// fn write_trajectory<Writer>(system: &System, file: impl AsRef<Path>)
    ///     where Writer: TrajWrite
    /// {
    ///     // open the trajectory file for writing
    ///     let mut writer = Writer::new(system, file).unwrap();
    ///
    ///     // write frame into the trajectory file
    ///     writer.write_frame().unwrap();
    /// }
    ///
    /// // `write_trajectory` can be then called to write different trajectory files
    /// let system = System::from_file("system.gro").unwrap();
    /// write_trajectory::<XtcWriter>(&system, "trajectory.xtc");
    /// write_trajectory::<TrrWriter>(&system, "trajectory.trr");
    /// ```
    fn new(system: &System, filename: impl AsRef<Path>) -> Result<Self, WriteTrajError>
    where
        Self: Sized;

    /// Write the current state of the system into an open trajectory file.
    fn write_frame(&mut self) -> Result<(), WriteTrajError>;
}

/// Any structure implementing `TrajGroupWrite` can be used to write a trajectory file.
pub trait TrajGroupWrite {
    /// Open a new xdr file for writing and associate a specific group from a specific system with it.
    ///
    /// ## Example
    /// Using `TrajGroupWrite::new` in a generic function.
    /// ```no_run
    /// use groan_rs::prelude::*;
    /// use std::path::Path;
    ///
    /// // this function can write any trajectory file implementing `TrajGroupWrite` trait
    /// fn write_trajectory_group<Writer>(system: &System, file: impl AsRef<Path>)
    ///     where Writer: TrajGroupWrite
    /// {
    ///     // open the trajectory file for writing
    ///     let mut writer = Writer::new(system, "Protein", file).unwrap();
    ///
    ///     // write frame into the trajectory file
    ///     writer.write_frame().unwrap();
    /// }
    ///
    /// // `write_trajectory_group` can be then called to write the group `Protein` into different trajectory files
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("Protein", "@protein").unwrap();
    /// write_trajectory_group::<XtcGroupWriter>(&system, "trajectory.xtc");
    /// write_trajectory_group::<TrrGroupWriter>(&system, "trajectory.trr");
    /// ```
    fn new(
        system: &System,
        group: &str,
        filename: impl AsRef<Path>,
    ) -> Result<Self, WriteTrajError>
    where
        Self: Sized;

    /// Write the current state of the specified group into an open trajectory file.
    fn write_frame(&mut self) -> Result<(), WriteTrajError>;
}
