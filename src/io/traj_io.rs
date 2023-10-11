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
fn sanity_check_timerange(start_time: f32, end_time: f32) -> Result<(), ReadTrajError> {
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

/*********************************************/
/*  TrajFile and supported trajectory files  */
/*********************************************/

/// Any trajectory file must implement this trait.
/// Note that the exact nature of the trajectory file is not relevant,
/// but the `FrameData::from_frame` function must be able to read it.
pub trait TrajFile {}
impl TrajFile for CXdrFile {}

/*****************************/
/*  TrajRead and TrajReader  */
/*****************************/

/// Structure storing data from a single trajectory frame.
pub trait FrameData {
    type TrajFile: TrajFile;

    /// Method specifying how a frame of the trajectory should be read and stored in the `FrameData` structure.
    fn from_frame(
        traj_file: &mut Self::TrajFile,
        system: &System,
    ) -> Option<Result<Self, ReadTrajError>>
    where
        Self: Sized;

    /// Method specifying how the `System` structure should be updated based on the data in the `FrameData` structure.
    fn update_system(self, system: &mut System);
}

/// Any structure implementing `TrajRead` can be used to read a trajectory file.
pub trait TrajRead<'a> {
    type FrameData: FrameData;

    /// Method specifying how to open the trajectory file.
    /// This function does not return the structure implementing TrajRead directly,
    /// instead it should return `TrajReader` wrapping the custom structure.
    /// You can use `TrajReader::wrap_traj` to wrap your structure into `TrajReader`.
    fn new(
        system: &'a mut System,
        filename: impl AsRef<Path>,
    ) -> Result<TrajReader<'a, Self>, ReadTrajError>
    where
        Self: Sized;

    /// Method specifying how to get a mutable pointer to the `System` structure.
    /// Mutable pointer to the `System` structure must be part of the trajectory reader.
    fn get_system(&mut self) -> *mut System;

    /// Method specifying how to get a mutable handle to the file containing the trajectory.
    /// The exact nature of the 'file handle' is actually not relevant, but it must be readable by
    /// the corresponding `FrameData::from_file` method.
    fn get_file_handle(
        &mut self,
    ) -> &mut <<Self as TrajRead<'a>>::FrameData as FrameData>::TrajFile;
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
    /// Wrap trajectory reader implementing `TrajRead` into `TrajReader` structure.
    /// This should only be used in the implementations of the `TrajRead::new` methods.
    pub fn wrap_traj(traj_reader: R) -> TrajReader<'a, R> {
        TrajReader {
            traj_reader,
            _phantom: &PhantomData,
        }
    }
}

/// Iterate the `TrajReader`.
impl<'a, R: TrajRead<'a>> Iterator for TrajReader<'a, R> {
    type Item = Result<&'a mut System, ReadTrajError>;

    /// Read the next frame in the trajectory file and update the `System` structure.
    ///
    /// ## Returns
    /// - `Some(Ok(&mut System))` if the frame has been succesfully read.
    /// - `Some(Err(ReadTrajError))` if the frame could not be read.
    /// - `None` if the end of the trajectory file has been reached.
    fn next(&mut self) -> Option<Self::Item> {
        unsafe {
            let system = self.traj_reader.get_system();

            match R::FrameData::from_frame(self.traj_reader.get_file_handle(), &*system) {
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
    /// Convert `TrajReader` into `TrajRangeReader` structure iterating only through a part of the
    /// trajectory specified using the provided time range.
    ///
    /// ## Details
    /// Depending on the specific implementation, the iteration using `with_range` can be very efficient.
    /// Note however that calling this function always involves some initial overhead as it needs to
    /// locate the starting point of the iteration.
    ///
    /// For the `xtc` and `trr` files, frames with time lower than `start_time` are skipped over, i.e., the
    /// properties of the atoms in these frames are not read at all. The above-mentioned overhead is therefore very small.
    ///
    /// Iteration is ended at the frame with time corresponding to `end_time` or once the end of the trajectory file is reached.
    /// The range is inclusive on both ends, i.e., frames with `time = start_time` and `time = end_time` will be included in the iteration.
    ///
    /// If the frame corresponding to the `start_time` doesn't exist in the trajectory file,
    /// the iterator starts at the frame closest in time to but greater than the `start_time`.
    ///
    /// If either the `start_time` or the `end_time` is negative, it results in a `ReadTrajError::TimeRangeNegative` error.
    ///
    /// If the `start_time` is greater than the `end_time`, it results in a `ReadTrajError::InvalidTimeRange` error.
    ///
    /// If the `start_time` exceeds the time of any frame in the xtc file, it results in a `ReadTrajError::StartNotFound` error.
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

/// Any structure implementing this trait can be used to construct `TrajRangeReader` structure from `TrajReader`,
/// i.e. any structure implementing `TrajRangeRead` can be used to iterate over a part of a trajectory file
/// specified by the provided time range.
/// `TrajRangeRead` can be constructed by using the `TrajReader::with_range` method.
pub trait TrajRangeRead<'a>: TrajRead<'a> {
    /// Method specifying how to navigate the trajectory file to reach the frame with `start_time`.
    /// Ideally, this should involve completely skipping over the frames of the trajectory that should not be read.
    fn jump_to_start(&mut self, start_time: f32) -> Result<(), ReadTrajError>;
}

/// Structure for partial reading of trajectory files using time ranges.
pub struct TrajRangeReader<'a, R: TrajRead<'a>> {
    pub traj_reader: R,
    pub start_time: f32,
    pub end_time: f32,
    _phantom: &'a PhantomData<R>,
}

/// Iterate the `TrajRangeReader`.
impl<'a, R> Iterator for TrajRangeReader<'a, R>
where
    R: TrajRead<'a>,
    R::FrameData: FrameData,
{
    type Item = Result<&'a mut System, ReadTrajError>;

    /// Read the next frame in the specified range of the trajectory and update the `System` structure.
    ///
    /// ## Returns
    /// - `Some(Ok(&mut System))` if the frame has been succesfully read.
    /// - `Some(Err(ReadTrajError))` if the frame could not be read.
    /// - `None` if the end of the range of the end of the trajectory file has been reached.
    fn next(&mut self) -> Option<Self::Item> {
        unsafe {
            let system = self.traj_reader.get_system();

            match R::FrameData::from_frame(self.traj_reader.get_file_handle(), &*system) {
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

/// ## Generic methods for iterating over trajectory files.
impl System {
    /// Iterate over any trajectory file implementing a trajectory reader.
    /// A 'trajectory reader' is any structure implementing the `TrajRead` trait.
    ///
    /// ## Returns
    /// `TrajReader<TrajRead>` if the trajectory file has been successfully opened.
    /// `ReadTrajError` in case of an error.
    ///
    /// ## Example
    /// Using `traj_iter` in a generic function.
    /// ```no_run
    /// use groan_rs::prelude::*;
    /// use groan_rs::errors::ReadTrajError;
    /// use std::path::Path;
    ///
    /// fn generic_iteration<'a, Reader>(
    ///     system: &'a mut System,
    ///     filename: impl AsRef<Path>,
    /// ) -> Result<(), ReadTrajError>
    /// where
    ///     Reader: TrajRead<'a> + 'a,
    /// {
    ///     // loop through the trajectory
    ///     for raw_frame in system.traj_iter::<Reader>(&filename)? {
    ///         let frame = raw_frame?;
    ///         println!("{:?}", frame.group_get_center("all"));
    ///     }
    ///
    ///     Ok(())
    /// }
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    /// generic_iteration::<XtcReader>(&mut system, "trajectory.xtc").unwrap();
    /// generic_iteration::<TrrReader>(&mut system, "trajectory.trr").unwrap();
    /// ```
    ///
    /// ## Warning
    /// - `XtcReader` and `TrrReader` currently only support orthogonal simulation boxes!
    ///
    /// ## Notes
    /// - The `System` structure is modified while iterating through the trajectory.
    /// - `xtc` and `trr` files also have their own specific functions implementing iteration.
    /// See `System::xtc_iter()` and `System::trr_iter()`.
    pub fn traj_iter<'a, Read>(
        &'a mut self,
        filename: impl AsRef<Path>,
    ) -> Result<TrajReader<'a, Read>, ReadTrajError>
    where
        Read: TrajRead<'a>,
    {
        Read::new(self, filename)
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

/**************************/
/*       UNIT TESTS       */
/**************************/

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::*;
    use float_cmp::assert_approx_eq;

    fn compare_atoms(atom1: &Atom, atom2: &Atom) {
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

    #[test]
    fn traj_iter_xtc() {
        let mut system_xtc = System::from_file("test_files/example.gro").unwrap();
        let mut system_traj = System::from_file("test_files/example.gro").unwrap();

        for (raw1, raw2) in system_xtc
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .zip(
                system_traj
                    .traj_iter::<XtcReader>("test_files/short_trajectory.xtc")
                    .unwrap(),
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
    fn traj_iter_trr() {
        let mut system_trr = System::from_file("test_files/example.gro").unwrap();
        let mut system_traj = System::from_file("test_files/example.gro").unwrap();

        for (raw1, raw2) in system_trr
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .zip(
                system_traj
                    .traj_iter::<TrrReader>("test_files/short_trajectory.trr")
                    .unwrap(),
            )
        {

            let frame1 = raw1.unwrap();
            let frame2 = raw2.unwrap();

            for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                compare_atoms(atom1, atom2);
            }
        }
    }
}
