// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Structures and functions for efficient trajectory concatenation.

use std::marker::PhantomData;
use std::path::Path;

use crate::{errors::ReadTrajError, io::traj_io::TrajRead, system::System};

use crate::io::traj_io::{
    FrameData, FrameDataTime, TrajRangeRead, TrajReadOpen, TrajStepRead, TrajStepTimeRead,
};

use super::traj_io::{TrajFile, TrajReader};

pub struct TrajCat<'a, R: TrajRead<'a>> {
    traj_readers: Vec<R>,
    curr_reader: usize,
    last_time: Option<f32>,
    last_skipped_time: Option<f32>,
    _phantom: &'a PhantomData<R>,
}

impl<'a, R> TrajFile for TrajCat<'a, R> where R: TrajRead<'a> {}

pub struct TrajCatFrame<'a, R>
where
    R: TrajRead<'a>,
{
    frame: R::FrameData,
    _phantom: &'a PhantomData<R>,
}

impl<'a, R> FrameData for TrajCatFrame<'a, R>
where
    R: TrajRead<'a>,
    R::FrameData: FrameDataTime,
{
    type TrajFile = TrajCat<'a, R>;

    fn from_frame(
        trajcat: &mut Self::TrajFile,
        system: &System,
    ) -> Option<Result<Self, ReadTrajError>> {
        let frame = match trajcat.traj_readers.get_mut(trajcat.curr_reader) {
            // all trajectories read
            None => None,
            Some(trajectory) => {
                match R::FrameData::from_frame(trajectory.get_file_handle(), system) {
                    // no frame left -> move to next trajectory
                    None => {
                        trajcat.curr_reader += 1;
                        // save the last time of the current trajectory
                        trajcat.last_time = Some(system.get_simulation_time());
                        Self::from_frame(trajcat, system)
                    }
                    // error while reading
                    Some(Err(e)) => Some(Err(e)),
                    // normal reading
                    Some(Ok(x)) => Some(Ok(TrajCatFrame {
                        frame: x,
                        _phantom: &PhantomData,
                    })),
                }
            }
        };

        // if we are at the boundary between two trajectory files,
        // check for matching times between the first and last frame
        match (frame, trajcat.last_time, trajcat.last_skipped_time) {
            (Some(Err(e)), _, _) => Some(Err(e)),
            (None, _, _) => None,
            (x, None, None) => x,
            (Some(Ok(curr)), Some(last), None) => {
                // skip this frame, if the time matches the time of the last read frame
                if last == curr.get_time() {
                    Self::from_frame(trajcat, system)
                } else {
                    trajcat.last_time = None;
                    Some(Ok(curr))
                }
            }
            (Some(Ok(curr)), None, Some(last)) => {
                // skip this frame, if the time matches the time of the last skipped frame
                if last == curr.get_time() {
                    Self::from_frame(trajcat, system)
                } else {
                    trajcat.last_time = None;
                    Some(Ok(curr))
                }
            }
            (Some(Ok(curr)), Some(last_read), Some(last_skipped)) => {
                // skip this frame, if the time matches the time of the last read frame or the last skipped frame
                if last_read == curr.get_time() || last_skipped == curr.get_time() {
                    Self::from_frame(trajcat, system)
                } else {
                    trajcat.last_time = None;
                    Some(Ok(curr))
                }
            }
        }
    }

    fn update_system(self, system: &mut System) {
        R::FrameData::update_system(self.frame, system)
    }
}

impl<'a, R> FrameDataTime for TrajCatFrame<'a, R>
where
    R: TrajRead<'a>,
    R::FrameData: FrameDataTime,
{
    fn get_time(&self) -> f32 {
        self.frame.get_time()
    }
}

pub struct TrajConcatenator<'a, R: TrajRead<'a>> {
    // safety: only this system pointer can be used!
    system: *mut System,
    traj_cat: TrajCat<'a, R>,
    _phantom: &'a PhantomData<R>,
}

impl<'a, R> TrajRead<'a> for TrajConcatenator<'a, R>
where
    R: TrajRead<'a>,
    R::FrameData: FrameDataTime,
{
    type FrameData = TrajCatFrame<'a, R>;

    fn get_system(&mut self) -> *mut System {
        self.system
    }

    fn get_file_handle(
        &mut self,
    ) -> &mut <<Self as TrajRead<'a>>::FrameData as FrameData>::TrajFile {
        &mut self.traj_cat
    }
}

impl<'a, R> TrajRangeRead<'a> for TrajConcatenator<'a, R>
where
    R: TrajRangeRead<'a>,
    R::FrameData: FrameDataTime,
{
    fn jump_to_start(&mut self, start_time: f32) -> Result<(), ReadTrajError> {
        // jump through all trajectories until the starting time is found
        for traj in self.traj_cat.traj_readers.iter_mut() {
            match traj.jump_to_start(start_time) {
                Ok(_) => return Ok(()),
                // starting time not found in the trajectory, continue with the next trajectory
                Err(_) => continue,
            }
        }

        // if starting time is not found in any trajectory, return an error
        Err(ReadTrajError::StartNotFound(start_time.to_string()))
    }
}

impl<'a, R> TrajStepRead<'a> for TrajConcatenator<'a, R>
where
    R: TrajStepTimeRead<'a>,
    R::FrameData: FrameDataTime,
{
    fn skip_frame(&mut self) -> Result<bool, ReadTrajError> {
        let traj = match self
            .traj_cat
            .traj_readers
            .get_mut(self.traj_cat.curr_reader)
        {
            // if all trajectories are read, end the iteration
            None => return Ok(false),
            Some(x) => x,
        };

        match (
            traj.skip_frame_time(),
            self.traj_cat.last_skipped_time,
            self.traj_cat.last_time,
        ) {
            // if the frame could not be skipped (EOF), move to another trajectory and attempt to skip a frame
            (Ok(None), _, _) => {
                // we are at the file boundary; store the simulation time from the last frame
                unsafe {
                    self.traj_cat.last_time = Some((*self.get_system()).get_simulation_time());
                }
                self.traj_cat.curr_reader += 1;
                self.skip_frame()
            }

            // if an error is returned, propagate it
            (Err(e), _, _) => Err(e),

            // if the frame was skipped and no `last_skipped_frame` is stored, just store
            // the time of the frame
            (Ok(Some(time)), None, None) => {
                self.traj_cat.last_skipped_time = Some(time);
                Ok(true)
            }

            (Ok(Some(time)), None, Some(last)) => {
                // the frame that was skipped has the same time as the last READ frame
                // skip one more frame
                if time == last {
                    self.skip_frame()
                } else {
                    self.traj_cat.last_skipped_time = Some(time);
                    Ok(true)
                }
            }

            (Ok(Some(time)), Some(last), None) => {
                // the frame was skipped but it had the same time as the last skipped frame
                // skip one more frame
                if time == last {
                    self.skip_frame()
                } else {
                    self.traj_cat.last_skipped_time = Some(time);
                    Ok(true)
                }
            }

            (Ok(Some(time)), Some(last_skipped), Some(last_read)) => {
                // the frame was skipped but it has the same time as the last skipped or read frame
                // skip one more frame
                if time == last_skipped || time == last_read {
                    self.skip_frame()
                } else {
                    self.traj_cat.last_skipped_time = Some(time);
                    Ok(true)
                }
            }
        }
    }
}

impl<'a, R> TrajConcatenator<'a, R>
where
    R: TrajRead<'a>,
    R::FrameData: FrameDataTime,
{
    pub(super) fn new(mut readers: Vec<R>) -> TrajReader<'a, Self> {
        if readers.is_empty() {
            panic!("FATAL GROAN ERROR | TrajConcatenator::new | Vector of trajectory readers should not be empty.");
        }

        let system = readers
            .get_mut(0)
            .expect("FATAL GROAN ERROR | TrajConcatenator::new | The first reader should exist.")
            .get_system();

        let traj_cat = TrajCat {
            traj_readers: readers,
            curr_reader: 0,
            last_time: None,
            last_skipped_time: None,
            _phantom: &PhantomData,
        };

        let concatenator = TrajConcatenator {
            system,
            traj_cat,
            _phantom: &PhantomData,
        };

        TrajReader::wrap_traj(concatenator)
    }
}

/***************************************/
/*        Generic System methods       */
/***************************************/

impl System {
    /// Iterate through multiple trajectory files.
    /// Any duplicate frames at the boundaries of the trajectories are skipped.
    ///
    /// ## Examples
    /// Iterate across several xtc trajectories.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// let trajectories = vec!["md0001.xtc", "md0002.xtc", "md0003.xtc"];
    ///
    /// for frame in system.traj_cat_iter::<XtcReader>(&trajectories).unwrap() {
    ///     let frame = frame.unwrap();
    ///    
    ///     // perform some operation with the frame
    /// }
    /// ```
    /// Iterate across frames in a specified time interval (100 to 300 ns) that may cross several trajectories.
    /// In the specified time interval, read only every 5th frame.
    /// Also print progress of the trajectory reading.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// let trajectories = vec!["md0001.xtc", "md0002.xtc", "md0003.xtc"];
    ///
    /// for frame in system
    ///     .traj_cat_iter::<XtcReader>(&trajectories)
    ///     .unwrap()
    ///     .print_progress(ProgressPrinter::default())
    ///     .with_range(100_000.0, 300_000.0)
    ///     .unwrap()
    ///     .with_step(5)
    ///     .unwrap()
    /// {
    ///     let frame = frame.unwrap();
    ///     // perform some operation with the frame
    /// }
    /// ```
    ///
    /// ## Notes
    /// - The trajectory files are not reordered in any way and are read in the order in which they are provided.
    /// - It is quite common that the last frame of a simulation trajectory is the same frame as the first frame
    /// of the following simulation trajectory. When concatenating trajectories using this function,
    /// any duplicate frames, i.e., frames with matching simulation time, at the trajectory boundaries are removed.
    /// - ONLY duplicate frames at the trajectory boundaries are removed, all other frames are kept and iterated through.
    /// - For instance, the following three trajectory files (each composed of 4 frames with specified time)
    /// will be concatenated into one file like this:
    /// ```text
    ///       trajectory 1            trajectory 2              trajectory 3
    /// [ 0 , 100 , 200 , 300 ] [ 300 , 400 , 500 , 600 ] [ 100 , 200 , 300 , 400 ]
    /// ----->
    ///                   concatenated trajectory
    /// [ 0 , 100 , 200 , 300 , 400 , 500 , 600 , 100 , 200 , 300 , 400 ]
    /// ```
    /// - All the concatenated trajectory files must be of the same type, i.e. it is not possible to concatenate
    /// e.g. `xtc` with `trr` files.
    /// - Only trajectories which contain information about the simulation time can be concatenated using this function.
    pub fn traj_cat_iter<'a, Read>(
        &'a mut self,
        filenames: &[impl AsRef<Path>],
    ) -> Result<TrajReader<'a, TrajConcatenator<'a, Read>>, ReadTrajError>
    where
        Read: TrajReadOpen<'a>,
        Read::FrameData: FrameDataTime,
    {
        if filenames.is_empty() {
            return Err(ReadTrajError::CatNoTrajectories);
        }

        let system = self as *mut System;

        let readers: Result<Vec<Read>, _> = filenames
            .iter()
            .map(|name| unsafe { Read::new(&mut *system, name) })
            .collect();

        readers.map(TrajConcatenator::new)
    }
}

/**************************/
/*       UNIT TESTS       */
/**************************/

#[cfg(test)]
mod tests {
    use std::fs::File;

    use float_cmp::assert_approx_eq;

    use crate::io::traj_io::TrajMasterRead;
    use crate::io::trr_io::TrrReader;
    use crate::io::xtc_io::XtcReader;
    use crate::progress::ProgressPrinter;
    use crate::test_utilities::*;

    use self::utilities::*;

    use super::*;

    #[test]
    fn cat_xtc_simple() {
        let mut system_single = System::from_file("test_files/example.gro").unwrap();
        let mut system_cat = System::from_file("test_files/example.gro").unwrap();

        let traj_single = system_single
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap();
        let traj_cat = system_cat
            .traj_cat_iter::<XtcReader>(&[
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
                frame_single.get_box_as_ref().unwrap(),
                frame_cat.get_box_as_ref().unwrap(),
            );

            for (atom_single, atom_cat) in frame_single.atoms_iter().zip(frame_cat.atoms_iter()) {
                compare_atoms(atom_single, atom_cat);
            }
        }
    }

    #[test]
    fn cat_xtc_with_ranges() {
        let mut system_single = System::from_file("test_files/example.gro").unwrap();
        let mut system_cat = System::from_file("test_files/example.gro").unwrap();

        let ranges = vec![(0.0, 570.0), (320.0, f32::MAX), (320.0, 500.0)];

        for (start, end) in ranges {
            let traj_single = system_single
                .xtc_iter("test_files/short_trajectory.xtc")
                .unwrap()
                .with_range(start, end)
                .unwrap();
            let traj_cat = system_cat
                .traj_cat_iter::<XtcReader>(&[
                    "test_files/split/traj1.xtc",
                    "test_files/split/traj2.xtc",
                    "test_files/split/traj3.xtc",
                    "test_files/split/traj4.xtc",
                    "test_files/split/traj5.xtc",
                    "test_files/split/traj6.xtc",
                ])
                .unwrap()
                .with_range(start, end)
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
                    frame_single.get_box_as_ref().unwrap(),
                    frame_cat.get_box_as_ref().unwrap(),
                );

                for (atom_single, atom_cat) in frame_single.atoms_iter().zip(frame_cat.atoms_iter())
                {
                    compare_atoms(atom_single, atom_cat);
                }
            }
        }
    }

    #[test]
    fn cat_xtc_steps() {
        for step in 2..=11 {
            let mut system_single = System::from_file("test_files/example.gro").unwrap();
            let mut system_cat = System::from_file("test_files/example.gro").unwrap();

            let traj_single = system_single
                .xtc_iter("test_files/short_trajectory.xtc")
                .unwrap()
                .with_step(step)
                .unwrap();
            let traj_cat = system_cat
                .traj_cat_iter::<XtcReader>(&[
                    "test_files/split/traj1.xtc",
                    "test_files/split/traj2.xtc",
                    "test_files/split/traj3.xtc",
                    "test_files/split/traj4.xtc",
                    "test_files/split/traj5.xtc",
                    "test_files/split/traj6.xtc",
                ])
                .unwrap()
                .with_step(step)
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
                    frame_single.get_box_as_ref().unwrap(),
                    frame_cat.get_box_as_ref().unwrap(),
                );

                for (atom_single, atom_cat) in frame_single.atoms_iter().zip(frame_cat.atoms_iter())
                {
                    compare_atoms(atom_single, atom_cat);
                }
            }
        }
    }

    #[test]
    fn cat_xtc_steps_with_ranges() {
        let ranges = vec![(0.0, 570.0), (320.0, f32::MAX), (220.0, 800.0)];

        for step in 2..=11 {
            for (start, end) in &ranges {
                let mut system_single = System::from_file("test_files/example.gro").unwrap();
                let mut system_cat = System::from_file("test_files/example.gro").unwrap();

                let traj_single = system_single
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_step(step)
                    .unwrap()
                    .with_range(*start, *end)
                    .unwrap();
                let traj_cat = system_cat
                    .traj_cat_iter::<XtcReader>(&[
                        "test_files/split/traj1.xtc",
                        "test_files/split/traj2.xtc",
                        "test_files/split/traj3.xtc",
                        "test_files/split/traj4.xtc",
                        "test_files/split/traj5.xtc",
                        "test_files/split/traj6.xtc",
                    ])
                    .unwrap()
                    .with_step(step)
                    .unwrap()
                    .with_range(*start, *end)
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
                        frame_single.get_box_as_ref().unwrap(),
                        frame_cat.get_box_as_ref().unwrap(),
                    );

                    for (atom_single, atom_cat) in
                        frame_single.atoms_iter().zip(frame_cat.atoms_iter())
                    {
                        compare_atoms(atom_single, atom_cat);
                    }
                }
            }
        }
    }

    #[test]
    fn cat_trr_simple() {
        let mut system_single = System::from_file("test_files/example.gro").unwrap();
        let mut system_cat = System::from_file("test_files/example.gro").unwrap();

        let traj_single = system_single
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap();
        let traj_cat = system_cat
            .traj_cat_iter::<TrrReader>(&[
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
                frame_single.get_box_as_ref().unwrap(),
                frame_cat.get_box_as_ref().unwrap(),
            );

            for (atom_single, atom_cat) in frame_single.atoms_iter().zip(frame_cat.atoms_iter()) {
                compare_atoms(atom_single, atom_cat);
            }
        }
    }

    #[test]
    fn cat_trr_with_ranges() {
        let mut system_single = System::from_file("test_files/example.gro").unwrap();
        let mut system_cat = System::from_file("test_files/example.gro").unwrap();

        let ranges = vec![(0.0, 400.0), (250.0, f32::MAX), (250.0, 400.0)];

        for (start, end) in ranges {
            let traj_single = system_single
                .trr_iter("test_files/short_trajectory.trr")
                .unwrap()
                .with_range(start, end)
                .unwrap();
            let traj_cat = system_cat
                .traj_cat_iter::<TrrReader>(&[
                    "test_files/split/traj1.trr",
                    "test_files/split/traj2.trr",
                    "test_files/split/traj3.trr",
                    "test_files/split/traj4.trr",
                    "test_files/split/traj5.trr",
                    "test_files/split/traj6.trr",
                ])
                .unwrap()
                .with_range(start, end)
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
                    frame_single.get_box_as_ref().unwrap(),
                    frame_cat.get_box_as_ref().unwrap(),
                );

                for (atom_single, atom_cat) in frame_single.atoms_iter().zip(frame_cat.atoms_iter())
                {
                    compare_atoms(atom_single, atom_cat);
                }
            }
        }
    }

    #[test]
    fn cat_trr_steps() {
        for step in 2..=11 {
            let mut system_single = System::from_file("test_files/example.gro").unwrap();
            let mut system_cat = System::from_file("test_files/example.gro").unwrap();

            let traj_single = system_single
                .trr_iter("test_files/short_trajectory.trr")
                .unwrap()
                .with_step(step)
                .unwrap();
            let traj_cat = system_cat
                .traj_cat_iter::<TrrReader>(&[
                    "test_files/split/traj1.trr",
                    "test_files/split/traj2.trr",
                    "test_files/split/traj3.trr",
                    "test_files/split/traj4.trr",
                    "test_files/split/traj5.trr",
                    "test_files/split/traj6.trr",
                ])
                .unwrap()
                .with_step(step)
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
                    frame_single.get_box_as_ref().unwrap(),
                    frame_cat.get_box_as_ref().unwrap(),
                );

                for (atom_single, atom_cat) in frame_single.atoms_iter().zip(frame_cat.atoms_iter())
                {
                    compare_atoms(atom_single, atom_cat);
                }
            }
        }
    }

    #[test]
    fn cat_trr_steps_with_ranges() {
        let ranges = vec![(0.0, 400.0), (250.0, f32::MAX), (250.0, 400.0)];

        for step in 2..=11 {
            for (start, end) in &ranges {
                let mut system_single = System::from_file("test_files/example.gro").unwrap();
                let mut system_cat = System::from_file("test_files/example.gro").unwrap();

                let traj_single = system_single
                    .trr_iter("test_files/short_trajectory.trr")
                    .unwrap()
                    .with_step(step)
                    .unwrap()
                    .with_range(*start, *end)
                    .unwrap();
                let traj_cat = system_cat
                    .traj_cat_iter::<TrrReader>(&[
                        "test_files/split/traj1.trr",
                        "test_files/split/traj2.trr",
                        "test_files/split/traj3.trr",
                        "test_files/split/traj4.trr",
                        "test_files/split/traj5.trr",
                        "test_files/split/traj6.trr",
                    ])
                    .unwrap()
                    .with_step(step)
                    .unwrap()
                    .with_range(*start, *end)
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
                        frame_single.get_box_as_ref().unwrap(),
                        frame_cat.get_box_as_ref().unwrap(),
                    );

                    for (atom_single, atom_cat) in
                        frame_single.atoms_iter().zip(frame_cat.atoms_iter())
                    {
                        compare_atoms(atom_single, atom_cat);
                    }
                }
            }
        }
    }

    #[test]
    fn cat_traj_print_progress() {
        let mut system_single = System::from_file("test_files/example.gro").unwrap();
        let mut system_cat = System::from_file("test_files/example.gro").unwrap();

        let output_single = File::create("tmp_cat_traj_print_progress_single.txt").unwrap();
        let output_cat = File::create("tmp_cat_traj_print_progress_cat.txt").unwrap();

        let single_printer = ProgressPrinter::new()
            .with_print_freq(3)
            .with_output(Box::from(output_single))
            .with_colored(false);
        let cat_printer = ProgressPrinter::new()
            .with_print_freq(3)
            .with_output(Box::from(output_cat))
            .with_colored(false);

        let traj_single = system_single
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .print_progress(single_printer);
        let traj_cat = system_cat
            .traj_cat_iter::<XtcReader>(&[
                "test_files/split/traj1.xtc",
                "test_files/split/traj2.xtc",
                "test_files/split/traj3.xtc",
                "test_files/split/traj4.xtc",
                "test_files/split/traj5.xtc",
                "test_files/split/traj6.xtc",
            ])
            .unwrap()
            .print_progress(cat_printer);

        for frame in traj_single {
            frame.unwrap();
        }

        for frame in traj_cat {
            frame.unwrap();
        }

        let mut result = File::open("tmp_cat_traj_print_progress_cat.txt").unwrap();
        let mut expected = File::open("tmp_cat_traj_print_progress_single.txt").unwrap();
        assert!(file_diff::diff_files(&mut result, &mut expected));

        std::fs::remove_file("tmp_cat_traj_print_progress_cat.txt").unwrap();
        std::fs::remove_file("tmp_cat_traj_print_progress_single.txt").unwrap();
    }

    #[test]
    fn cat_traj_duplicates() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let times = [
            0.0, 100.0, 200.0, 300.0, 400.0, 0.0, 100.0, 200.0, 300.0, 400.0, 900.0, 1000.0, 0.0,
            100.0, 200.0,
        ];

        for (i, frame) in system
            .traj_cat_iter::<XtcReader>(&[
                "test_files/split/traj1.xtc",
                "test_files/split/traj2.xtc",
                "test_files/split/traj3.xtc",
                "test_files/split/traj1.xtc",
                "test_files/split/traj3.xtc",
                "test_files/split/traj6.xtc",
                "test_files/split/traj1.xtc",
                "test_files/split/traj2.xtc",
            ])
            .unwrap()
            .enumerate()
        {
            let frame = frame.unwrap();
            assert_approx_eq!(f32, frame.get_simulation_time(), times[i]);
        }
    }

    #[test]
    fn cat_traj_empty() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let empty: Vec<&str> = vec![];

        match system.traj_cat_iter::<XtcReader>(&empty) {
            Ok(_) => panic!("Function should have failed."),
            Err(e) => assert_eq!(ReadTrajError::CatNoTrajectories, e),
        }
    }

    #[test]
    fn cat_traj_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let empty: Vec<&str> = vec![
            "test_files/split/traj1.trr",
            "test_files/split/traj_nonexistent.trr",
            "test_files/split/traj3.trr",
        ];

        match system.traj_cat_iter::<TrrReader>(&empty) {
            Ok(_) => panic!("Function should have failed."),
            Err(ReadTrajError::FileNotFound(file)) => assert_eq!(
                file.to_str().unwrap(),
                "test_files/split/traj_nonexistent.trr"
            ),
            Err(e) => panic!("Incorrect error type returned `{}`", e),
        }
    }
}
