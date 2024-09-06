// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Basic utilities for parallel implementations of functions.

use std::{ops::Add, path::Path};

use crate::{
    io::traj_io::{
        FrameDataTime, TrajMasterRead, TrajRangeRead, TrajRead, TrajReadOpen, TrajStepRead,
    },
    progress::{ProgressPrinter, ProgressStatus},
    structures::container::AtomContainer,
    system::System,
};

impl System {
    /// This method performs embarrassingly parallel iteration over a trajectory (xtc or trr) file,
    /// using the MapReduce parallel scheme.
    ///
    /// The trajectory is divided into `n_threads` threads, each analyzing a portion of the data independently (`map` phase).
    /// The results from each thread are subsequently merged to produce a single result (`reduce` phase).
    ///
    /// ## Panics
    /// Panics if `n_threads` is set to zero.
    ///
    /// ## Generic Parameters
    /// The method uses the following generic parameters:
    /// - `Reader`: The type of trajectory reader to utilize (e.g., `XtcReader`, `TrrReader`).
    /// - `Data`: The data structure for storing analysis results, which must implement `Add` and `Default`.
    ///   The `Add` trait defines how `Data` instances are combined. `Default` defines the initial state of the data structure.
    /// - `Error`: The error type that the `body` function may return if an error occurs.
    ///
    /// ## Arguments
    /// - `trajectory_file`: The path to the trajectory file to read.
    /// - `n_threads`: The number of threads to spawn.
    /// - `body`: A function or closure to apply to each trajectory frame, storing results in the `Data` structure.
    /// - `init_data`: Initial state of the `Data` structure used during the calculation.
    /// - `start_time`: The starting time for the iteration.
    /// - `end_time`: The ending time for the iteration.
    /// - `step`: The step interval for frame analysis (i.e., analyze every `step`th frame).
    /// - `progress_printer`: A specification for printing the progress of trajectory reading, used only by the master thread.
    ///
    /// ## Example
    /// Calculating the average number of phosphorus atoms above and below the global membrane center of geometry.
    /// ```no_run
    /// #[cfg(feature = "parallel")]
    /// # {
    /// # use groan_rs::prelude::*;
    /// # use std::ops::Add;
    /// # use groan_rs::errors::GroupError;
    /// #
    /// // definition of the `Data` structure
    /// // we derive `Default` which will correspond to an empty structure
    /// #[derive(Debug, Clone, Default)]
    /// struct LeafletComposition {
    ///     // number of phosphorus atoms in the upper leaflet for each trajectory frame
    ///     upper: Vec<usize>,
    ///     // number of phosphorus atoms in the lower leaflet for each trajectory frame
    ///     lower: Vec<usize>,
    /// }
    ///
    /// // implementation of the Add operator allowing us to merge the results from different threads
    /// impl Add for LeafletComposition {
    ///     type Output = LeafletComposition;
    ///
    ///     fn add(self, rhs: Self) -> Self::Output {
    ///         let mut new = self.clone();
    ///
    ///         new.upper.extend(rhs.upper.iter());
    ///         new.lower.extend(rhs.lower.iter());
    ///
    ///         new
    ///     }
    /// }
    ///
    /// // assignment function that will be performed for each trajectory frame
    /// fn assign_lipids(
    ///     frame: &System,
    ///     composition: &mut LeafletComposition,
    /// ) -> Result<(), GroupError> {
    ///     let mut upper = 0;
    ///     let mut lower = 0;
    ///
    ///     let membrane_com = frame.group_get_center("Membrane")?;
    ///     for atom in frame.group_iter("Phosphori")? {
    ///         let z = atom
    ///             .distance_from_point(
    ///                 &membrane_com,
    ///                 Dimension::Z,
    ///                 frame.get_box_as_ref().unwrap(),
    ///             )
    ///             .unwrap();
    ///
    ///         // increase the upper leaflet counter if above membrane center
    ///         if z > 0.0 {
    ///             upper += 1;
    ///         // increase the lower leaflet counter if below membrane center
    ///         } else {
    ///             lower += 1;
    ///         }
    ///     }
    ///
    ///     composition.upper.push(upper);
    ///     composition.lower.push(lower);
    ///
    ///     Ok(())
    /// }
    ///
    /// // preparing the system to use
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("Membrane", "@membrane").unwrap();
    /// system
    ///     .group_create("Phosphori", "Membrane and name P")
    ///     .unwrap();
    ///
    /// // performing the parallel iteration using 4 threads
    /// // we are only interested in a time range of 200-500 ns
    /// // and we want to analyze every 5th trajectory frame
    /// let result = system
    ///     .traj_iter_map_reduce::<XtcReader, LeafletComposition, GroupError>(
    ///         "md.xtc",
    ///         4,
    ///         assign_lipids,
    ///         LeafletComposition::default(),
    ///         Some(200_000.0),
    ///         Some(500_000.0),
    ///         Some(5),
    ///         None,
    ///     )
    ///     .unwrap();
    ///
    /// // post-processing the result (collected from all threads)
    /// let av_upper_leaflet =
    ///     result.upper.iter().sum::<usize>() as f32 / result.upper.len() as f32;
    /// let av_lower_leaflet =
    ///     result.lower.iter().sum::<usize>() as f32 / result.lower.len() as f32;
    /// # }
    /// ```
    ///
    /// ## Notes
    /// - The order of iteration through trajectory frames is completely undefined!
    /// - The original `System` structure is not modified and remains in the same state as at the start of the iteration.
    ///   That is different from the standard (serial) iteration over trajectories.
    /// - If a single thread encounters an error during the iteration, the entire function returns an error.
    ///   However, this error is propagated only after all the other threads finish their work.
    pub fn traj_iter_map_reduce<'a, Reader, Data, Error>(
        &self,
        trajectory_file: impl AsRef<Path> + Send + Clone,
        n_threads: usize,
        body: impl Fn(&System, &mut Data) -> Result<(), Error> + Send + Clone,
        init_data: Data,
        start_time: Option<f32>,
        end_time: Option<f32>,
        step: Option<usize>,
        progress_printer: Option<ProgressPrinter>,
    ) -> Result<Data, Box<dyn std::error::Error + Send + Sync>>
    where
        Reader: TrajReadOpen<'a> + TrajRangeRead<'a> + TrajStepRead<'a> + TrajRead<'a> + 'a,
        <Reader as TrajRead<'a>>::FrameData: FrameDataTime,
        Data: Add<Output = Data> + Clone + Send,
        Error: std::error::Error + Send + Sync + 'static,
    {
        if n_threads == 0 {
            panic!("FATAL GROAN ERROR | System::traj_iter_map_reduce | Number of threads to spawn must be > 0.");
        }

        std::thread::scope(
            |s| -> Result<Data, Box<dyn std::error::Error + Send + Sync>> {
                let mut handles = Vec::new();

                for n in 0..n_threads {
                    let system_clone = self.clone();
                    let body_clone = body.clone();
                    let filename = trajectory_file.clone();
                    let progress_clone = progress_printer.clone();
                    let mut data = init_data.clone();

                    let handle = s.spawn(
                        move || -> (Result<Data, Box<dyn std::error::Error + Send + Sync>>, u64, f32) {
                            match system_clone.thread_iter::<Reader, Data, Error>(
                                &mut data,
                                filename,
                                n,
                                n_threads,
                                body_clone,
                                start_time,
                                end_time,
                                step,
                                progress_clone,
                            ) {
                                (Ok(_), step, time) => (Ok(data), step, time),
                                (Err(e), step, time) => (Err(e), step, time),
                            }
                        },
                    );
                    handles.push(handle);
                }

                // unpack data from the threads
                let mut data = Vec::new();
                let mut last_step = 0u64;
                let mut last_time = 0.0f32;
                for handle in handles {
                    let (result, step, time) = handle.join().expect(
                        "FATAL GROAN ERROR | System::traj_iter_map_reduce | A thread panicked!",
                    );

                    // get the last frame of the trajectory that was analyzed by any thread
                    if time > last_time {
                        last_step = step;
                        last_time = time;
                    }

                    match result {
                        Ok(x) => data.push(x),
                        Err(e) => {
                            // set the progress printer to FAILED if an error is detected
                            if let Some(mut progress) = progress_printer {
                                progress.set_status(ProgressStatus::Failed);
                                // print information about the frame where the trajectory reading failed
                                progress.print(0, step, time);
                            }
                            // propagate the error
                            return Err(e);
                        }
                    }
                }

                // print information about the actual last frame read, not about the last frame from the master thread
                if let Some(mut progress) = progress_printer {
                    progress.set_status(ProgressStatus::Completed);
                    // frame number can be any => printing will be done anyway because of the Completed status
                    progress.print(0, last_step, last_time);
                }

                // reduce data
                // TODO! implement parallel reduction
                Ok(data.into_iter().fold(init_data, |acc, x| acc + x))
            },
        )
    }

    /// Iterate over a part of the trajectory in a single thread.
    fn thread_iter<'a, Reader, Data, Error>(
        mut self,
        data: &mut Data,
        filename: impl AsRef<Path>,
        thread_number: usize,
        n_threads: usize,
        body: impl Fn(&System, &mut Data) -> Result<(), Error>,
        start_time: Option<f32>,
        end_time: Option<f32>,
        step: Option<usize>,
        progress_printer: Option<ProgressPrinter>,
    ) -> (
        Result<(), Box<dyn std::error::Error + Send + Sync>>,
        u64,
        f32,
    )
    where
        Reader: TrajReadOpen<'a> + TrajRangeRead<'a> + TrajStepRead<'a> + TrajRead<'a> + 'a,
        <Reader as TrajRead<'a>>::FrameData: FrameDataTime,
        Data: Add<Output = Data>,
        Error: std::error::Error + Send + Sync + 'static,
    {
        let start = start_time.unwrap_or(0.0);
        let end = end_time.unwrap_or(f32::MAX);
        let step = step.unwrap_or(1);

        // prepare the iterator
        // TODO! remove this unsafe
        let mut iterator =
            match unsafe { (*(&mut self as *mut Self)).traj_iter::<Reader>(&filename) } {
                Ok(x) => x,
                Err(e) => return (Err(Box::from(e)), 0, 0.0),
            };

        // associate the progress printer with the iterator only in the master thread
        if thread_number == 0 {
            if let Some(printer) = progress_printer {
                iterator = iterator.print_progress(printer.clone().with_newline_at_end(false));
            }
        }

        // find the start of the reading
        let mut iterator = match iterator.with_range(start, end) {
            Ok(x) => x,
            Err(e) => return (Err(Box::from(e)), 0, 0.0),
        };

        // prepare the iterator for reading (skip N frames)
        for _ in 0..(thread_number * step) {
            match iterator.next() {
                Some(Ok(_)) => (),
                Some(Err(e)) => {
                    return (
                        Err(Box::from(e)),
                        self.get_simulation_step(),
                        self.get_simulation_time(),
                    )
                }
                // if the iterator has nothing to read, then just finish with Ok
                None => return (Ok(()), 0, 0.0),
            }
        }

        // the iterator should step by N * n_threads frames
        let iterator = match iterator.with_step(step * n_threads) {
            Ok(x) => x,
            Err(e) => return (Err(Box::from(e)), 0, 0.0),
        };

        for frame in iterator {
            let frame = match frame {
                Ok(x) => x,
                Err(e) => {
                    return (
                        Err(Box::from(e)),
                        self.get_simulation_step(),
                        self.get_simulation_time(),
                    )
                }
            };

            if let Err(e) = body(&frame, data) {
                return (
                    Err(Box::from(e)),
                    frame.get_simulation_step(),
                    frame.get_simulation_time(),
                );
            }
        }

        // we return the final simulation step and time so we can print correct
        // information about the completion of the trajectory reading
        (
            Ok(()),
            self.get_simulation_step(),
            self.get_simulation_time(),
        )
    }

    /// Distribute all atoms of the system between threads.
    ///
    /// `n_atoms` must be >= `n_threads`.
    /// `n_threads` must be > 0.
    pub(crate) fn distribute_atoms(&self, n_atoms: usize, n_threads: usize) -> Vec<(usize, usize)> {
        let atoms_per_thread = n_atoms / n_threads;
        let mut extra_atoms = n_atoms % n_threads;

        (0..n_threads)
            .scan(0, move |start, _| {
                let additional_atom = if extra_atoms > 0 { 1 } else { 0 };

                extra_atoms = extra_atoms.saturating_sub(1);

                let end = *start + atoms_per_thread + additional_atom;
                let range = (*start, end);
                *start = end;

                Some(range)
            })
            .collect()
    }

    /// Distribute atoms of a group between threads.
    ///
    /// ## Panics
    /// Panics if the `group` does not exist.
    /// `n_threads` must be > 0.
    /// `n_threads` must be higher than or equal to the number of atoms in the group.
    #[allow(dead_code)]
    pub(crate) fn distribute_atoms_group(
        &self,
        group: &str,
        n_threads: usize,
    ) -> Vec<AtomContainer> {
        let group_atoms = self
            .get_groups_as_ref()
            .get(group)
            .expect("FATAL GROAN ERROR | System::distribute_atoms_group | Group should exist.");

        let n_atoms = group_atoms.get_n_atoms();
        let n_atoms_per_thread = n_atoms / n_threads;
        let mut extra_atoms = n_atoms % n_threads;

        let mut iterator = group_atoms.get_atoms().iter();
        let mut distribution = Vec::new();

        for _ in 0..n_threads {
            let mut thread = Vec::new();
            let atoms_per_this_thread = if extra_atoms > 0 {
                extra_atoms -= 1;
                n_atoms_per_thread + 1
            } else {
                n_atoms_per_thread
            };

            for _ in 0..atoms_per_this_thread {
                if let Some(index) = iterator.next() {
                    thread.push(index);
                }
            }

            distribution.push(AtomContainer::from_indices(thread, usize::MAX));
        }

        distribution
    }

    /*pub(crate) fn distribution2atoms_mut(
        &mut self,
        distribution: Vec<AtomContainer>,
    ) -> Vec<Vec<Mutex<&mut Atom>>> {
        let mut atom_distribution = Vec::new();
        let atoms = unsafe { self.get_atoms_as_ref_mut() as *mut Vec<Atom> };

        for block in distribution {
            let mut atom_pointers = Vec::new();

            for index in block.iter() {
                atom_pointers.push(Mutex::new(unsafe { (*atoms).get_unchecked_mut(index) }));
            }

            atom_distribution.push(atom_pointers)
        }

        atom_distribution
    }*/
}

#[cfg(test)]
mod tests {

    use std::{error::Error, fs::File, io::BufRead};

    use tempfile::NamedTempFile;

    use crate::{
        errors::AtomError,
        io::{trr_io::TrrReader, xtc_io::XtcReader},
    };

    use super::*;

    #[test]
    fn distribute_atoms() {
        for file_name in [
            "test_files/aa_peptide.pdb",
            "test_files/aa_membrane_peptide.gro",
            "test_files/example.gro",
        ] {
            let system = System::from_file(file_name).unwrap();

            for n_threads in 1..=128 {
                let rem = if system.get_n_atoms() % n_threads != 0 {
                    1
                } else {
                    0
                };
                let atoms_per_block = system.get_n_atoms() / n_threads;

                let distribution = system.distribute_atoms(system.get_n_atoms(), n_threads);
                assert_eq!(distribution.len(), n_threads);
                let mut sum = 0;
                for range in distribution {
                    assert!(
                        range.1 - range.0 == atoms_per_block
                            || range.1 - range.0 == atoms_per_block + rem
                    );
                    sum += range.1 - range.0;
                }

                assert_eq!(sum, system.get_n_atoms());
            }
        }
    }

    #[test]
    fn distribute_atoms_group() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system
            .group_create("Group", "name C1B C2B C3B PO4 D2A")
            .unwrap();

        let distribution = system.distribute_atoms_group("Group", 7);

        let mut sum = 0;
        for block in distribution {
            let atoms = block.iter().count();
            sum += atoms;
            assert!(atoms == 365 || atoms == 366);
        }

        assert_eq!(sum, 2560);
    }

    #[derive(Debug, Clone, Default)]
    struct Steps(Vec<u64>);

    impl Add for Steps {
        type Output = Steps;

        fn add(self, rhs: Self) -> Self::Output {
            let mut new = Steps(self.0);
            new.0.extend(rhs.0.into_iter());
            new
        }
    }

    fn frame_get_number(frame: &System, data: &mut Steps) -> Result<(), AtomError> {
        data.0.push(frame.get_simulation_step());
        Ok(())
    }

    fn frame_get_number_may_fail(frame: &System, data: &mut Steps) -> Result<(), AtomError> {
        if frame.get_simulation_step() == 35000 {
            Err(AtomError::OutOfRange(10)) // arbitrary error
        } else {
            data.0.push(frame.get_simulation_step());
            Ok(())
        }
    }

    fn run_traj_iter_single_threaded<'a, Reader>(
        system: &'a mut System,
        filename: &str,
        start: Option<f32>,
        end: Option<f32>,
        step: Option<usize>,
    ) -> Steps
    where
        Reader: TrajReadOpen<'a> + TrajRangeRead<'a> + TrajStepRead<'a> + TrajRead<'a> + 'a,
        <Reader as TrajRead<'a>>::FrameData: FrameDataTime,
    {
        let mut steps = Steps::default();

        let start = start.unwrap_or(0.0);
        let end = end.unwrap_or(f32::MAX);
        let step = step.unwrap_or(1);

        for frame in system
            .traj_iter::<Reader>(filename)
            .unwrap()
            .with_range(start, end)
            .unwrap()
            .with_step(step)
            .unwrap()
        {
            let frame = frame.unwrap();

            steps.0.push(frame.get_simulation_step());
        }

        steps.0.sort();
        steps
    }

    fn run_traj_iter_map_reduce<'a, Reader>(
        system: &'a System,
        filename: &str,
        n_threads: usize,
        start: Option<f32>,
        end: Option<f32>,
        step: Option<usize>,
        progress: Option<ProgressPrinter>,
    ) -> Steps
    where
        Reader: TrajReadOpen<'a> + TrajRangeRead<'a> + TrajStepRead<'a> + TrajRead<'a> + 'a,
        <Reader as TrajRead<'a>>::FrameData: FrameDataTime,
    {
        let mut steps = system
            .traj_iter_map_reduce::<Reader, Steps, AtomError>(
                filename,
                n_threads,
                frame_get_number,
                Steps::default(),
                start,
                end,
                step,
                progress,
            )
            .unwrap();

        steps.0.sort();
        steps
    }

    fn run_traj_iter_map_reduce_fails<'a, Reader>(
        system: &'a System,
        filename: &str,
        n_threads: usize,
        start: Option<f32>,
        end: Option<f32>,
        step: Option<usize>,
        progress: Option<ProgressPrinter>,
    ) -> Result<Steps, Box<dyn Error + Send + Sync>>
    where
        Reader: TrajReadOpen<'a> + TrajRangeRead<'a> + TrajStepRead<'a> + TrajRead<'a> + 'a,
        <Reader as TrajRead<'a>>::FrameData: FrameDataTime,
    {
        let mut steps = system.traj_iter_map_reduce::<Reader, Steps, AtomError>(
            filename,
            n_threads,
            frame_get_number_may_fail,
            Steps::default(),
            start,
            end,
            step,
            progress,
        )?;

        steps.0.sort();
        Ok(steps)
    }

    #[test]
    fn xtc_iter_map_reduce_basic() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let result_singlethreaded = run_traj_iter_single_threaded::<XtcReader>(
            &mut system,
            "test_files/short_trajectory.xtc",
            None,
            None,
            None,
        );

        for n_threads in 1..=16 {
            let result_multithreaded = run_traj_iter_map_reduce::<XtcReader>(
                &system,
                "test_files/short_trajectory.xtc",
                n_threads,
                None,
                None,
                None,
                None,
            );

            assert_ne!(result_multithreaded.0.len(), 0);
            assert_eq!(result_singlethreaded.0.len(), result_multithreaded.0.len());

            for (item1, item2) in result_singlethreaded
                .0
                .iter()
                .zip(result_multithreaded.0.iter())
            {
                assert_eq!(item1, item2);
            }
        }
    }

    #[test]
    fn xtc_iter_map_reduce_start() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let result_singlethreaded = run_traj_iter_single_threaded::<XtcReader>(
            &mut system,
            "test_files/short_trajectory.xtc",
            Some(200.0),
            None,
            None,
        );

        for n_threads in 1..=16 {
            let result_multithreaded = run_traj_iter_map_reduce::<XtcReader>(
                &system,
                "test_files/short_trajectory.xtc",
                n_threads,
                Some(200.0),
                None,
                None,
                None,
            );

            assert_ne!(result_multithreaded.0.len(), 0);
            assert_eq!(result_singlethreaded.0.len(), result_multithreaded.0.len());

            for (item1, item2) in result_singlethreaded
                .0
                .iter()
                .zip(result_multithreaded.0.iter())
            {
                assert_eq!(item1, item2);
            }
        }
    }

    #[test]
    fn xtc_iter_map_reduce_end() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let result_singlethreaded = run_traj_iter_single_threaded::<XtcReader>(
            &mut system,
            "test_files/short_trajectory.xtc",
            None,
            Some(750.0),
            None,
        );

        for n_threads in 1..=16 {
            let result_multithreaded = run_traj_iter_map_reduce::<XtcReader>(
                &system,
                "test_files/short_trajectory.xtc",
                n_threads,
                None,
                Some(750.0),
                None,
                None,
            );

            assert_ne!(result_multithreaded.0.len(), 0);
            assert_eq!(result_singlethreaded.0.len(), result_multithreaded.0.len());

            for (item1, item2) in result_singlethreaded
                .0
                .iter()
                .zip(result_multithreaded.0.iter())
            {
                assert_eq!(item1, item2);
            }
        }
    }

    #[test]
    fn xtc_iter_map_reduce_start_end() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let result_singlethreaded = run_traj_iter_single_threaded::<XtcReader>(
            &mut system,
            "test_files/short_trajectory.xtc",
            Some(200.0),
            Some(750.0),
            None,
        );

        for n_threads in 1..=16 {
            let result_multithreaded = run_traj_iter_map_reduce::<XtcReader>(
                &system,
                "test_files/short_trajectory.xtc",
                n_threads,
                Some(200.0),
                Some(750.0),
                None,
                None,
            );

            assert_ne!(result_multithreaded.0.len(), 0);
            assert_eq!(result_singlethreaded.0.len(), result_multithreaded.0.len());

            for (item1, item2) in result_singlethreaded
                .0
                .iter()
                .zip(result_multithreaded.0.iter())
            {
                assert_eq!(item1, item2);
            }
        }
    }

    #[test]
    fn xtc_iter_map_reduce_step() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        for step in [1, 2, 3, 5, 7, 20].into_iter() {
            let result_singlethreaded = run_traj_iter_single_threaded::<XtcReader>(
                &mut system,
                "test_files/short_trajectory.xtc",
                None,
                None,
                Some(step),
            );

            for n_threads in 1..=16 {
                let result_multithreaded = run_traj_iter_map_reduce::<XtcReader>(
                    &system,
                    "test_files/short_trajectory.xtc",
                    n_threads,
                    None,
                    None,
                    Some(step),
                    None,
                );

                assert_ne!(result_multithreaded.0.len(), 0);
                assert_eq!(result_singlethreaded.0.len(), result_multithreaded.0.len());

                for (item1, item2) in result_singlethreaded
                    .0
                    .iter()
                    .zip(result_multithreaded.0.iter())
                {
                    assert_eq!(item1, item2);
                }
            }
        }
    }

    #[test]
    fn xtc_iter_map_reduce_start_end_step() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let result_singlethreaded = run_traj_iter_single_threaded::<XtcReader>(
            &mut system,
            "test_files/short_trajectory.xtc",
            Some(200.0),
            Some(750.0),
            Some(3),
        );

        for n_threads in 1..=16 {
            let result_multithreaded = run_traj_iter_map_reduce::<XtcReader>(
                &system,
                "test_files/short_trajectory.xtc",
                n_threads,
                Some(200.0),
                Some(750.0),
                Some(3),
                None,
            );

            assert_ne!(result_multithreaded.0.len(), 0);
            assert_eq!(result_singlethreaded.0.len(), result_multithreaded.0.len());

            for (item1, item2) in result_singlethreaded
                .0
                .iter()
                .zip(result_multithreaded.0.iter())
            {
                assert_eq!(item1, item2);
            }
        }
    }

    #[test]
    fn xtc_iter_map_reduce_progress_print() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let result_singlethreaded = run_traj_iter_single_threaded::<XtcReader>(
            &mut system,
            "test_files/short_trajectory.xtc",
            None,
            None,
            None,
        );

        let output = File::create("xtc_iter_map_reduce_progress_print.txt").unwrap();

        let printer = ProgressPrinter::new()
            .with_output(Box::from(output))
            .with_print_freq(1)
            .with_colored(false)
            .with_terminating("\n");

        let result_multithreaded = run_traj_iter_map_reduce::<XtcReader>(
            &system,
            "test_files/short_trajectory.xtc",
            4,
            None,
            None,
            None,
            Some(printer),
        );

        assert_ne!(result_multithreaded.0.len(), 0);
        assert_eq!(result_singlethreaded.0.len(), result_multithreaded.0.len());

        for (item1, item2) in result_singlethreaded
            .0
            .iter()
            .zip(result_multithreaded.0.iter())
        {
            assert_eq!(item1, item2);
        }

        let mut result = File::open("xtc_iter_map_reduce_progress_print.txt").unwrap();
        let mut expected = File::open("test_files/progress_multithreaded_4.txt").unwrap();
        assert!(file_diff::diff_files(&mut result, &mut expected));

        std::fs::remove_file("xtc_iter_map_reduce_progress_print.txt").unwrap();
    }

    #[test]
    fn xtc_iter_map_reduce_progress_print_many_threads() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let result_singlethreaded = run_traj_iter_single_threaded::<XtcReader>(
            &mut system,
            "test_files/short_trajectory.xtc",
            None,
            None,
            None,
        );

        let output = File::create("xtc_iter_map_reduce_progress_print_many_threads.txt").unwrap();

        let printer = ProgressPrinter::new()
            .with_output(Box::from(output))
            .with_print_freq(1)
            .with_colored(false)
            .with_terminating("\n");

        let result_multithreaded = run_traj_iter_map_reduce::<XtcReader>(
            &system,
            "test_files/short_trajectory.xtc",
            21,
            None,
            None,
            None,
            Some(printer),
        );

        assert_ne!(result_multithreaded.0.len(), 0);
        assert_eq!(result_singlethreaded.0.len(), result_multithreaded.0.len());

        for (item1, item2) in result_singlethreaded
            .0
            .iter()
            .zip(result_multithreaded.0.iter())
        {
            assert_eq!(item1, item2);
        }

        let mut result = File::open("xtc_iter_map_reduce_progress_print_many_threads.txt").unwrap();
        let mut expected = File::open("test_files/progress_multithreaded_many.txt").unwrap();
        assert!(file_diff::diff_files(&mut result, &mut expected));

        std::fs::remove_file("xtc_iter_map_reduce_progress_print_many_threads.txt").unwrap();
    }

    #[test]
    fn trr_iter_map_reduce_basic() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let result_singlethreaded = run_traj_iter_single_threaded::<TrrReader>(
            &mut system,
            "test_files/short_trajectory.trr",
            None,
            None,
            None,
        );

        for n_threads in 1..=16 {
            let result_multithreaded = run_traj_iter_map_reduce::<TrrReader>(
                &system,
                "test_files/short_trajectory.trr",
                n_threads,
                None,
                None,
                None,
                None,
            );

            assert_ne!(result_multithreaded.0.len(), 0);
            assert_eq!(result_singlethreaded.0.len(), result_multithreaded.0.len());

            for (item1, item2) in result_singlethreaded
                .0
                .iter()
                .zip(result_multithreaded.0.iter())
            {
                assert_eq!(item1, item2);
            }
        }
    }

    #[test]
    fn trr_iter_map_reduce_start() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let result_singlethreaded = run_traj_iter_single_threaded::<TrrReader>(
            &mut system,
            "test_files/short_trajectory.trr",
            Some(200.0),
            None,
            None,
        );

        for n_threads in 1..=16 {
            let result_multithreaded = run_traj_iter_map_reduce::<TrrReader>(
                &system,
                "test_files/short_trajectory.trr",
                n_threads,
                Some(200.0),
                None,
                None,
                None,
            );

            assert_ne!(result_multithreaded.0.len(), 0);
            assert_eq!(result_singlethreaded.0.len(), result_multithreaded.0.len());

            for (item1, item2) in result_singlethreaded
                .0
                .iter()
                .zip(result_multithreaded.0.iter())
            {
                assert_eq!(item1, item2);
            }
        }
    }

    #[test]
    fn trr_iter_map_reduce_end() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let result_singlethreaded = run_traj_iter_single_threaded::<TrrReader>(
            &mut system,
            "test_files/short_trajectory.trr",
            None,
            Some(510.0),
            None,
        );

        for n_threads in 1..=16 {
            let result_multithreaded = run_traj_iter_map_reduce::<TrrReader>(
                &system,
                "test_files/short_trajectory.trr",
                n_threads,
                None,
                Some(510.0),
                None,
                None,
            );

            assert_ne!(result_multithreaded.0.len(), 0);
            assert_eq!(result_singlethreaded.0.len(), result_multithreaded.0.len());

            for (item1, item2) in result_singlethreaded
                .0
                .iter()
                .zip(result_multithreaded.0.iter())
            {
                assert_eq!(item1, item2);
            }
        }
    }

    #[test]
    fn trr_iter_map_reduce_start_end() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let result_singlethreaded = run_traj_iter_single_threaded::<TrrReader>(
            &mut system,
            "test_files/short_trajectory.trr",
            Some(200.0),
            Some(510.0),
            None,
        );

        for n_threads in 1..=16 {
            let result_multithreaded = run_traj_iter_map_reduce::<TrrReader>(
                &system,
                "test_files/short_trajectory.trr",
                n_threads,
                Some(200.0),
                Some(510.0),
                None,
                None,
            );

            assert_ne!(result_multithreaded.0.len(), 0);
            assert_eq!(result_singlethreaded.0.len(), result_multithreaded.0.len());

            for (item1, item2) in result_singlethreaded
                .0
                .iter()
                .zip(result_multithreaded.0.iter())
            {
                assert_eq!(item1, item2);
            }
        }
    }

    #[test]
    fn trr_iter_map_reduce_step() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        for step in [1, 2, 3, 5, 7, 20].into_iter() {
            let result_singlethreaded = run_traj_iter_single_threaded::<TrrReader>(
                &mut system,
                "test_files/short_trajectory.trr",
                None,
                None,
                Some(step),
            );

            for n_threads in 1..=16 {
                let result_multithreaded = run_traj_iter_map_reduce::<TrrReader>(
                    &system,
                    "test_files/short_trajectory.trr",
                    n_threads,
                    None,
                    None,
                    Some(step),
                    None,
                );

                assert_ne!(result_multithreaded.0.len(), 0);
                assert_eq!(result_singlethreaded.0.len(), result_multithreaded.0.len());

                for (item1, item2) in result_singlethreaded
                    .0
                    .iter()
                    .zip(result_multithreaded.0.iter())
                {
                    assert_eq!(item1, item2);
                }
            }
        }
    }

    #[test]
    fn trr_iter_map_reduce_start_end_step() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let result_singlethreaded = run_traj_iter_single_threaded::<TrrReader>(
            &mut system,
            "test_files/short_trajectory.trr",
            Some(200.0),
            Some(510.0),
            Some(3),
        );

        for n_threads in 1..=16 {
            let result_multithreaded = run_traj_iter_map_reduce::<TrrReader>(
                &system,
                "test_files/short_trajectory.trr",
                n_threads,
                Some(200.0),
                Some(510.0),
                Some(3),
                None,
            );

            assert_ne!(result_multithreaded.0.len(), 0);
            assert_eq!(result_singlethreaded.0.len(), result_multithreaded.0.len());

            for (item1, item2) in result_singlethreaded
                .0
                .iter()
                .zip(result_multithreaded.0.iter())
            {
                assert_eq!(item1, item2);
            }
        }
    }

    #[test]
    fn xtc_iter_map_reduce_fail() {
        let system = System::from_file("test_files/example.gro").unwrap();

        for n_threads in 1..=16 {
            assert!(run_traj_iter_map_reduce_fails::<XtcReader>(
                &system,
                "test_files/short_trajectory.xtc",
                n_threads,
                None,
                None,
                None,
                None,
            )
            .is_err());
        }
    }

    #[test]
    fn xtc_iter_map_reduce_fail_progress_printing() {
        let system = System::from_file("test_files/example.gro").unwrap();

        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_owned();

        let progress = ProgressPrinter::new()
            .with_output(Box::from(output))
            .with_terminating("\n")
            .with_colored(false);
        for n_threads in [1, 2, 3, 4, 8, 16] {
            assert!(run_traj_iter_map_reduce_fails::<XtcReader>(
                &system,
                "test_files/short_trajectory.xtc",
                n_threads,
                None,
                None,
                None,
                Some(progress.clone()),
            )
            .is_err());
        }

        let file = File::open(path_to_output).unwrap();
        let reader = std::io::BufReader::new(file);
        let pattern = "[ FAILED! ]   Step        35000 | Time          700 ps";

        let count = reader
            .lines()
            .filter_map(|line| line.ok())
            .filter(|line| line.contains(pattern))
            .count();

        assert_eq!(count, 6);
    }
}
