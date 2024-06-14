// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Basic utilities for parallel implementations of functions.

use std::{ops::Add, path::Path};

use crate::{
    io::traj_io::{
        FrameDataTime, TrajMasterRead, TrajRangeRead, TrajRead, TrajReadOpen, TrajStepRead,
    },
    progress::ProgressPrinter,
    structures::container::AtomContainer,
    system::System,
};

impl System {
    /// This method performs embarrassingly parallel iteration over a trajectory (xtc or trr) file,
    /// using the MapReduce parallel scheme.
    ///
    /// The trajectory is divided into `n_thread` threads, each analyzing a portion of the data independently (`map` phase).
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
    /// - `start_time`: The starting time for the iteration.
    /// - `end_time`: The ending time for the iteration.
    /// - `step`: The step interval for frame analysis (i.e., analyze every `step`th frame).
    /// - `progress_printer`: A specification for printing the progress of trajectory reading, used only by the master thread.
    ///
    /// ## Example
    /// Calculating the average number of phosphorus atoms above and below the global membrane center of geometry.
    /// ```no_run
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
    ///
    /// ```
    ///
    /// ## Notes
    /// - The order of iteration through trajectory frames is completely undefined!
    /// - The original `System` structure is not modified and remains in the same state as at the start of the iteration.
    /// That is different from the standard (serial) iteration over trajectories.
    /// - If a single thread encounters an error during the iteration, the entire function returns an error.
    /// However, this error is propagated only after all the other threads finish their work.
    pub fn traj_iter_map_reduce<'a, Reader, Data, Error>(
        &self,
        trajectory_file: impl AsRef<Path> + Send + Clone,
        n_threads: usize,
        body: impl Fn(&System, &mut Data) -> Result<(), Error> + Send + Clone,
        start_time: Option<f32>,
        end_time: Option<f32>,
        step: Option<usize>,
        progress_printer: Option<ProgressPrinter>,
    ) -> Result<Data, Box<dyn std::error::Error + Send + Sync>>
    where
        Reader: TrajReadOpen<'a> + TrajRangeRead<'a> + TrajStepRead<'a> + TrajRead<'a> + 'a,
        <Reader as TrajRead<'a>>::FrameData: FrameDataTime,
        Data: Add<Output = Data> + Default + Send,
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
                    let mut data = Data::default();

                    let handle = s.spawn(
                        move || -> Result<Data, Box<dyn std::error::Error + Send + Sync>> {
                            system_clone.thread_iter::<Reader, Data, Error>(
                                &mut data,
                                filename,
                                n,
                                n_threads,
                                body_clone,
                                start_time,
                                end_time,
                                step,
                                progress_clone,
                            )?;

                            Ok(data)
                        },
                    );
                    handles.push(handle);
                }

                let mut data = Vec::new();
                for handle in handles {
                    data.push(handle.join().expect(
                        "FATAL GROAN ERROR | System::traj_iter_map_reduce | A thread panicked!",
                    )?)
                }

                // reduce data
                // TODO! implement parallel reduction
                Ok(data.into_iter().fold(Data::default(), |acc, x| acc + x))
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
    ) -> Result<(), Box<dyn std::error::Error + Send + Sync>>
    where
        Reader: TrajReadOpen<'a> + TrajRangeRead<'a> + TrajStepRead<'a> + TrajRead<'a> + 'a,
        <Reader as TrajRead<'a>>::FrameData: FrameDataTime,
        Data: Add<Output = Data> + Default,
        Error: std::error::Error + Send + Sync + 'static,
    {
        let start = start_time.unwrap_or(0.0);
        let end = end_time.unwrap_or(f32::MAX);
        let step = step.unwrap_or(1);

        // prepare the iterator
        // TODO! remove this unsafe
        let mut iterator = unsafe { (*(&mut self as *mut Self)).traj_iter::<Reader>(&filename)? };

        // associate the progress printer with the iterator only in the master thread
        if thread_number == 0 {
            if let Some(printer) = progress_printer {
                iterator = iterator.print_progress(printer.clone());
            }
        }

        // find the start of the reading
        let mut iterator = iterator.with_range(start, end)?;

        // prepare the iterator for reading (skip N frames)
        for _ in 0..(thread_number * step) {
            match iterator.next() {
                Some(Ok(_)) => (),
                Some(Err(e)) => return Err(Box::from(e)),
                // if the iterator has nothing to read, then just finish with Ok
                None => return Ok(()),
            }
        }

        // the iterator should step by N * n_threads frames
        let iterator = iterator.with_step(step * n_threads)?;

        for frame in iterator {
            let frame = frame?;
            body(&frame, data)?;
        }

        Ok(())
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
}
