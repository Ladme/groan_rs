// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Basic utilities for parallel implementations of functions.

use crate::{structures::container::AtomContainer, system::general::System};

impl System {
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

    use crate::system::general::System;

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
