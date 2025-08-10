// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Methods for iterating over atoms of the System.

use std::collections::{HashSet, VecDeque};

use crate::structures::atom::Atom;
use crate::structures::group::Group;
use crate::structures::iterators::{
    AtomIterator, AtomPairIterator, MoleculeIterator, MutAtomIterator, MutAtomPairIterator,
    MutMoleculeIterator, OwnedAtomIterator, OwnedMutAtomIterator,
};
use crate::structures::simbox::SimBox;
use crate::system::System;

use crate::errors::{AtomError, GroupError, SelectError};

/// ## Methods for iterating over atoms of the system.
impl System {
    /// Create an iterator over a group of atoms. The atoms are immutable.
    ///
    /// ## Returns
    /// `AtomIterator` or `GroupError::NotFound` in case the group does not exist.
    ///
    /// ## Example
    /// Printing the atoms of group "Protein".
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("Protein", "resid 1 to 29").unwrap();
    ///
    /// match system.group_iter("Protein") {
    ///     Ok(iterator) => {
    ///         for atom in iterator {
    ///             println!("{:?}", atom);
    ///         }
    ///     },
    ///     Err(e) => eprintln!("{}", e),
    /// };
    /// ```
    pub fn group_iter(&self, name: &str) -> Result<AtomIterator<'_>, GroupError> {
        let group = self.groups.get(name)?;

        Ok(AtomIterator::new(
            self.get_atoms(),
            group.get_atoms(),
            self.get_box(),
        ))
    }

    /// Create an iterator over a group of atoms. The atoms are mutable.
    ///
    /// ## Returns
    /// `MutAtomIterator` or `GroupError::NotFound` in case the group does not exist.
    ///
    /// ## Example
    /// Translating the atoms of the group "Protein".
    /// Note that using `system.group_translate()` may be faster.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("Protein", "resid 1 to 29").unwrap();
    /// let simulation_box = system.get_box_copy().unwrap();
    ///
    /// match system.group_iter_mut("Protein") {
    ///     Ok(iterator) => {
    ///         for atom in iterator {
    ///             atom.translate(&[1.0, 2.0, -1.0].into(), &simulation_box);
    ///         }
    ///     },
    ///     Err(e) => eprintln!("{}", e),
    /// };
    /// ```
    pub fn group_iter_mut(&mut self, name: &str) -> Result<MutAtomIterator<'_>, GroupError> {
        let simbox = self.get_box().map(|x| x as *const SimBox);

        let group = { self.groups.get_mut(name)? as *mut Group };

        unsafe {
            Ok(MutAtomIterator::new(
                self.get_atoms_mut(),
                (*group).get_atoms(),
                simbox.map(|x| &*x),
            ))
        }
    }

    /// Create an iterator over all atoms in the system. The atoms are immutable.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// for atom in system.atoms_iter() {
    ///     println!("{:?}", atom);
    /// }
    /// ```
    ///
    /// ## Note on performance
    /// It might be slightly faster to get reference to the atoms [`System::get_atoms`](`crate::system::System::get_atoms`)
    /// and iterate through it if you do not care about the additional methods [`AtomIterator`](`crate::structures::iterators::AtomIterator`) implements.
    #[inline(always)]
    pub fn atoms_iter(&self) -> AtomIterator<'_> {
        self.group_iter("all")
            .expect("FATAL GROAN ERROR | System::atoms_iter | Default group `all` does not exist.")
    }

    /// Create an iterator over all atoms in the system. The atoms are mutable.
    ///
    /// ## Example
    /// Translating all the atoms in the system by a specified vector.
    /// Note that using [`System::atoms_translate`](`crate::system::System::atoms_translate`) might be faster.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// let simulation_box = system.get_box_copy().unwrap();
    ///
    /// for atom in system.atoms_iter_mut() {
    ///     atom.translate(&[1.0, -1.0, 2.5].into(), &simulation_box);
    /// }
    /// ```
    #[inline(always)]
    pub fn atoms_iter_mut(&mut self) -> MutAtomIterator<'_> {
        self.group_iter_mut("all").expect(
            "FATAL GROAN ERROR | System::atoms_iter_mut | Default group `all` does not exist.",
        )
    }

    /// Create an iterator over atoms that are bonded to atom with target `index`.
    /// The atoms are immutable. Atoms are indexed starting from 0.
    ///
    /// ## Returns
    /// `AtomIterator` if the index is valid. `AtomError` if the index is out of range.
    ///
    /// ## Example
    /// Calculating distances between an atom and atoms that are bonded to it.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.pdb").unwrap();
    /// system.add_bonds_from_pdb("system.pdb").unwrap();
    ///
    /// // get target atom
    /// let target_atom = match system.get_atom(15) {
    ///     Ok(atom) => atom,
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;
    ///     }
    /// };
    ///
    /// let mut distances = Vec::new();
    /// // iterate over atoms bonded to an atom indexed as 15
    /// match system.bonded_atoms_iter(15) {
    ///     Ok(iterator) => {
    ///         for bonded_atom in iterator {
    ///             distances.push(target_atom.distance(bonded_atom, Dimension::XYZ, system.get_box().unwrap()));
    ///         }
    ///         println!("{:?}", distances);
    ///     }
    ///     Err(e) => eprintln!("{}", e),
    /// }
    /// ```
    pub fn bonded_atoms_iter(&self, index: usize) -> Result<AtomIterator<'_>, AtomError> {
        let atom = self.get_atom(index)?;

        Ok(AtomIterator::new(
            self.get_atoms(),
            atom.get_bonded(),
            self.get_box(),
        ))
    }

    /// Create an iterator over atoms that are bonded to atom with target `index`.
    /// The atoms are mutable. Atoms are indexed starting from 0.
    ///
    /// ## Returns
    /// `MutAtomIterator` if the index is valid. `AtomError` if the index is out of range.
    ///
    /// ## Example
    /// Rename atoms that are bonded to target atom.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.pdb").unwrap();
    /// system.add_bonds_from_pdb("system.pdb").unwrap();
    ///
    /// // iterate over atoms bonded to an atom indexed as 15
    /// match system.bonded_atoms_iter_mut(15) {
    ///     Ok(iterator) => {
    ///         for atom in iterator {
    ///             atom.set_atom_name("ATM");
    ///         }
    ///     }
    ///     Err(e) => eprintln!("{}", e),
    /// }
    /// ```
    pub fn bonded_atoms_iter_mut(&mut self, index: usize) -> Result<MutAtomIterator<'_>, AtomError> {
        let simbox = self.get_box().map(|x| x as *const SimBox);

        unsafe {
            let atom = self.get_atom_mut(index)? as *mut Atom;

            Ok(MutAtomIterator::new(
                self.get_atoms_mut(),
                (*atom).get_bonded(),
                simbox.map(|x| &*x),
            ))
        }
    }

    /// Perform breadth-first iteration through a molecule starting from atom with `index`.
    /// The atoms are immutable. Atoms are indexed starting from 0.
    ///
    /// ## Returns
    /// `MoleculeIterator` if the index is valid. `AtomError` otherwise.
    ///
    /// ## Notes
    /// - This function performs breadth-first traversal. Given a molecule
    /// ```text
    /// 1 --- 2 --- 3
    /// |            \___
    /// |                \
    /// 4 --- 5 --- 6 --- 7
    /// |
    /// |
    /// 8 --- 9
    /// ```
    ///   and starting from atom 1, the method will traverse the atoms in this order:
    ///   1, 2, 4, 3, 5, 8, 7, 6, 9
    /// - I.e. atoms closer to the starting atom will be visited before atoms that are further way.
    pub fn molecule_iter(&self, index: usize) -> Result<MoleculeIterator<'_>, AtomError> {
        let indices = get_molecule_indices(self, index)?;
        Ok(MoleculeIterator::new(
            self.get_atoms(),
            indices,
            self.get_box(),
        ))
    }

    /// Perform breadth-first iteration through a molecule starting from atom with `index`.
    /// The atoms are mutable. Atoms are indexed starting from 0.
    ///
    /// ## Returns
    /// `MutableMoleculeIterator` if the index is valid. `AtomError` otherwise.
    ///
    /// ## Notes
    /// - This function performs breadth-first traversal. Given a molecule
    /// ```text
    /// 1 --- 2 --- 3
    /// |            \___
    /// |                \
    /// 4 --- 5 --- 6 --- 7
    /// |
    /// |
    /// 8 --- 9
    /// ```
    ///   and starting from atom 1, the method will traverse the atoms in this order:
    ///   1, 2, 4, 3, 5, 8, 7, 6, 9
    /// - I.e. atoms closer to the starting atom will be visited before atoms that are further way.
    pub fn molecule_iter_mut(&mut self, index: usize) -> Result<MutMoleculeIterator<'_>, AtomError> {
        let indices = get_molecule_indices(self, index)?;
        let simbox = self.get_box().map(|x| x as *const SimBox);

        unsafe {
            Ok(MutMoleculeIterator::new(
                self.get_atoms_mut(),
                indices,
                simbox.map(|x| &*x),
            ))
        }
    }

    /// Perform breadth-first iteration through a molecule starting from the atom with `index` and
    /// collect all bonded pairs of atoms in the molecule.
    /// The atoms are immutable. Atoms are indexed starting from 0.
    ///
    /// ## Returns
    /// `AtomPairIterator` if the index is valid. `AtomError` otherwise.
    ///
    /// ## Notes
    /// - This function performs breadth-first traveral. Given a molecule
    /// ```text
    /// 1 --- 2 --- 3
    /// |            \___
    /// |                \
    /// 4 --- 5 --- 6 --- 7
    /// |
    /// |
    /// 8 --- 9
    /// ```
    ///   and starting from the atom 1, the metod will traverse the pairs of atoms in this order:
    ///   (1, 2), (1, 4), (2, 3), (4, 5), (4, 8), **(3, 7)**, (5, 6), (8, 9), (6, 7).
    /// - Inside the pair, the atom with the lower index is always first.
    pub fn molecule_bonds_iter(&self, index: usize) -> Result<AtomPairIterator<'_>, AtomError> {
        let pairs = get_molecule_bonds(self, index)?;
        Ok(AtomPairIterator::new(self.get_atoms(), pairs))
    }

    /// Perform breadth-first iteration through a molecule starting from the atom with `index` and
    /// collect all bonded pairs of atoms in the molecule.
    /// The atoms are mutable. Atoms are indexed starting from 0.
    ///
    /// ## Returns
    /// `MutAtomPairIterator` if the index is valid. `AtomError` otherwise.
    ///
    /// ## Notes
    /// - This function performs breadth-first traveral. Given a molecule
    /// ```text
    /// 1 --- 2 --- 3
    /// |            \___
    /// |                \
    /// 4 --- 5 --- 6 --- 7
    /// |
    /// |
    /// 8 --- 9
    /// ```
    ///   and starting from the atom 1, the metod will traverse the pairs of atoms in this order:
    ///   (1, 2), (1, 4), (2, 3), (4, 5), (4, 8), **(3, 7)**, (5, 6), (8, 9), (6, 7).
    /// - Inside the pair, the atom with the lower index is always first.
    pub fn molecule_bonds_iter_mut(
        &mut self,
        index: usize,
    ) -> Result<MutAtomPairIterator<'_>, AtomError> {
        let pairs = get_molecule_bonds(self, index)?;
        Ok(MutAtomPairIterator::new(self.get_atoms_mut(), pairs))
    }

    /// Create an iterator over atoms specified using a Groan Selection Language query.
    /// This allows selecting atoms without adding a group into the system.
    ///
    /// ## Returns
    /// `OwnedAtomIterator` if the query is valid. Otherwise `SelectError`.
    ///
    /// ## Example
    /// Iterate over atoms with name `CA` or residue name `LYS`.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let system = System::from_file("system.gro").unwrap();
    ///
    /// for atom in system.selection_iter("name CA or resname LYS").unwrap() {
    ///     // perform some operation with the atom
    /// }
    /// ```
    pub fn selection_iter(&self, query: &str) -> Result<OwnedAtomIterator<'_>, SelectError> {
        let group = crate::structures::group::Group::from_query(query, self)?;

        Ok(OwnedAtomIterator::new(
            self.get_atoms(),
            group.get_atoms().to_owned(),
            self.get_box(),
        ))
    }

    /// Create a mutable iterator over atoms specified using a Groan Selection Language query.
    /// This allows selecting atoms without adding a group into the system.
    ///
    /// ## Returns
    /// `OwnedMutAtomIterator` if the query is valid. Otherwise `SelectError`.
    ///
    /// ## Example
    /// Iterate over atoms with name `CA` or residue name `LYS`.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// for atom in system.selection_iter_mut("name CA or resname LYS").unwrap() {
    ///     // perform some operation with the atom
    /// }
    /// ```
    pub fn selection_iter_mut(&mut self, query: &str) -> Result<OwnedMutAtomIterator<'_>, SelectError> {
        let group = crate::structures::group::Group::from_query(query, self)?;

        let simbox = self.get_box().map(|x| x as *const SimBox);

        unsafe {
            Ok(OwnedMutAtomIterator::new(
                self.get_atoms_mut(),
                group.get_atoms().to_owned(),
                simbox.map(|x| &*x),
            ))
        }
    }
}

/// Perform breadth-first traversal through a molecule.
///
/// ## Returns
/// Vector of indices of atoms of the molecule.
/// `AtomError` in case the index is out of range.
pub(crate) fn get_molecule_indices(system: &System, index: usize) -> Result<Vec<usize>, AtomError> {
    if index >= system.get_n_atoms() {
        return Err(AtomError::OutOfRange(index));
    }

    let mut indices = Vec::new();
    let mut queue = VecDeque::new();
    let mut visited = HashSet::new();

    queue.push_back(index);
    visited.insert(index);

    // perform BFS
    while !queue.is_empty() {
        let index = queue
            .pop_front()
            .expect("FATAL GROAN ERROR | iterators::get_molecule_indices | Attempted to dequeue element from an empty queue.");

        let atom = system.get_atom(index).expect(
            "FATAL GROAN ERROR | iterators::get_molecule_indices | Atom index does not exist.",
        );

        // add atom to list of indices
        indices.push(index);

        for bonded_index in atom.get_bonded().iter() {
            if visited.insert(bonded_index) {
                queue.push_back(bonded_index);
            }
        }
    }

    Ok(indices)
}

/// Perform breadth-first traversal through a molecule collecting all bonds.
///
/// ## Returns
/// Vector of pairs of bonded atom indices that are part of the molecule.
/// `AtomError` in case the index is out of range.
pub(crate) fn get_molecule_bonds(
    system: &System,
    index: usize,
) -> Result<Vec<(usize, usize)>, AtomError> {
    if index >= system.get_n_atoms() {
        return Err(AtomError::OutOfRange(index));
    }

    let mut pairs = Vec::new();
    let mut queue = VecDeque::new();
    let mut visited = HashSet::new();

    queue.push_back(index);
    visited.insert(index);

    // perform BFS
    while !queue.is_empty() {
        let index = queue
            .pop_front()
            .expect("FATAL GROAN ERROR | iterators::get_molecule_bonds | Attempted to dequeue element from an empty queue.");

        let atom = system.get_atom(index).expect(
            "FATAL GROAN ERROR | iterators::get_molecule_bonds | Atom index does not exist.",
        );

        for bonded_index in atom.get_bonded().iter() {
            let pair = if index < bonded_index {
                (index, bonded_index)
            } else {
                (bonded_index, index)
            };

            // O(n) operation but we can't use HashSet since it does not maintain the order of items
            if !pairs.contains(&pair) {
                pairs.push(pair);
            }

            if visited.insert(bonded_index) {
                queue.push_back(bonded_index);
            }
        }
    }

    Ok(pairs)
}

/**************************/
/*       UNIT TESTS       */
/**************************/

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        structures::{
            dimension::Dimension,
            iterators::{AtomIteratorWithBox, MutAtomIteratorWithBox},
            shape::*,
            vector3d::Vector3D,
        },
        test_utilities::utilities::compare_atoms,
    };

    use float_cmp::assert_approx_eq;

    #[test]
    fn group_iter_protein() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        for (group_atom, system_atom) in system
            .group_iter("Protein")
            .unwrap()
            .zip(system.get_atoms().iter().take(61))
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }
    }

    #[test]
    fn group_iter_membrane() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        for (group_atom, system_atom) in system
            .group_iter("Membrane")
            .unwrap()
            .zip(system.get_atoms().iter().skip(61))
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }
    }

    #[test]
    fn group_iter_ion() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        for (group_atom, system_atom) in system
            .group_iter("ION")
            .unwrap()
            .zip(system.get_atoms().iter().skip(16604))
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }
    }

    #[test]
    fn atoms_iter() {
        let system = System::from_file("test_files/example.gro").unwrap();

        let atoms = system.get_atoms();

        for (extracted_atom, system_atom) in atoms.iter().zip(system.atoms_iter()) {
            assert_eq!(
                system_atom.get_atom_number(),
                extracted_atom.get_atom_number()
            );
        }
    }

    #[test]
    fn group_iter_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let extracted = system.get_atoms_copy();

        for (group_atom, system_atom) in system
            .group_iter_mut("Protein")
            .unwrap()
            .zip(extracted.iter().take(61))
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());

            group_atom
                .translate_nopbc(&Vector3D::new(0.5, -1.1, 2.4))
                .unwrap();
            assert_approx_eq!(
                f32,
                group_atom.get_position().unwrap().x,
                system_atom.get_position().unwrap().x + 0.5
            );
            assert_approx_eq!(
                f32,
                group_atom.get_position().unwrap().y,
                system_atom.get_position().unwrap().y - 1.1
            );
            assert_approx_eq!(
                f32,
                group_atom.get_position().unwrap().z,
                system_atom.get_position().unwrap().z + 2.4
            );
        }
    }

    #[test]
    fn atoms_iter_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let extracted = system.get_atoms_copy();

        for (group_atom, system_atom) in system.atoms_iter_mut().zip(extracted.iter()) {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());

            group_atom
                .translate_nopbc(&Vector3D::new(0.5, -1.1, 2.4))
                .unwrap();
            assert_approx_eq!(
                f32,
                group_atom.get_position().unwrap().x,
                system_atom.get_position().unwrap().x + 0.5
            );
            assert_approx_eq!(
                f32,
                group_atom.get_position().unwrap().y,
                system_atom.get_position().unwrap().y - 1.1
            );
            assert_approx_eq!(
                f32,
                group_atom.get_position().unwrap().z,
                system_atom.get_position().unwrap().z + 2.4
            );
        }
    }

    #[test]
    fn filter_sphere() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let sphere_pos = system.group_get_center("Protein").unwrap();
        let sphere = Sphere::new(sphere_pos.clone(), 5.0);

        for atom in system.atoms_iter().filter_geometry(sphere.clone()) {
            assert!(
                atom.get_position().unwrap().distance(
                    &sphere_pos,
                    Dimension::XYZ,
                    system.get_box().unwrap()
                ) < 5.0
            );
        }

        let sphere_pos2 = Vector3D::new(9.0, 3.0, 1.0);
        let sphere2 = Sphere::new(sphere_pos2.clone(), 4.0);

        for atom in system
            .atoms_iter()
            .filter_geometry(sphere)
            .filter_geometry(sphere2)
        {
            assert!(
                atom.get_position().unwrap().distance(
                    &sphere_pos,
                    Dimension::XYZ,
                    system.get_box().unwrap()
                ) < 5.0
            );
            assert!(
                atom.get_position().unwrap().distance(
                    &sphere_pos2,
                    Dimension::XYZ,
                    system.get_box().unwrap()
                ) < 4.0
            );
        }
    }

    #[test]
    fn filter_sphere_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let sphere = Sphere::new(system.group_get_center("Protein").unwrap(), 1.5);
        let sphere2 = Sphere::new(system.group_get_center("Protein").unwrap(), 0.4);

        let sbox = system.get_box_copy();

        let orig_atoms = system.get_atoms_copy();

        for atom in system
            .atoms_iter_mut()
            .filter_geometry(sphere.clone())
            .filter_geometry(sphere2.clone())
        {
            atom.set_atom_name("XYZ");
        }

        for (i, atom) in system.atoms_iter().enumerate() {
            if sphere.inside(atom.get_position().unwrap(), sbox.as_ref().unwrap())
                && sphere2.inside(atom.get_position().unwrap(), sbox.as_ref().unwrap())
            {
                assert_eq!(atom.get_atom_name(), "XYZ");
            } else {
                assert_eq!(atom.get_atom_name(), orig_atoms[i].get_atom_name());
            }
        }
    }

    #[test]
    fn filter_sphere_count_water() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let sphere = Sphere::new(system.group_get_center("Protein").unwrap(), 4.5);
        let count = system
            .group_iter("W")
            .unwrap()
            .filter_geometry(sphere)
            .count();

        assert_eq!(count, 1305);
    }

    #[test]
    fn filter_sphere_count_phosphates() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.group_create("Phosphates", "name PO4").unwrap();

        let sphere = Sphere::new(system.group_get_center("Protein").unwrap(), 2.1);
        let count = system
            .group_iter("Phosphates")
            .unwrap()
            .filter_geometry(sphere)
            .count();

        assert_eq!(count, 6);
    }

    #[test]
    fn filter_xcylinder_count_water() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let cylinder = Cylinder::new(
            system.group_get_center("Protein").unwrap(),
            2.0,
            3.0,
            Dimension::X,
        );
        let count = system
            .group_iter("W")
            .unwrap()
            .filter_geometry(cylinder)
            .count();

        assert_eq!(count, 29);
    }

    #[test]
    fn filter_ycylinder_count_water() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let cylinder = Cylinder::new(
            system.group_get_center("Protein").unwrap(),
            2.0,
            3.0,
            Dimension::Y,
        );
        let count = system
            .group_iter("W")
            .unwrap()
            .filter_geometry(cylinder)
            .count();

        assert_eq!(count, 17);
    }

    #[test]
    fn filter_zcylinder_count_water() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let cylinder = Cylinder::new(
            system.group_get_center("Protein").unwrap(),
            2.0,
            3.0,
            Dimension::Z,
        );
        let count = system
            .group_iter("W")
            .unwrap()
            .filter_geometry(cylinder)
            .count();

        assert_eq!(count, 79);
    }

    #[test]
    fn filter_rectangular_count_water() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let rect = Rectangular::new(system.group_get_center("Protein").unwrap(), 2.0, 3.0, 4.0);
        let count = system
            .group_iter("W")
            .unwrap()
            .filter_geometry(rect)
            .count();

        assert_eq!(count, 92);
    }

    #[test]
    fn filter_rectangular_full_count_water() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let rect = Rectangular::new([0.0, 0.0, 0.0].into(), 100.0, 100.0, 100.0);
        let count = system
            .group_iter("W")
            .unwrap()
            .filter_geometry(rect)
            .count();

        assert_eq!(count, system.group_get_n_atoms("W").unwrap());
    }

    #[test]
    fn filter_atoms_without_positions() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        // remove positions
        for atom in system.atoms_iter_mut() {
            atom.reset_position();
        }

        let sphere = Sphere::new([0.0, 0.0, 0.0].into(), 30.0);
        let count = system
            .group_iter("W")
            .unwrap()
            .filter_geometry(sphere)
            .count();

        assert_eq!(count, 0);
    }

    #[test]
    fn bonded_atoms_iter() {
        let mut system = System::from_file("test_files/example.pdb").unwrap();
        system
            .add_bonds_from_pdb("test_files/bonds_for_example.pdb")
            .unwrap();

        let expected_numbers = [28, 30, 32, 36, 38, 42, 48];

        for (i, bonded) in system.bonded_atoms_iter(28).unwrap().enumerate() {
            assert_eq!(bonded.get_atom_number(), expected_numbers[i]);
        }

        assert_eq!(system.bonded_atoms_iter(28).unwrap().count(), 7);
        assert_eq!(system.bonded_atoms_iter(49).unwrap().count(), 0);
    }

    #[test]
    fn bonded_atoms_iter_mut() {
        let mut system = System::from_file("test_files/example.pdb").unwrap();
        system
            .add_bonds_from_pdb("test_files/bonds_for_example.pdb")
            .unwrap();

        let expected_numbers = [28, 30, 32, 36, 38, 42, 48];

        for (i, bonded) in system.bonded_atoms_iter_mut(28).unwrap().enumerate() {
            bonded.set_atom_name("ATM");
            assert_eq!(bonded.get_atom_number(), expected_numbers[i]);
            assert_eq!(bonded.get_atom_name(), "ATM");
        }

        assert_eq!(system.bonded_atoms_iter_mut(28).unwrap().count(), 7);
        assert_eq!(system.bonded_atoms_iter_mut(49).unwrap().count(), 0);
    }

    #[test]
    fn bonded_atoms_iter_fail() {
        let mut system = System::from_file("test_files/example.pdb").unwrap();
        system
            .add_bonds_from_pdb("test_files/bonds_for_example.pdb")
            .unwrap();

        match system.bonded_atoms_iter(50) {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(AtomError::OutOfRange(e)) => assert_eq!(e, 50),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn bonded_atoms_iter_mut_fail() {
        let mut system = System::from_file("test_files/example.pdb").unwrap();
        system
            .add_bonds_from_pdb("test_files/bonds_for_example.pdb")
            .unwrap();

        match system.bonded_atoms_iter_mut(50) {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(AtomError::OutOfRange(e)) => assert_eq!(e, 50),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn molecule_iter_index0() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        let expected_numbers = vec![
            1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 11, 13, 12, 14, 15, 16, 17, 18, 19, 20, 21, 24, 22, 23,
            25, 26, 27, 28, 29, 30, 32, 36, 38, 42, 48, 31, 33, 34, 35, 37, 39, 41, 45, 49, 40, 43,
            46, 44, 47,
        ];

        for (index, atom) in system.molecule_iter(0).unwrap().enumerate() {
            assert_eq!(atom.get_atom_number(), expected_numbers[index]);
        }
    }

    #[test]
    fn molecule_iter_index28() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        let expected_numbers = vec![
            29, 28, 30, 32, 36, 38, 42, 48, 26, 31, 33, 34, 35, 37, 39, 41, 45, 49, 24, 27, 40, 43,
            46, 20, 25, 44, 47, 18, 21, 16, 19, 22, 23, 15, 17, 13, 14, 10, 8, 11, 6, 9, 12, 3, 7,
            1, 4, 2, 5,
        ];

        for (index, atom) in system.molecule_iter(28).unwrap().enumerate() {
            assert_eq!(atom.get_atom_number(), expected_numbers[index]);
        }
    }

    #[test]
    fn molecule_iter_index49() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        assert_eq!(system.molecule_iter(49).unwrap().count(), 1);
        assert_eq!(
            system
                .molecule_iter(49)
                .unwrap()
                .next()
                .unwrap()
                .get_atom_number(),
            50
        );
    }

    #[test]
    fn molecule_iter_invalid() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        match system.molecule_iter(50) {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(AtomError::OutOfRange(number)) => assert_eq!(number, 50),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn molecule_iter_index0_filter() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        let sphere = Sphere::new(
            system.get_atom(0).unwrap().get_position().unwrap().clone(),
            2.0,
        );

        let expected_numbers = vec![
            1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 11, 13, 12, 14, 15, 16, 17, 18, 19, 20, 21, 24, 22, 23,
            25, 26, 27,
        ];

        for (index, atom) in system
            .molecule_iter(0)
            .unwrap()
            .filter_geometry(sphere)
            .enumerate()
        {
            assert_eq!(atom.get_atom_number(), expected_numbers[index]);
        }
    }

    #[test]
    fn molecule_iter_mut_index0() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        let expected_numbers = vec![
            1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 11, 13, 12, 14, 15, 16, 17, 18, 19, 20, 21, 24, 22, 23,
            25, 26, 27, 28, 29, 30, 32, 36, 38, 42, 48, 31, 33, 34, 35, 37, 39, 41, 45, 49, 40, 43,
            46, 44, 47,
        ];

        for (index, atom) in system.molecule_iter_mut(0).unwrap().enumerate() {
            atom.set_residue_name("MOL");

            assert_eq!(atom.get_atom_number(), expected_numbers[index]);
            assert_eq!(atom.get_residue_name(), "MOL");
        }
    }

    #[test]
    fn molecule_iter_mut_index28() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        let expected_numbers = vec![
            29, 28, 30, 32, 36, 38, 42, 48, 26, 31, 33, 34, 35, 37, 39, 41, 45, 49, 24, 27, 40, 43,
            46, 20, 25, 44, 47, 18, 21, 16, 19, 22, 23, 15, 17, 13, 14, 10, 8, 11, 6, 9, 12, 3, 7,
            1, 4, 2, 5,
        ];

        for (index, atom) in system.molecule_iter_mut(28).unwrap().enumerate() {
            atom.set_residue_name("MOL");

            assert_eq!(atom.get_atom_number(), expected_numbers[index]);
            assert_eq!(atom.get_residue_name(), "MOL");
        }
    }

    #[test]
    fn molecule_iter_mut_invalid() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        match system.molecule_iter_mut(50) {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(AtomError::OutOfRange(number)) => assert_eq!(number, 50),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn molecule_iter_mut_index0_filter() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        let sphere = Sphere::new(
            system.get_atom(0).unwrap().get_position().unwrap().clone(),
            2.0,
        );

        let expected_numbers = vec![
            1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 11, 13, 12, 14, 15, 16, 17, 18, 19, 20, 21, 24, 22, 23,
            25, 26, 27,
        ];

        for (index, atom) in system
            .molecule_iter_mut(0)
            .unwrap()
            .filter_geometry(sphere)
            .enumerate()
        {
            atom.set_residue_name("MOL");

            assert_eq!(atom.get_atom_number(), expected_numbers[index]);
            assert_eq!(atom.get_residue_name(), "MOL");
        }
    }

    #[test]
    fn selection_iter() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let mut count = 0;
        for (atom1, atom2) in system
            .selection_iter("@protein")
            .unwrap()
            .zip(system.group_iter("Protein").unwrap())
        {
            compare_atoms(atom1, atom2);
            count += 1;
        }

        assert_eq!(count, 61);
    }

    #[test]
    fn selection_iter_filter_geometry() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let cylinder = Cylinder::new(
            system.group_get_center("Protein").unwrap(),
            2.0,
            3.0,
            Dimension::X,
        );

        let mut count = 0;
        for (atom1, atom2) in system
            .selection_iter("resname W")
            .unwrap()
            .filter_geometry(cylinder.clone())
            .zip(system.group_iter("W").unwrap().filter_geometry(cylinder))
        {
            compare_atoms(atom1, atom2);
            count += 1;
        }

        assert_eq!(count, 29);
    }

    #[test]
    fn selection_iter_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let mut count = 0;
        for atom in system.selection_iter_mut("@protein").unwrap() {
            atom.set_position(Vector3D::new(1.0, 2.0, 3.0));
            count += 1;
        }
        assert_eq!(count, 61);

        for atom in system.group_iter("Protein").unwrap() {
            let pos = atom.get_position().unwrap();
            assert_approx_eq!(f32, pos.x, 1.0);
            assert_approx_eq!(f32, pos.y, 2.0);
            assert_approx_eq!(f32, pos.z, 3.0);
        }
    }

    #[test]
    fn selection_iter_mut_filter_geometry() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let cylinder = Cylinder::new(
            system.group_get_center("Protein").unwrap(),
            2.0,
            3.0,
            Dimension::X,
        );

        let mut count = 0;
        for atom in system
            .selection_iter_mut("resname W")
            .unwrap()
            .filter_geometry(cylinder.clone())
        {
            atom.set_velocity(Vector3D::new(1.0, 2.0, 3.0));
            count += 1;
        }
        assert_eq!(count, 29);

        for atom in system.group_iter("W").unwrap().filter_geometry(cylinder) {
            let vel = atom.get_velocity().unwrap();
            assert_approx_eq!(f32, vel.x, 1.0);
            assert_approx_eq!(f32, vel.y, 2.0);
            assert_approx_eq!(f32, vel.z, 3.0);
        }
    }

    #[test]
    fn molecule_bonds_iter_artificial() {
        let atom1 = Atom::new(1, "LYS", 1, "BB");
        let atom2 = Atom::new(1, "LYS", 2, "SC1");
        let atom3 = Atom::new(1, "LYS", 3, "SC2");
        let atom4 = Atom::new(1, "LYS", 4, "SC3");

        let mut system = System::new("Test", vec![atom1, atom2, atom3, atom4], None);

        // no bonds
        let mut iterator = system.molecule_bonds_iter(1).unwrap();
        assert!(iterator.next().is_none());

        // linear molecule
        system.add_bond(0, 1).unwrap();
        system.add_bond(1, 2).unwrap();
        system.add_bond(2, 3).unwrap();

        let iterator = system.molecule_bonds_iter(2).unwrap();
        let expected = [(1, 2), (2, 3), (0, 1)];
        for ((a, b), (i, j)) in iterator.zip(expected.into_iter()) {
            assert_eq!(a.get_index(), i);
            assert_eq!(b.get_index(), j);
        }

        // cyclic molecule
        system.add_bond(0, 3).unwrap();
        let iterator = system.molecule_bonds_iter(2).unwrap();
        let expected = [(1, 2), (2, 3), (0, 1), (0, 3)];
        for ((a, b), (i, j)) in iterator.zip(expected.into_iter()) {
            assert_eq!(a.get_index(), i);
            assert_eq!(b.get_index(), j);
        }

        // all atoms bonded to all atoms
        system.add_bond(0, 2).unwrap();
        system.add_bond(1, 3).unwrap();

        let iterator = system.molecule_bonds_iter(2).unwrap();
        let expected = [(0, 2), (1, 2), (2, 3), (0, 1), (0, 3), (1, 3)];
        for ((a, b), (i, j)) in iterator.zip(expected.into_iter()) {
            assert_eq!(a.get_index(), i);
            assert_eq!(b.get_index(), j);
        }

        // larger molecule
        // 1 --- 2 --- 3
        // |            \___
        // |                \
        // 4 --- 5 --- 6 --- 7
        // |
        // |
        // 8 --- 9
        let atom1 = Atom::new(1, "MOL", 1, "A");
        let atom2 = Atom::new(1, "MOL", 2, "B");
        let atom3 = Atom::new(1, "MOL", 3, "C");
        let atom4 = Atom::new(1, "MOL", 4, "D");
        let atom5 = Atom::new(1, "MOL", 5, "E");
        let atom6 = Atom::new(1, "MOL", 6, "F");
        let atom7 = Atom::new(1, "MOL", 7, "G");
        let atom8 = Atom::new(1, "MOL", 8, "H");
        let atom9 = Atom::new(1, "MOL", 9, "I");

        let mut system = System::new(
            "Test",
            vec![
                atom1, atom2, atom3, atom4, atom5, atom6, atom7, atom8, atom9,
            ],
            None,
        );

        system.add_bond(0, 1).unwrap();
        system.add_bond(0, 3).unwrap();
        system.add_bond(1, 2).unwrap();
        system.add_bond(2, 6).unwrap();
        system.add_bond(3, 4).unwrap();
        system.add_bond(3, 7).unwrap();
        system.add_bond(7, 8).unwrap();
        system.add_bond(4, 5).unwrap();
        system.add_bond(5, 6).unwrap();

        let iterator = system.molecule_bonds_iter(0).unwrap();
        let expected = [
            (0, 1),
            (0, 3),
            (1, 2),
            (3, 4),
            (3, 7),
            (2, 6),
            (4, 5),
            (7, 8),
            (5, 6),
        ];
        for ((a, b), (i, j)) in iterator.zip(expected.into_iter()) {
            assert_eq!(a.get_index(), i);
            assert_eq!(b.get_index(), j);
        }
    }

    #[test]
    fn molecule_bonds_iter() {
        let system = System::from_file("test_files/example.tpr").unwrap();

        let iterator = system.molecule_bonds_iter(69).unwrap();
        let expected = [
            (64, 69),
            (69, 70),
            (63, 64),
            (70, 71),
            (62, 63),
            (63, 65),
            (71, 72),
            (61, 62),
            (65, 66),
            (66, 67),
            (67, 68),
        ];

        for ((a, b), (i, j)) in iterator.zip(expected.into_iter()) {
            assert_eq!(a.get_index(), i);
            assert_eq!(b.get_index(), j);
        }
    }

    #[test]
    fn molecule_bonds_iter_mut_artificial() {
        let atom1 = Atom::new(1, "LYS", 1, "BB");
        let atom2 = Atom::new(1, "LYS", 2, "SC1");
        let atom3 = Atom::new(1, "LYS", 3, "SC2");
        let atom4 = Atom::new(1, "LYS", 4, "SC3");

        let mut system = System::new("Test", vec![atom1, atom2, atom3, atom4], None);

        // no bonds
        let mut iterator = system.molecule_bonds_iter_mut(1).unwrap();
        assert!(iterator.next().is_none());

        // linear molecule
        system.add_bond(0, 1).unwrap();
        system.add_bond(1, 2).unwrap();
        system.add_bond(2, 3).unwrap();

        let iterator = system.molecule_bonds_iter_mut(2).unwrap();
        let expected = [(1, 2), (2, 3), (0, 1)];
        for ((a, b), (i, j)) in iterator.zip(expected.into_iter()) {
            assert_eq!(a.get_index(), i);
            assert_eq!(b.get_index(), j);
        }

        // cyclic molecule
        system.add_bond(0, 3).unwrap();
        let iterator = system.molecule_bonds_iter_mut(2).unwrap();
        let expected = [(1, 2), (2, 3), (0, 1), (0, 3)];
        for ((a, b), (i, j)) in iterator.zip(expected.into_iter()) {
            assert_eq!(a.get_index(), i);
            assert_eq!(b.get_index(), j);
        }

        // all atoms bonded to all atoms
        system.add_bond(0, 2).unwrap();
        system.add_bond(1, 3).unwrap();

        let iterator = system.molecule_bonds_iter_mut(2).unwrap();
        let expected = [(0, 2), (1, 2), (2, 3), (0, 1), (0, 3), (1, 3)];
        for ((a, b), (i, j)) in iterator.zip(expected.into_iter()) {
            assert_eq!(a.get_index(), i);
            assert_eq!(b.get_index(), j);
        }

        // larger molecule
        // 1 --- 2 --- 3
        // |            \___
        // |                \
        // 4 --- 5 --- 6 --- 7
        // |
        // |
        // 8 --- 9
        let atom1 = Atom::new(1, "MOL", 1, "A");
        let atom2 = Atom::new(1, "MOL", 2, "B");
        let atom3 = Atom::new(1, "MOL", 3, "C");
        let atom4 = Atom::new(1, "MOL", 4, "D");
        let atom5 = Atom::new(1, "MOL", 5, "E");
        let atom6 = Atom::new(1, "MOL", 6, "F");
        let atom7 = Atom::new(1, "MOL", 7, "G");
        let atom8 = Atom::new(1, "MOL", 8, "H");
        let atom9 = Atom::new(1, "MOL", 9, "I");

        let mut system = System::new(
            "Test",
            vec![
                atom1, atom2, atom3, atom4, atom5, atom6, atom7, atom8, atom9,
            ],
            None,
        );

        system.add_bond(0, 1).unwrap();
        system.add_bond(0, 3).unwrap();
        system.add_bond(1, 2).unwrap();
        system.add_bond(2, 6).unwrap();
        system.add_bond(3, 4).unwrap();
        system.add_bond(3, 7).unwrap();
        system.add_bond(7, 8).unwrap();
        system.add_bond(4, 5).unwrap();
        system.add_bond(5, 6).unwrap();

        let iterator = system.molecule_bonds_iter_mut(0).unwrap();
        let expected = [
            (0, 1),
            (0, 3),
            (1, 2),
            (3, 4),
            (3, 7),
            (2, 6),
            (4, 5),
            (7, 8),
            (5, 6),
        ];
        for ((a, b), (i, j)) in iterator.zip(expected.into_iter()) {
            assert_eq!(a.get_index(), i);
            assert_eq!(b.get_index(), j);
        }
    }

    #[test]
    fn molecule_bonds_iter_mut() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();

        let iterator = system.molecule_bonds_iter_mut(69).unwrap();
        let expected = [
            (64, 69),
            (69, 70),
            (63, 64),
            (70, 71),
            (62, 63),
            (63, 65),
            (71, 72),
            (61, 62),
            (65, 66),
            (66, 67),
            (67, 68),
        ];

        for ((a, b), (i, j)) in iterator.zip(expected.into_iter()) {
            assert_eq!(a.get_index(), i);
            assert_eq!(b.get_index(), j);
        }
    }

    #[test]
    fn molecule_bonds_iter_fail() {
        let system = System::from_file("test_files/example.tpr").unwrap();

        match system.molecule_bonds_iter(16844) {
            Ok(_) => panic!("Function should have failed."),
            Err(AtomError::OutOfRange(x)) => assert_eq!(x, 16844),
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn molecule_bonds_iter_mut_fail() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();

        match system.molecule_bonds_iter_mut(16844) {
            Ok(_) => panic!("Function should have failed."),
            Err(AtomError::OutOfRange(x)) => assert_eq!(x, 16844),
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }
}
