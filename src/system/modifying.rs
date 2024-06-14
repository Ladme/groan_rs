// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of System methods for modifying the system.

use std::collections::HashSet;

use crate::errors::{AtomError, GroupError, PositionError};
use crate::structures::iterators::MasterMutAtomIterator;
use crate::structures::{
    atom::Atom,
    simbox::{simbox_check, SimBox},
    vector3d::Vector3D,
};
use crate::system::System;

/// ## Methods for modifying the properties of the system.
impl System {
    /// Translate all atoms of a group by target vector.
    ///
    /// ## Returns
    /// - `Ok` if everything was successful.
    /// - `GroupError::NotFound` in case the group does not exist.
    /// - `GroupError::InvalidSimBox` if the system has no simulation box
    /// or the simulation box is not orthogonal.
    /// - `GroupError::InvalidPosition` if any of the atoms of the group
    /// has an undefined position.
    ///
    /// ## Example
    /// Translating the atoms of the group "Protein".
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("Protein", "resid 1 to 29");
    ///
    /// match system.group_translate("Protein", &[1.0, 2.0, -1.0].into()) {
    ///     Err(e) => eprintln!("{}", e),
    ///     Ok(_) => (),   
    /// }
    /// ```
    #[inline(always)]
    pub fn group_translate(&mut self, name: &str, vector: &Vector3D) -> Result<(), GroupError> {
        match self.group_iter_mut(name)?.translate(vector) {
            Ok(_) => Ok(()),
            Err(AtomError::InvalidPosition(e)) => Err(GroupError::InvalidPosition(e)),
            Err(AtomError::InvalidSimBox(e)) => Err(GroupError::InvalidSimBox(e)),
            _ => panic!("FATAL GROAN ERROR | System::group_translate | Invalid error type returned from `MasterMutAtomIterator::translate`.")
        }
    }

    /// Translate all atoms in the system by target vector.
    ///
    /// ## Returns
    /// - `Ok` if everything was successful.
    /// - `AtomError::InvalidSimBox` if the system has no simulation box
    /// or the simulation box is not orthogonal.
    /// - `AtomError::InvalidPosition` if any of the atoms of the system
    /// has an undefined position.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// system.atoms_translate(&[1.0, 2.0, -1.0].into());
    /// ```
    #[inline(always)]
    pub fn atoms_translate(&mut self, vector: &Vector3D) -> Result<(), AtomError> {
        self.atoms_iter_mut().translate(vector)
    }

    /// Renumber all atoms of the system. This function will give a new atom number
    /// to each atom depending on the index of the atom. The atom numbers start with 1.
    ///
    /// ## Example
    /// Constructing a new system containing a dimer
    /// of a protein from the original system.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// // load system and ndx groups from files
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// // copy protein atoms to a new vector
    /// let mut protein = system.group_extract("Protein").unwrap();
    /// // copy protein atoms again (second protomer)
    /// let mut protein2 = protein.clone();
    ///
    /// // translate atoms of the second protomer
    /// let translate = Vector3D::new(2.0, 0.0, 0.0);
    /// for atom in protein2.iter_mut() {
    ///     atom.translate_nopbc(&translate);
    /// }
    ///
    /// // add atoms of the second protomer to the first protomer
    /// protein.extend(protein2);
    ///
    /// // create new system
    /// let mut new_system = System::new("New system", protein, system.get_box_copy());
    /// // the atom numbers in this system will not be unique...
    ///
    /// // give new (correct) numbers to the atoms
    /// new_system.atoms_renumber();
    ///
    /// // write a new gro file with correct atom numbers
    /// new_system.write_gro("output.gro", true).unwrap();
    /// ```
    pub fn atoms_renumber(&mut self) {
        for (i, atom) in self.atoms_iter_mut().enumerate() {
            atom.set_atom_number(i + 1);
        }
    }

    /// Renumber all residues of the system. This function will give a new residue number
    /// to each atom based on the position of the residue in the list of atoms.
    /// The residue numbers start with 1.
    ///
    /// ## Example
    /// Constructing a new system containing a dimer
    /// of a protein from the original system.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// // load system and ndx groups from files
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// // copy protein atoms to a new vector
    /// let mut protein = system.group_extract("Protein").unwrap();
    /// // copy protein atoms again (second protomer)
    /// let mut protein2 = protein.clone();
    ///
    /// // translate atoms of the second protomer
    /// let translate = Vector3D::new(2.0, 0.0, 0.0);
    /// for atom in protein2.iter_mut() {
    ///     atom.translate_nopbc(&translate);
    /// }
    ///
    /// // add atoms of the second protomer to the first protomer
    /// protein.extend(protein2);
    ///
    /// // create new system
    /// let mut new_system = System::new("New system", protein, system.get_box_copy());
    ///
    /// // give new (correct) numbers to the atoms and residues
    /// new_system.atoms_renumber();
    /// new_system.residues_renumber();
    ///
    /// // write a new gro file with correct atom and residue numbers
    /// new_system.write_gro("output.gro", true).unwrap();
    /// ```
    ///
    /// ## Notes
    /// - In case the residues are 'broken' meaning that atoms of one residue do not follow each other,
    /// this function will not be able to renumber the residues correctly.
    /// In other words, if your gro file looks like this:
    /// ```text
    /// 1GLN      N    1   5.349   9.908   1.871 -0.3054  0.4903 -0.0291
    /// 1GLN    HT1    2   5.293   9.941   1.951  0.0970  1.0823  0.0129
    /// 1GLN    HT2    3   5.293   9.881   1.788 -1.7290  0.8706  0.8040
    /// 2GLU      N    4   5.642   9.890   2.010 -0.1857  0.7126 -0.0478
    /// 2GLU     HN    5   5.677   9.906   1.918  1.8246  1.2580  0.8023
    /// 1GLN    HT3    6   5.419   9.983   1.854 -0.7520  0.8288 -0.3938
    /// 1GLN     CA    7   5.432   9.785   1.924  0.2711 -0.5576 -0.7626
    /// 1GLN     HA    8   5.363   9.708   1.957 -0.1090  0.8599  1.9361
    /// ```
    /// After renumbering, atoms 1-3 will have a residue number 1, atoms 4-5 will have a residue number 2,
    /// and atoms 6-8 will have a residue number 3.
    pub fn residues_renumber(&mut self) {
        let mut current_res = 0;
        let mut renumbered_res = 0;

        for atom in self.atoms_iter_mut() {
            if atom.get_residue_number() != current_res {
                current_res = atom.get_residue_number();
                renumbered_res += 1;
                atom.set_residue_number(renumbered_res);
            } else {
                atom.set_residue_number(renumbered_res);
            }
        }
    }

    /// Wrap atoms of the system into the simulation box.
    ///
    /// ## Returns
    /// `Ok` if everything was successful.
    /// `AtomError::InvalidSimBox` if the system has no simulation box
    /// or the simulation box is not orthogonal.
    /// `AtomError::InvalidPosition` if any of the atoms in the system
    /// has an undefined position.
    #[inline(always)]
    pub fn atoms_wrap(&mut self) -> Result<(), AtomError> {
        self.atoms_iter_mut().wrap()
    }

    /// Wrap atoms of a given group into the simulation box.
    ///
    /// ## Returns
    /// - `Ok` if everything was successful.
    /// - `GroupError::NotFound` if the group does not exist.
    /// - `GroupError::InvalidSimBox` if the system has no simulation box
    /// or the simulation box is not orthogonal.
    /// - `GroupError::InvalidPosition` if any of the atoms of the group
    /// has an undefined position.
    #[inline(always)]
    pub fn group_wrap(&mut self, name: &str) -> Result<(), GroupError> {
        match self.group_iter_mut(name)?.wrap() {
            Ok(_) => Ok(()),
            Err(AtomError::InvalidPosition(e)) => Err(GroupError::InvalidPosition(e)),
            Err(AtomError::InvalidSimBox(e)) => Err(GroupError::InvalidSimBox(e)),
            _ => panic!("FATAL GROAN ERROR | System::group_wrap | Invalid error type returned from `MasterMutAtomIterator::wrap`.")
        }
    }

    /// Add bond connecting two atoms with target indices. Atoms are indexed from 0.
    ///
    /// ## Returns
    /// `Ok` if both atoms exist and the bond was created.
    /// `AtomError::OutOfRange` if any of the atoms does not exist.
    /// `AtomError::InvalidBond` if the bond is invalid.
    ///
    /// In case the bond already exists, nothing happens.
    ///
    /// ## Notes
    /// - This function resets the reference atoms for molecules (`mol_references`) in the system.
    pub fn add_bond(&mut self, index1: usize, index2: usize) -> Result<(), AtomError> {
        if index1 == index2 {
            return Err(AtomError::InvalidBond(index1, index2));
        }

        unsafe {
            let atom1 = self.get_atom_as_mut(index1)? as *mut Atom;
            let atom2 = self.get_atom_as_mut(index2)? as *mut Atom;

            (*atom1).add_bonded(index2);
            (*atom2).add_bonded(index1);
        }

        // reset molecule references
        self.reset_mol_references();

        Ok(())
    }

    /// Analyze topology of the system and select reference atoms
    /// for all polyatomic molecules in the system.
    /// (Reference atom is generally the first atom of the molecule.)
    /// This is used to make `System::make_molecules_whole` more efficient.
    fn create_mol_references(&mut self) {
        let mut visited = HashSet::new();
        let mut new_mol_refs = Vec::new();

        for (a, atom) in self.atoms_iter().enumerate() {
            if visited.contains(&a) {
                continue;
            }

            // ignore monoatomic molecules
            if atom.get_n_bonded() == 0 {
                continue;
            }

            new_mol_refs.push(a);

            visited.insert(a);
            for a2 in crate::system::iterating::get_molecule_indices(self, a).expect(
                "FATAL GROAN ERROR | System::create_mol_references | Atom index does not exist.",
            ) {
                visited.insert(a2);
            }
        }

        self.set_mol_references(new_mol_refs);
    }

    /// Make molecules whole in the simulation box.
    ///
    /// ## Returns
    /// - `Ok` if everything was successful.
    /// - `AtomError::InvalidSimBox` if the system has no simulation box
    /// or the simulation box is not orthogonal.
    /// - `AtomError::InvalidPosition` if any atom that is part of any
    /// polyatomic molecule has no position.
    ///
    /// ## Notes
    /// - Assume you have a system composed of two molecules:
    /// ```text
    ///
    ///   ╔════.═══════════╗
    ///   ║    |     <R>   ║
    ///   ║   <R>     |    ║
    ///   .—o         o——o—.
    ///   ║              | ║
    ///   ║ o            o ║
    ///   ║ |              ║
    ///   .—o——o         o—.
    ///   ╚════.═══════════╝
    ///
    ///
    /// ```
    /// All atoms are indicated by `o` except for the
    /// reference atom of the molecule which is indicated by `<R>`.
    ///
    /// All the atoms are nicely wrapped to the inside of the box,
    /// but the molecules are broken on the periodic boundaries.
    /// This method makes molecules whole, i.e.
    /// it transforms the above system into this:
    /// ```text
    ///     o
    ///     |
    ///  o——o——o═══════════╗
    ///   ║    |      <R>  ║
    ///   ║   <R>      |   ║
    ///   ║           o——o——o
    ///   ║              | ║
    ///   ║              o ║
    ///   ║                ║
    ///   ║                ║
    ///   ╚════════════════╝
    ///
    ///
    /// ```
    /// The reference atom is wrapped into the simulation box, while other
    /// atoms of the molecule are positioned based on the reference atom.
    /// - This function uses `mol_references` from the `System` structure as
    /// the reference atoms for the polyatomic molecules.
    /// In case `mol_references` do not exist, they are generated and stored in the `System` structure.
    /// Note that all functions changing the topology of the `System` MUST reset `mol_references`.
    pub fn make_molecules_whole(&mut self) -> Result<(), AtomError> {
        if self.get_mol_references().is_none() {
            self.create_mol_references();
        }

        let starts = self
            .get_mol_references()
            .expect("FATAL GROAN ERROR | System::make_molecules_whole (1) | `mol_starts` should be `Some` but it is `None`.") 
            as *const Vec<usize>;

        let simbox =
            simbox_check(self.get_box_as_ref()).map_err(AtomError::InvalidSimBox)? as *const SimBox;

        unsafe {
            for index in (*starts).iter() {
                let atom = self.get_atom_as_mut(*index).expect(
                    "FATAL GROAN ERROR | System::make_molecules_whole (2) | Atom index does not exist.",
                ) as *mut Atom;

                // wrap reference atom to the simulation box
                (*atom).wrap(&*simbox)?;

                let ref_atom_position = match (*atom).get_position() {
                    Some(x) => x,
                    None => {
                        return Err(AtomError::InvalidPosition(PositionError::NoPosition(
                            (*atom).get_atom_number(),
                        )))
                    }
                };

                // iterate through other atoms of the molecule
                for atom2 in self
                    .molecule_iter_mut(*index)
                    .expect("FATAL GROAN ERROR | System::make_molecules_whole (4) | Atom index does not exist.") 
                    .skip(1)
                {
                    let atom2_pos = match atom2.get_position() {
                        Some(x) => x,
                        None => return Err(AtomError::InvalidPosition(PositionError::NoPosition(atom2.get_atom_number()))),
                    };

                    // get the shortest vector between the reference atom and the target atom
                    let vector = ref_atom_position.vector_to(atom2_pos, &*simbox);

                    // place the target atom to position based on the shortest vector
                    let new_position = Vector3D::new(
                        ref_atom_position.x + vector.x,
                        ref_atom_position.y + vector.y,
                        ref_atom_position.z + vector.z
                    );

                    atom2.set_position(new_position);
                }
            }
        }

        Ok(())
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use crate::errors::SimBoxError;

    use super::*;
    use float_cmp::assert_approx_eq;
    use std::fs::File;
    use tempfile::NamedTempFile;

    #[test]
    fn atoms_translate() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system
            .atoms_translate(&Vector3D::new(3.5, -1.1, 5.4))
            .unwrap();

        let first = system.atoms_iter().next().unwrap();
        let last = system.atoms_iter().last().unwrap();

        let first_pos = first.get_position();
        let last_pos = last.get_position();

        assert_approx_eq!(f32, first_pos.unwrap().x, 12.997);
        assert_approx_eq!(f32, first_pos.unwrap().y, 0.889);
        assert_approx_eq!(f32, first_pos.unwrap().z, 1.64453);

        assert_approx_eq!(f32, last_pos.unwrap().x, 12.329);
        assert_approx_eq!(f32, last_pos.unwrap().y, 10.086);
        assert_approx_eq!(f32, last_pos.unwrap().z, 7.475);
    }

    #[test]
    fn atoms_translate_fail_simbox() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let shift = Vector3D::new(3.5, -1.1, 5.4);

        system.reset_box();

        match system.atoms_translate(&shift) {
            Err(AtomError::InvalidSimBox(SimBoxError::DoesNotExist)) => (),
            Ok(_) => panic!("Funtion should have failed, but it succeeded."),
            Err(e) => panic!(
                "Function successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn atoms_translate_fail_position() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let shift = Vector3D::new(3.5, -1.1, 5.4);

        system.get_atom_as_mut(15).unwrap().reset_position();

        match system.atoms_translate(&shift) {
            Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => assert_eq!(x, 16),
            Ok(_) => panic!("Funtion should have failed, but it succeeded."),
            Err(e) => panic!(
                "Function successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn group_translate() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system
            .group_translate("all", &Vector3D::new(3.5, -1.1, 5.4))
            .unwrap();

        let first = system.atoms_iter().next().unwrap();
        let last = system.atoms_iter().last().unwrap();

        let first_pos = first.get_position();
        let last_pos = last.get_position();

        assert_approx_eq!(f32, first_pos.unwrap().x, 12.997);
        assert_approx_eq!(f32, first_pos.unwrap().y, 0.889);
        assert_approx_eq!(f32, first_pos.unwrap().z, 1.64453);

        assert_approx_eq!(f32, last_pos.unwrap().x, 12.329);
        assert_approx_eq!(f32, last_pos.unwrap().y, 10.086);
        assert_approx_eq!(f32, last_pos.unwrap().z, 7.475);
    }

    #[test]
    fn group_translate_fail_group() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let shift = Vector3D::new(3.5, -1.1, 5.4);

        match system.group_translate("Nonexistent", &shift) {
            Err(GroupError::NotFound(x)) => assert_eq!(x, "Nonexistent"),
            Ok(_) => panic!("Funtion should have failed, but it succeeded."),
            Err(e) => panic!(
                "Function successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn group_translate_fail_simbox() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let shift = Vector3D::new(3.5, -1.1, 5.4);

        system.reset_box();

        match system.group_translate("all", &shift) {
            Err(GroupError::InvalidSimBox(SimBoxError::DoesNotExist)) => (),
            Ok(_) => panic!("Funtion should have failed, but it succeeded."),
            Err(e) => panic!(
                "Function successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn group_translate_fail_position() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        let shift = Vector3D::new(3.5, -1.1, 5.4);

        system.get_atom_as_mut(15).unwrap().reset_position();

        match system.group_translate("all", &shift) {
            Err(GroupError::InvalidPosition(PositionError::NoPosition(x))) => assert_eq!(x, 16),
            Ok(_) => panic!("Funtion should have failed, but it succeeded."),
            Err(e) => panic!(
                "Function successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn atoms_renumber() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        for atom in system.atoms_iter_mut() {
            atom.set_atom_number(1);
        }

        system.atoms_renumber();

        for (i, atom) in system.atoms_iter().enumerate() {
            assert_eq!(atom.get_atom_number(), i + 1);
        }
    }

    #[test]
    fn residues_renumber_1() {
        let system1 = System::from_file("test_files/example_novelocities.gro").unwrap();
        let mut system2 = System::from_file("test_files/example_novelocities.gro").unwrap();

        system2.get_atom_as_mut(0).unwrap().set_residue_number(3);
        system2.get_atom_as_mut(1).unwrap().set_residue_number(3);

        system2.residues_renumber();

        for (a1, a2) in system1.atoms_iter().zip(system2.atoms_iter()) {
            assert_eq!(a1.get_atom_number(), a2.get_atom_number());
            assert_eq!(a1.get_residue_number(), a2.get_residue_number());
        }
    }

    #[test]
    fn residues_renumber_2() {
        let system = System::from_file("test_files/example_novelocities.gro").unwrap();
        let mut atoms = system.atoms_extract();
        let atoms2 = atoms.clone();

        atoms.extend(atoms2);

        let mut new_system = System::new("New system", atoms, system.get_box_copy());
        new_system.residues_renumber();

        let first_atom = new_system.get_atom_as_ref(0).unwrap();
        let middle_atom = new_system.get_atom_as_ref(50).unwrap();
        let last_atom = new_system.get_atom_as_ref(99).unwrap();

        assert_eq!(first_atom.get_residue_number(), 1);
        assert_eq!(middle_atom.get_residue_number(), 22);
        assert_eq!(last_atom.get_residue_number(), 42);

        assert_eq!(first_atom.get_atom_number(), 1);
        assert_eq!(middle_atom.get_atom_number(), 1);
        assert_eq!(last_atom.get_atom_number(), 50);
    }

    #[test]
    fn atoms_wrap() {
        let system_orig = System::from_file("test_files/example.gro").unwrap();
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let simbox = system.get_box_copy().unwrap();
        let translate1 = Vector3D::new(simbox.x * 3.0, -simbox.y, 0.0);

        for index in [154, 1754, 12345, 4, 37, 0] {
            system
                .get_atom_as_mut(index)
                .unwrap()
                .translate_nopbc(&translate1)
                .unwrap();
        }

        let translate2 = Vector3D::new(0.0, simbox.y, -simbox.z * 2.0);
        for index in [13, 65, 9853, 16843, 7832, 489] {
            system
                .get_atom_as_mut(index)
                .unwrap()
                .translate_nopbc(&translate2)
                .unwrap();
        }

        system.atoms_wrap().unwrap();

        for (a1, a2) in system_orig.atoms_iter().zip(system.atoms_iter()) {
            assert_approx_eq!(
                f32,
                a1.get_position().unwrap().x,
                a2.get_position().unwrap().x,
                epsilon = 0.00001
            );
            assert_approx_eq!(
                f32,
                a1.get_position().unwrap().y,
                a2.get_position().unwrap().y,
                epsilon = 0.00001
            );
            assert_approx_eq!(
                f32,
                a1.get_position().unwrap().z,
                a2.get_position().unwrap().z,
                epsilon = 0.00001
            );
        }
    }

    #[test]
    fn atoms_wrap_fail_simbox() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.reset_box();

        match system.atoms_wrap() {
            Err(AtomError::InvalidSimBox(SimBoxError::DoesNotExist)) => (),
            Ok(_) => panic!("Funtion should have failed, but it succeeded."),
            Err(e) => panic!(
                "Function successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn atoms_wrap_fail_position() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.get_atom_as_mut(15).unwrap().reset_position();

        match system.atoms_wrap() {
            Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => assert_eq!(x, 16),
            Ok(_) => panic!("Funtion should have failed, but it succeeded."),
            Err(e) => panic!(
                "Function successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn group_wrap() {
        let system_orig = System::from_file("test_files/example.gro").unwrap();
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.read_ndx("test_files/index.ndx").unwrap();

        let simbox = system.get_box_copy().unwrap();
        let translate1 = Vector3D::new(simbox.x * 3.0, -simbox.y, 0.0);

        for index in [154, 1754, 12345, 4, 37, 0] {
            system
                .get_atom_as_mut(index)
                .unwrap()
                .translate_nopbc(&translate1)
                .unwrap();
        }

        let translate2 = Vector3D::new(0.0, simbox.y, -simbox.z * 2.0);
        for index in [13, 65, 9853, 16843, 7832, 489] {
            system
                .get_atom_as_mut(index)
                .unwrap()
                .translate_nopbc(&translate2)
                .unwrap();
        }

        system.group_wrap("Protein").unwrap();

        let nonprotein_translated1 = [154, 1754, 12345];
        let nonprotein_translated2 = [65, 9853, 16843, 7832, 489];

        for (index, (a1, a2)) in system_orig
            .atoms_iter()
            .zip(system.atoms_iter())
            .enumerate()
        {
            if nonprotein_translated1.contains(&index) {
                assert_approx_eq!(
                    f32,
                    a1.get_position().unwrap().x + translate1.x,
                    a2.get_position().unwrap().x,
                    epsilon = 0.00001
                );
                assert_approx_eq!(
                    f32,
                    a1.get_position().unwrap().y + translate1.y,
                    a2.get_position().unwrap().y,
                    epsilon = 0.00001
                );
                assert_approx_eq!(
                    f32,
                    a1.get_position().unwrap().z + translate1.z,
                    a2.get_position().unwrap().z,
                    epsilon = 0.00001
                );
            } else if nonprotein_translated2.contains(&index) {
                assert_approx_eq!(
                    f32,
                    a1.get_position().unwrap().x + translate2.x,
                    a2.get_position().unwrap().x,
                    epsilon = 0.00001
                );
                assert_approx_eq!(
                    f32,
                    a1.get_position().unwrap().y + translate2.y,
                    a2.get_position().unwrap().y,
                    epsilon = 0.00001
                );
                assert_approx_eq!(
                    f32,
                    a1.get_position().unwrap().z + translate2.z,
                    a2.get_position().unwrap().z,
                    epsilon = 0.00001
                );
            } else {
                assert_approx_eq!(
                    f32,
                    a1.get_position().unwrap().x,
                    a2.get_position().unwrap().x,
                    epsilon = 0.00001
                );
                assert_approx_eq!(
                    f32,
                    a1.get_position().unwrap().y,
                    a2.get_position().unwrap().y,
                    epsilon = 0.00001
                );
                assert_approx_eq!(
                    f32,
                    a1.get_position().unwrap().z,
                    a2.get_position().unwrap().z,
                    epsilon = 0.00001
                );
            }
        }
    }

    #[test]
    fn group_wrap_fail_group() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.group_wrap("Protein") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(GroupError::NotFound(g)) => assert_eq!(g, "Protein"),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{:?}` was returned",
                e
            ),
        }
    }

    #[test]
    fn group_wrap_fail_simbox() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.reset_box();

        match system.group_wrap("all") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(GroupError::InvalidSimBox(SimBoxError::DoesNotExist)) => (),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{:?}` was returned",
                e
            ),
        }
    }

    #[test]
    fn group_wrap_fail_position() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.get_atom_as_mut(15).unwrap().reset_position();

        match system.group_wrap("all") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(GroupError::InvalidPosition(PositionError::NoPosition(x))) => assert_eq!(x, 16),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{:?}` was returned",
                e
            ),
        }
    }

    #[test]
    fn add_bond() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        for a in 1..system.get_n_atoms() {
            system.add_bond(0, a).unwrap();
        }

        assert!(system.has_bonds());

        for atom in system.atoms_iter().skip(1) {
            assert_eq!(atom.get_n_bonded(), 1);
            assert!(atom.get_bonded().isin(0));
        }

        system.add_bond(1, 3).unwrap();

        let atom1 = system.get_atom_as_ref(1).unwrap();
        assert_eq!(atom1.get_n_bonded(), 2);
        assert!(atom1.get_bonded().isin(0));
        assert!(atom1.get_bonded().isin(3));

        let atom3 = system.get_atom_as_ref(3).unwrap();
        assert_eq!(atom3.get_n_bonded(), 2);
        assert!(atom3.get_bonded().isin(0));
        assert!(atom3.get_bonded().isin(1));
    }

    #[test]
    fn add_bond_fail_outofrange_1() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.add_bond(15, 102743) {
            Err(AtomError::OutOfRange(e)) => assert_eq!(e, 102743),
            Ok(_) => panic!("Funtion should have failed, but it succeeded."),
            Err(e) => panic!(
                "Function successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn add_bond_fail_outofrange_2() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.add_bond(102743, 15) {
            Err(AtomError::OutOfRange(e)) => assert_eq!(e, 102743),
            Ok(_) => panic!("Funtion should have failed, but it succeeded."),
            Err(e) => panic!(
                "Function successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn add_bond_fail_selfbonding() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.add_bond(15, 15) {
            Err(AtomError::InvalidBond(i, j)) => assert_eq!((i, j), (15, 15)),
            Ok(_) => panic!("Funtion should have failed, but it succeeded."),
            Err(e) => panic!(
                "Function successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn prepare_topology() {
        let mut system = System::from_file("test_files/multiple_molecules_conect.pdb").unwrap();
        system
            .add_bonds_from_pdb("test_files/multiple_molecules_conect.pdb")
            .unwrap();

        assert_eq!(system.get_mol_references(), None);

        system.create_mol_references();

        assert_eq!(system.get_mol_references(), Some(&vec![0, 5, 33]));
    }

    #[test]
    fn add_bond_topology() {
        let mut system = System::from_file("test_files/multiple_molecules_conect.pdb").unwrap();
        system
            .add_bonds_from_pdb("test_files/multiple_molecules_conect.pdb")
            .unwrap();

        system.create_mol_references();

        system.add_bond(10, 15).unwrap();

        assert_eq!(system.get_mol_references(), None);
    }

    #[test]
    fn make_molecules_whole_basic() {
        let atom1 = Atom::new(1, "RES", 1, "ATM").with_position([6.0, 6.0, 2.0].into());

        let atom2 = Atom::new(1, "RES", 2, "ATM").with_position([1.0, 4.0, 2.0].into());

        let atom3 = Atom::new(1, "RES", 2, "ATM").with_position([4.0, 1.0, 2.0].into());

        let atoms = vec![atom1, atom2, atom3];

        let mut system = System::new("System", atoms.clone(), Some([5.0, 5.0, 5.0].into()));

        system.add_bond(0, 1).unwrap();
        system.add_bond(0, 2).unwrap();
        system.make_molecules_whole().unwrap();

        let atom1 = system.atoms_iter().next().unwrap();
        let atom2 = system.atoms_iter().nth(1).unwrap();
        let atom3 = system.atoms_iter().nth(2).unwrap();

        assert_eq!(atom1.get_position().unwrap().x, 1.0);
        assert_eq!(atom1.get_position().unwrap().y, 1.0);
        assert_eq!(atom1.get_position().unwrap().z, 2.0);

        assert_eq!(atom2.get_position().unwrap().x, 1.0);
        assert_eq!(atom2.get_position().unwrap().y, -1.0);
        assert_eq!(atom2.get_position().unwrap().z, 2.0);

        assert_eq!(atom3.get_position().unwrap().x, -1.0);
        assert_eq!(atom3.get_position().unwrap().y, 1.0);
        assert_eq!(atom3.get_position().unwrap().z, 2.0);

        let mut system = System::new("System", atoms.clone(), Some([5.0, 5.0, 5.0].into()));

        system.add_bond(1, 2).unwrap();
        system.make_molecules_whole().unwrap();

        let atom1 = system.atoms_iter().next().unwrap();
        let atom2 = system.atoms_iter().nth(1).unwrap();
        let atom3 = system.atoms_iter().nth(2).unwrap();

        assert_eq!(atom1.get_position().unwrap().x, 6.0);
        assert_eq!(atom1.get_position().unwrap().y, 6.0);
        assert_eq!(atom1.get_position().unwrap().z, 2.0);

        assert_eq!(atom2.get_position().unwrap().x, 1.0);
        assert_eq!(atom2.get_position().unwrap().y, 4.0);
        assert_eq!(atom2.get_position().unwrap().z, 2.0);

        assert_eq!(atom3.get_position().unwrap().x, -1.0);
        assert_eq!(atom3.get_position().unwrap().y, 6.0);
        assert_eq!(atom3.get_position().unwrap().z, 2.0);
    }

    #[test]
    fn make_molecules_whole() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();
        system.atoms_translate(&[3.5, 4.5, -3.0].into()).unwrap();
        system.make_molecules_whole().unwrap();

        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path();

        system.write_gro(path_to_output, false).unwrap();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/whole_molecules_expected.gro").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn make_molecules_whole_fail_simbox() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();
        system.reset_box();

        match system.make_molecules_whole() {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(AtomError::InvalidSimBox(SimBoxError::DoesNotExist)) => (),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{:?}` was returned",
                e
            ),
        }
    }

    #[test]
    fn make_molecules_whole_fail_position() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();
        system.get_atom_as_mut(15).unwrap().reset_position();

        match system.make_molecules_whole() {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            // atoms are renumbered, therefore we use 17
            Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => assert_eq!(x, 17),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{:?}` was returned",
                e
            ),
        }
    }
}
