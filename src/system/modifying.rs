// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of System methods for modifying the system.

use crate::errors::{AtomError, GroupError};
use crate::structures::{atom::Atom, simbox::SimBox, vector3d::Vector3D};
use crate::system::general::System;

/// ## Methods for modifying the properties of the system.
impl System {
    /// Translate all atoms of a group by target vector.
    ///
    /// ## Returns
    /// `Ok` or `GroupError::NotFound` in case the group does not exist.
    ///
    /// ## Example
    /// Translating the atoms of the group "Protein".
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("Protein", "resid 1 to 29");
    ///
    /// match system.group_translate("Protein", &[1.0, 2.0, -1.0].into()) {
    ///     Err(e) => eprintln!("{}", e),
    ///     Ok(_) => (),   
    /// }
    /// ```
    pub fn group_translate(&mut self, name: &str, vector: &Vector3D) -> Result<(), GroupError> {
        unsafe {
            let simbox = self.get_box_as_ref() as *const SimBox;

            for atom in self.group_iter_mut(name)? {
                atom.translate(
                    vector,
                    simbox
                        .as_ref()
                        .expect("FATAL GROAN ERROR | System::group_translate | SimBox is NULL which should not happen.")
                );
            }
        }

        Ok(())
    }

    /// Translate all atoms in the system by target vector.
    ///
    /// ## Example
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// system.atoms_translate(&[1.0, 2.0, -1.0].into());
    /// ```
    pub fn atoms_translate(&mut self, vector: &Vector3D) {
        unsafe {
            let simbox = self.get_box_as_ref() as *const SimBox;

            for atom in self.get_atoms_as_ref_mut().iter_mut() {
                atom.translate(
                    vector,
                    simbox
                        .as_ref()
                        .expect("FATAL GROAN ERROR | System::atoms_translate | SimBox is NULL which should not happen.")
                );
            }
        }
    }

    /// Renumber all atoms of the system. This function will give a new atom number
    /// to each atom depending on the index of the atom. The atom numbers start with 1.
    ///
    /// ## Example
    /// Constructing a new system containing a dimer
    /// of a protein from the original system.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
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
    /// let translate = Vector3D::from([2.0, 0.0, 0.0]);
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
        unsafe {
            for (i, atom) in self.get_atoms_as_ref_mut().iter_mut().enumerate() {
                atom.set_atom_number(i + 1);
            }
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
    /// use groan_rs::prelude::*;
    ///
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
    /// let translate = Vector3D::from([2.0, 0.0, 0.0]);
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
    /// ```ignore
    /// 1GLN      N    1   5.349   9.908   1.871 -0.3054  0.4903 -0.0291
    /// 1GLN    HT1    2   5.293   9.941   1.951  0.0970  1.0823  0.0129
    /// 1GLN    HT2    3   5.293   9.881   1.788 -1.7290  0.8706  0.8040
    /// 2GLU      N    4   5.642   9.890   2.010 -0.1857  0.7126 -0.0478
    /// 2GLU     HN    5   5.677   9.906   1.918  1.8246  1.2580  0.8023
    /// 1GLN    HT3    6   5.419   9.983   1.854 -0.7520  0.8288 -0.3938
    /// 1GLN     CA    7   5.432   9.785   1.924  0.2711 -0.5576 -0.7626
    /// 1GLN     HA    8   5.363   9.708   1.957 -0.1090  0.8599  1.9361
    /// ```
    /// Atoms 1-3 will have a residue number 1, atoms 4-5 will have a residue number 2,
    /// and atoms 6-8 will have a residue number 3 after renumbering.
    pub fn residues_renumber(&mut self) {
        unsafe {
            let mut current_res = 0;
            let mut renumbered_res = 0;

            for atom in self.get_atoms_as_ref_mut().iter_mut() {
                if atom.get_residue_number() != current_res {
                    current_res = atom.get_residue_number();
                    renumbered_res += 1;
                    atom.set_residue_number(renumbered_res);
                } else {
                    atom.set_residue_number(renumbered_res);
                }
            }
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
    /// ## Safety
    /// This method is unsafe because adding bonds can break the topology of the system
    /// such as the number and properties of molecules.
    pub unsafe fn add_bond(&mut self, index1: usize, index2: usize) -> Result<(), AtomError> {
        if index1 == index2 {
            return Err(AtomError::InvalidBond(index1, index2));
        }

        unsafe {
            let atom1 = self.get_atom_as_ref_mut(index1)? as *mut Atom;
            let atom2 = self.get_atom_as_ref_mut(index2)? as *mut Atom;

            (*atom1).add_bonded(index2);
            (*atom2).add_bonded(index1);
        }

        Ok(())
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn atoms_translate() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.atoms_translate(&Vector3D::from([3.5, -1.1, 5.4]));

        let first = system.atoms_iter().next().unwrap();
        let last = system.atoms_iter().last().unwrap();

        let first_pos = first.get_position();
        let last_pos = last.get_position();

        assert_approx_eq!(f32, first_pos.x, 12.997);
        assert_approx_eq!(f32, first_pos.y, 0.889);
        assert_approx_eq!(f32, first_pos.z, 1.64453);

        assert_approx_eq!(f32, last_pos.x, 12.329);
        assert_approx_eq!(f32, last_pos.y, 10.086);
        assert_approx_eq!(f32, last_pos.z, 7.475);
    }

    #[test]
    fn group_translate() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system
            .group_translate("all", &Vector3D::from([3.5, -1.1, 5.4]))
            .unwrap();

        let first = system.atoms_iter().next().unwrap();
        let last = system.atoms_iter().last().unwrap();

        let first_pos = first.get_position();
        let last_pos = last.get_position();

        assert_approx_eq!(f32, first_pos.x, 12.997);
        assert_approx_eq!(f32, first_pos.y, 0.889);
        assert_approx_eq!(f32, first_pos.z, 1.64453);

        assert_approx_eq!(f32, last_pos.x, 12.329);
        assert_approx_eq!(f32, last_pos.y, 10.086);
        assert_approx_eq!(f32, last_pos.z, 7.475);
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

        system2
            .get_atom_as_ref_mut(0)
            .unwrap()
            .set_residue_number(3);
        system2
            .get_atom_as_ref_mut(1)
            .unwrap()
            .set_residue_number(3);

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
    fn add_bond() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        for a in 1..system.get_n_atoms() {
            unsafe {
                system.add_bond(0, a).unwrap();
            }
        }

        assert!(system.has_bonds());

        for atom in system.atoms_iter().skip(1) {
            assert_eq!(atom.get_n_bonded(), 1);
            assert!(atom.get_bonded().isin(0));
        }

        unsafe {
            system.add_bond(1, 3).unwrap();
        }

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

        unsafe {
            match system.add_bond(15, 102743) {
                Err(AtomError::OutOfRange(e)) => assert_eq!(e, 102743),
                Ok(_) => panic!("Funtion should have failed, but it succeeded."),
                Err(e) => panic!(
                    "Function successfully failed but incorrect error type `{:?}` was returned.",
                    e
                ),
            }
        }
    }

    #[test]
    fn add_bond_fail_outofrange_2() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        unsafe {
            match system.add_bond(102743, 15) {
                Err(AtomError::OutOfRange(e)) => assert_eq!(e, 102743),
                Ok(_) => panic!("Funtion should have failed, but it succeeded."),
                Err(e) => panic!(
                    "Function successfully failed but incorrect error type `{:?}` was returned.",
                    e
                ),
            }
        }
    }

    #[test]
    fn add_bond_fail_selfbonding() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        unsafe {
            match system.add_bond(15, 15) {
                Err(AtomError::InvalidBond(i, j)) => assert_eq!((i, j), (15, 15)),
                Ok(_) => panic!("Funtion should have failed, but it succeeded."),
                Err(e) => panic!(
                    "Function successfully failed but incorrect error type `{:?}` was returned.",
                    e
                ),
            }
        }
    }
}
