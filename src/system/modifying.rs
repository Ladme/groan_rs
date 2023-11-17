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
            assert!(atom.get_bonded().unwrap().isin(0));
        }

        unsafe {
            system.add_bond(1, 3).unwrap();
        }

        let atom1 = system.get_atom_as_ref(1).unwrap();
        assert_eq!(atom1.get_n_bonded(), 2);
        assert!(atom1.get_bonded().unwrap().isin(0));
        assert!(atom1.get_bonded().unwrap().isin(3));

        let atom3 = system.get_atom_as_ref(3).unwrap();
        assert_eq!(atom3.get_n_bonded(), 2);
        assert!(atom3.get_bonded().unwrap().isin(0));
        assert!(atom3.get_bonded().unwrap().isin(1));
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
