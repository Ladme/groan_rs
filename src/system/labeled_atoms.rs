// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of labeled atoms.

use hashbrown::HashMap;

use crate::{
    errors::AtomLabelError,
    structures::{atom::Atom, group::Group},
    system::System,
};

/// ## Methods for labeling the atoms and selecting labeled atoms.
impl System {
    /// Assign a label to an atom with target index.
    ///
    /// ## Returns
    /// - `Ok` if the atom was successfully labeled.
    /// - `AtomLabelError::AlreadyExistsWarning` if the label already existed and has been now reassigned to a new atom.
    /// - `AtomLabelError::InvalidLabel` if the label is invalid (no label assigned).
    /// - `AtomLabelError::IndexOutOfRange` if the atom index does not exist.
    ///
    /// ## Example
    /// Assign a label `Reference Atom` to atom with an index of 17.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// if let Err(e) = system.label_atom("Reference Atom", 17) {
    ///     eprintln!("{}", e);
    ///     return;
    /// }
    /// ```
    ///
    /// ## Notes
    /// - Atoms are indexed starting from 0.
    /// - In case the same label is already assigned, it is reassigned to the new atom and a warning is raised.
    /// - The following characters are not allowed in labels: '"&|!@()<>=
    pub fn label_atom(&mut self, label: &str, index: usize) -> Result<(), AtomLabelError> {
        if index >= self.get_n_atoms() {
            return Err(AtomLabelError::IndexOutOfRange(index));
        }

        if !crate::aux::name_is_valid(label) {
            return Err(AtomLabelError::InvalidLabel(label.to_owned()));
        }

        self.insert_label(label, index)
    }

    /// Select an atom using a groan selection language query and label the atom.
    /// If the query does not select **a single atom**, an error is returned.
    ///
    /// ## Returns
    /// - `Ok` if a single atom was successfully selected and labeled.
    /// - `AtomLabelError::AlreadyExistsWarning` if the label already existed and has been now reassigned to a new atom.
    /// - `AtomLabelError::InvalidLabel` if the label is invalid (no label assigned).
    /// - `AtomLabelError::InvalidQuery` if the query could not be understood.
    /// - `AtomLabelError::InvalidNumberOfAtoms` if the query selects any other number of atoms than 1.
    ///
    /// ## Example
    /// Assign a label `Reference Atom` to atom with an atom name CA.
    /// The atom will be selected and labeled only if there is only one atom
    /// named CA in the system. Otherwise an error is returned.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// if let Err(e) = system.select_and_label("Reference Atom", "name CA") {
    ///     eprintln!("{}", e);
    ///     return;
    /// }
    /// ```
    ///
    /// ## Notes
    /// - In case the same label is already assigned, it is reassigned to the new atom and a warning is raised.
    /// - The following characters are not allowed in labels: '"&|!@()<>=
    pub fn select_and_label(&mut self, label: &str, query: &str) -> Result<(), AtomLabelError> {
        if !crate::aux::name_is_valid(label) {
            return Err(AtomLabelError::InvalidLabel(label.to_owned()));
        }

        let group = Group::from_query(query, self).map_err(AtomLabelError::InvalidQuery)?;
        let n_atoms = group.get_n_atoms();
        if n_atoms != 1 {
            return Err(AtomLabelError::InvalidNumberOfAtoms(n_atoms));
        }
        let index = group.get_atoms().iter().next().expect(
            "FATAL GROAN ERROR | System::select_and_label | Group should contain exactly one atom.",
        );

        self.insert_label(label, index)
    }

    /// Insert label into the system and return a warning if the label was reassigned.
    fn insert_label(&mut self, label: &str, index: usize) -> Result<(), AtomLabelError> {
        if let Some(x) = self.labeled_atoms.insert(label.to_owned(), index) {
            // do not raise a warning if the label was not reassigned
            if x == index {
                Ok(())
            } else {
                Err(AtomLabelError::AlreadyExistsWarning(
                    label.to_owned(),
                    x,
                    index,
                ))
            }
        } else {
            Ok(())
        }
    }

    /// Check whether a label is present inside the system.
    ///
    /// ## Returns
    /// - `true` if the label exists and is associated to an atom.
    /// - `false` otherwise.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.label_atom("Reference Atom", 17).unwrap();
    ///  
    /// assert!(system.label_exists("Reference Atom"));
    /// ```
    pub fn label_exists(&self, label: &str) -> bool {
        self.labeled_atoms.contains_key(label)
    }

    /// Get reference to an atom with target label.
    ///
    /// ## Returns
    /// Reference to the `Atom` if successful. Otherwise `AtomLabelError::NotFound`.
    ///
    /// ## Example
    /// Assign a label `Reference Atom` to atom with an index of 17.
    /// Then get a reference to the labeled atom.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.label_atom("Reference Atom", 17).unwrap();
    ///
    /// let atom = system.get_labeled_atom_as_ref("Reference Atom").unwrap();
    /// ```
    pub fn get_labeled_atom_as_ref(&self, label: &str) -> Result<&Atom, AtomLabelError> {
        match self.labeled_atoms.get(label) {
            // safety: all labels are checked at construction that they point to a valid atom
            // and atoms cannot be removed from the system
            Some(&i) => Ok(unsafe { self.atoms.get_unchecked(i) }),
            None => Err(AtomLabelError::NotFound(label.to_owned())),
        }
    }

    /// Get mutable reference to an atom with target label.
    ///
    /// ## Returns
    /// Mutable reference to the `Atom` if successful. Otherwise `AtomLabelError::NotFound`.
    ///
    /// ## Example
    /// Assign a label `Reference Atom` to atom with an index of 17.
    /// Then get a mutable reference to the labeled atom.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.label_atom("Reference Atom", 17).unwrap();
    ///
    /// let atom = system.get_labeled_atom_as_mut("Reference Atom").unwrap();
    /// ```
    pub fn get_labeled_atom_as_mut(&mut self, label: &str) -> Result<&mut Atom, AtomLabelError> {
        match self.labeled_atoms.get(label) {
            // safety: all labels are checked at construction that they point to a valid atom
            // and atoms cannot be removed from the system
            Some(&i) => Ok(unsafe { self.atoms.get_unchecked_mut(i) }),
            None => Err(AtomLabelError::NotFound(label.to_owned())),
        }
    }

    /// Get copy of an atom with target label.
    ///
    /// ## Returns
    /// Copy of the `Atom` if successful. Otherwise `AtomLabelError::NotFound`.
    pub fn get_labeled_atom_copy(&self, label: &str) -> Result<Atom, AtomLabelError> {
        match self.labeled_atoms.get(label) {
            // safety: all labels are checked at construction that they point to a valid atom
            // and atoms cannot be removed from the system
            Some(&i) => Ok(unsafe { self.atoms.get_unchecked(i).clone() }),
            None => Err(AtomLabelError::NotFound(label.to_owned())),
        }
    }

    /// Iterate over all labeled atoms. The atoms are immutable.
    ///
    /// ## Returns
    /// An iterator over pairs of label : reference to `Atom`.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// system.label_atom("Reference Atom", 17).unwrap();
    /// // ...assign more labels to atoms...
    ///
    /// for (label, atom) in system.labeled_atoms_iter() {
    ///     // perform some operation with the atom and/or label
    /// }
    /// ```
    ///
    /// ## Notes
    /// - The order of iteration through the pairs is undefined.
    pub fn labeled_atoms_iter(&self) -> impl Iterator<Item = (&str, &Atom)> {
        self.labeled_atoms
            .iter()
            // safety: all labels are checked at construction that they point to a valid atom
            // and atoms cannot be removed from the system
            .map(|(label, &index)| unsafe { (label.as_str(), self.atoms.get_unchecked(index)) })
    }

    /// Iterate over all labeled atoms. The atoms are mutable.
    ///
    /// ## Returns
    /// An iterator over pairs of label : mutable reference to `Atom`.
    ///
    /// ## Notes
    /// - The order of iteration through the pairs is undefined.
    pub fn labeled_atoms_iter_mut(&mut self) -> impl Iterator<Item = (&str, &mut Atom)> {
        let atoms_ptr = self.atoms.as_mut_ptr();
        self.labeled_atoms.iter().map(move |(label, &index)| {
            // safety: all labels are checked at construction that they point to a valid atom
            // and atoms cannot be removed from the system
            unsafe { (label.as_str(), &mut *atoms_ptr.add(index)) }
        })
    }

    /// Get the reference to labeled atoms in a raw form.
    pub(crate) fn get_labeled_atoms(&self) -> &HashMap<String, usize> {
        &self.labeled_atoms
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use float_cmp::assert_approx_eq;

    use crate::{errors::SelectError, test_utilities::utilities::compare_atoms};

    use super::*;

    #[test]
    fn label_atom_pass() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.label_atom("labeled atom", 174).unwrap();
        system.label_atom("first atom", 0).unwrap();
        system.label_atom("last atom", 16843).unwrap();

        assert_eq!(*system.labeled_atoms.get("labeled atom").unwrap(), 174);
        assert!(system.label_exists("labeled atom"));
        assert_eq!(*system.labeled_atoms.get("first atom").unwrap(), 0);
        assert!(system.label_exists("first atom"));
        assert_eq!(*system.labeled_atoms.get("last atom").unwrap(), 16843);
        assert!(system.label_exists("last atom"));
    }

    #[test]
    fn label_atom_overwrite() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.label_atom("labeled atom", 174).unwrap();
        match system.label_atom("labeled atom", 7564) {
            Err(AtomLabelError::AlreadyExistsWarning(label, original, new)) => {
                assert_eq!(&label, "labeled atom");
                assert_eq!(original, 174);
                assert_eq!(new, 7564);
            }
            Ok(_) => panic!("Warning should have been returned but it was not."),
            Err(e) => panic!("Incorrect error type '{}' returned.", e),
        }
        assert_eq!(*system.labeled_atoms.get("labeled atom").unwrap(), 7564);
        assert!(system.label_exists("labeled atom"));

        // does not overwrite
        assert!(system.label_atom("labeled atom", 7564).is_ok());
        assert_eq!(*system.labeled_atoms.get("labeled atom").unwrap(), 7564);
    }

    #[test]
    fn label_atom_invalid_label() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.label_atom("labeled at>m", 174) {
            Err(AtomLabelError::InvalidLabel(label)) => assert_eq!(&label, "labeled at>m"),
            Ok(_) => panic!("Error should have been returned but it was not."),
            Err(e) => panic!("Incorrect error type '{}' returned.", e),
        }

        assert!(!system.label_exists("labeled at>m"));
    }

    #[test]
    fn label_atom_out_of_range() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.label_atom("labeled atom", 16844) {
            Err(AtomLabelError::IndexOutOfRange(index)) => assert_eq!(index, 16844),
            Ok(_) => panic!("Error should have been returned but it was not."),
            Err(e) => panic!("Incorrect error type '{}' returned.", e),
        }

        assert!(!system.label_exists("labeled atom"));
    }

    #[test]
    fn select_and_label_pass() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system
            .select_and_label("labeled atom", "serial 175")
            .unwrap();
        system
            .select_and_label("first atom", "resname GLY")
            .unwrap();
        system.select_and_label("last atom", "resid 11180").unwrap();

        assert_eq!(*system.labeled_atoms.get("labeled atom").unwrap(), 174);
        assert!(system.label_exists("labeled atom"));
        assert_eq!(*system.labeled_atoms.get("first atom").unwrap(), 0);
        assert!(system.label_exists("first atom"));
        assert_eq!(*system.labeled_atoms.get("last atom").unwrap(), 16843);
        assert!(system.label_exists("last atom"));
    }

    #[test]
    fn select_and_label_overwrite() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.label_atom("labeled atom", 174).unwrap();
        match system.select_and_label("labeled atom", "resname W and resid 1901") {
            Err(AtomLabelError::AlreadyExistsWarning(label, original, new)) => {
                assert_eq!(&label, "labeled atom");
                assert_eq!(original, 174);
                assert_eq!(new, 7564);
            }
            Ok(_) => panic!("Warning should have been returned but it was not."),
            Err(e) => panic!("Incorrect error type '{}' returned.", e),
        }
        assert_eq!(*system.labeled_atoms.get("labeled atom").unwrap(), 7564);
        assert!(system.label_exists("labeled atom"));

        // does not overwrite
        assert!(system
            .select_and_label("labeled atom", "serial 7565")
            .is_ok());
        assert_eq!(*system.labeled_atoms.get("labeled atom").unwrap(), 7564);
    }

    #[test]
    fn select_and_label_invalid_label() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.select_and_label("labeled at>m", "serial 175") {
            Err(AtomLabelError::InvalidLabel(label)) => assert_eq!(&label, "labeled at>m"),
            Ok(_) => panic!("Error should have been returned but it was not."),
            Err(e) => panic!("Incorrect error type '{}' returned.", e),
        }

        assert!(!system.label_exists("labeled at>m"));
    }

    #[test]
    fn select_and_label_invalid_query() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.select_and_label("labeled atom", "serial 1S5") {
            Err(AtomLabelError::InvalidQuery(SelectError::InvalidNumber(query))) => {
                assert_eq!(&query, "serial 1S5")
            }
            Ok(_) => panic!("Error should have been returned but it was not."),
            Err(e) => panic!("Incorrect error type '{}' returned.", e),
        }

        assert!(!system.label_exists("labeled atom"));
    }

    #[test]
    fn select_and_label_too_many_atoms() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.select_and_label("labeled atom", "resname CYS") {
            Err(AtomLabelError::InvalidNumberOfAtoms(number)) => assert_eq!(number, 2),
            Ok(_) => panic!("Error should have been returned but it was not."),
            Err(e) => panic!("Incorrect error type '{}' returned.", e),
        }

        assert!(!system.label_exists("labeled atom"));
    }

    #[test]
    fn select_and_label_no_atoms() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.select_and_label("labeled atom", "resname IJR") {
            Err(AtomLabelError::InvalidNumberOfAtoms(number)) => assert_eq!(number, 0),
            Ok(_) => panic!("Error should have been returned but it was not."),
            Err(e) => panic!("Incorrect error type '{}' returned.", e),
        }

        assert!(!system.label_exists("labeled atom"));
    }

    #[test]
    fn get_labeled_atom_as_ref() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.label_atom("labeled atom", 174).unwrap();

        let atom_label = system.get_labeled_atom_as_ref("labeled atom").unwrap();
        let atom_index = system.get_atom_as_ref(174).unwrap();

        compare_atoms(atom_label, atom_index);

        // nonexistent label fail
        match system.get_labeled_atom_as_ref("nonexistent label") {
            Err(AtomLabelError::NotFound(label)) => assert_eq!(label, "nonexistent label"),
            Ok(_) => panic!("Error should have been returned but it was not."),
            Err(e) => panic!("Incorrect error type '{}' returned.", e),
        }
    }

    #[test]
    fn get_labeled_atom_as_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.label_atom("labeled atom", 174).unwrap();

        let mut atom_index = system.get_atom_as_ref(174).unwrap().clone();
        atom_index.set_chain('D');
        let atom_label = system.get_labeled_atom_as_mut("labeled atom").unwrap();
        atom_label.set_chain('D');

        compare_atoms(atom_label, &atom_index);

        // nonexistent label fail
        match system.get_labeled_atom_as_mut("nonexistent label") {
            Err(AtomLabelError::NotFound(label)) => assert_eq!(label, "nonexistent label"),
            Ok(_) => panic!("Error should have been returned but it was not."),
            Err(e) => panic!("Incorrect error type '{}' returned.", e),
        }
    }

    #[test]
    fn get_labeled_atom_as_copy() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.label_atom("labeled atom", 174).unwrap();

        let atom_label = system.get_labeled_atom_copy("labeled atom").unwrap();
        let atom_index = system.get_atom_copy(174).unwrap();

        compare_atoms(&atom_label, &atom_index);

        // nonexistent label fail
        match system.get_labeled_atom_copy("nonexistent label") {
            Err(AtomLabelError::NotFound(label)) => assert_eq!(label, "nonexistent label"),
            Ok(_) => panic!("Error should have been returned but it was not."),
            Err(e) => panic!("Incorrect error type '{}' returned.", e),
        }
    }

    #[test]
    fn labeled_atoms_iter() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.label_atom("labeled atom", 174).unwrap();
        system.label_atom("first atom", 0).unwrap();
        system.label_atom("last atom", 16843).unwrap();

        let mut collection = system.labeled_atoms_iter().collect::<Vec<(&str, &Atom)>>();
        collection.sort_by_key(|x| x.0);

        assert_eq!(collection.len(), 3);

        assert_eq!(collection[0].0, "first atom");
        compare_atoms(collection[0].1, system.get_atom_as_ref(0).unwrap());
        assert_eq!(collection[1].0, "labeled atom");
        compare_atoms(collection[1].1, system.get_atom_as_ref(174).unwrap());
        assert_eq!(collection[2].0, "last atom");
        compare_atoms(collection[2].1, system.get_atom_as_ref(16843).unwrap());

        system.label_atom("one more", 1144).unwrap();

        let collection = system.labeled_atoms_iter().collect::<Vec<(&str, &Atom)>>();
        assert_eq!(collection.len(), 4);
    }

    #[test]
    fn labeled_atoms_iter_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.label_atom("labeled atom", 174).unwrap();
        system.label_atom("first atom", 0).unwrap();
        system.label_atom("last atom", 16843).unwrap();

        let mut counter = 0;
        system.labeled_atoms_iter_mut().for_each(|x| {
            x.1.set_mass(10.5);
            counter += 1;
        });

        assert_eq!(counter, 3);

        system
            .labeled_atoms_iter()
            .for_each(|x| assert_approx_eq!(f32, x.1.get_mass().unwrap(), 10.5));
    }
}
