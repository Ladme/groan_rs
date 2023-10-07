// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of System methods for working with groups.

use std::collections::HashSet;

use indexmap::IndexMap;

use crate::errors::GroupError;
use crate::structures::group::Group;
use crate::structures::shape::Shape;
use crate::system::general::System;

/// ## Methods for working with groups of atoms.
impl System {
    /// Make a group with a given name from the given Groan selection language query.
    ///
    /// ## Returns
    /// - `Ok` if the group was successfully created.
    /// - `GroupError::AlreadyExistsWarning` if the new group has overwritten a previous group.
    /// - `GroupError::InvalidName` if the name of the group is invalid (no group created).
    /// - `GroupError::InvalidQuery` if the query could not be parsed.
    ///
    /// ## Example
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// if let Err(e) = system.group_create("Phosphori", "resname POPE POPG and name P") {
    ///     eprintln!("{}", e);
    ///     return;
    /// }
    /// ```
    ///
    /// ## Notes
    /// - In case a group with the given name already exists, it is replaced with the new group.
    /// - The following characters are not allowed in group names: '"&|!@()<>=
    /// - The group will be created even if the query selects no atoms.
    pub fn group_create(&mut self, name: &str, query: &str) -> Result<(), GroupError> {
        if !Group::name_is_valid(name) {
            return Err(GroupError::InvalidName(name.to_string()));
        }

        let group = match Group::from_query(query, self) {
            Ok(x) => x,
            Err(e) => return Err(GroupError::InvalidQuery(e)),
        };

        unsafe {
            match self.get_groups_as_ref_mut().insert(name.to_string(), group) {
                None => Ok(()),
                Some(_) => Err(GroupError::AlreadyExistsWarning(name.to_string())),
            }
        }
    }

    /// Make a group with a given name from the given Groan selection language query and geometry specification.
    /// The group is NOT dynamically updated.
    ///
    /// ## Returns
    /// - `Ok` if the group was successfully created.
    /// - `GroupError::AlreadyExistsWarning` if the new group has overwritten a previous group.
    /// - `GroupError::InvalidName` if the name of the group is invalid (no group created).
    /// - `GroupError::InvalidQuery` if the query could not be parsed.
    ///
    /// ## Example
    /// Select phosphori atoms which are inside a z-axis oriented cylinder
    /// with a radius of 2 nm and height of 4 nm located at coordinates x = 5 nm, y = 5 nm, z = 3 nm.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// let cylinder = Cylinder::new([5.0, 5.0, 3.0].into(), 2.0, 4.0, Dimension::Z);
    ///
    /// if let Err(e) = system.group_create_from_geometry("Phosphori", "name P", cylinder) {
    ///     eprintln!("{}", e);
    ///     return;
    /// }
    /// ```
    ///
    /// ## Warning
    /// - If you construct the group and then iterate through a trajectory, the group will still contain
    /// the same atoms as initially. In other words, the group is NOT dynamically updated.
    /// - If you want to choose atoms dynamically, it is better to use `AtomIterator` and `filter_geometry` function
    /// while iterating through the trajectory.
    ///
    /// ## Notes
    /// - In case a group with the given name already exists, it is replaced with the new group.
    /// - The following characters are not allowed in group names: '"&|!@()<>=
    /// - The group will be created even if no atoms are selected.
    pub fn group_create_from_geometry(
        &mut self,
        name: &str,
        query: &str,
        geometry: impl Shape,
    ) -> Result<(), GroupError> {
        if !Group::name_is_valid(name) {
            return Err(GroupError::InvalidName(name.to_string()));
        }

        let group = match Group::from_query_and_geometry(query, geometry, self) {
            Ok(x) => x,
            Err(e) => return Err(GroupError::InvalidQuery(e)),
        };

        unsafe {
            match self.get_groups_as_ref_mut().insert(name.to_string(), group) {
                None => Ok(()),
                Some(_) => Err(GroupError::AlreadyExistsWarning(name.to_string())),
            }
        }
    }

    /// Create a group with a given name from the provided atom indices.
    ///
    /// ## Returns
    /// - `Ok` if the group was successfully created without any warnings.
    /// - `GroupError::AlreadyExistsWarning` if the new group has overwritten a previous group.
    /// - `GroupError::InvalidName` if the name of the group is invalid (no group created).
    ///
    /// ## Example
    /// Creating a group "My Group" containing the atoms indexed as 0, 1, 2, 3, 10, 11, and 12.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// match system.group_create_from_indices("My Group", vec![0, 1, 2, 3, 12, 10, 11]) {
    ///     Err(e) => eprintln!("{}", e),
    ///     Ok(_) => (),
    /// }
    /// ```
    ///
    /// ## Notes
    /// - In case a group with the given name already exists, it is replaced with the new group.
    /// - The following characters are not allowed in group names: '"&|!@()<>=
    /// - The created Group will be valid even if the input `atom_indices` vector contains multiple identical indices.
    pub fn group_create_from_indices(
        &mut self,
        name: &str,
        atom_indices: Vec<usize>,
    ) -> Result<(), GroupError> {
        if !Group::name_is_valid(name) {
            return Err(GroupError::InvalidName(name.to_string()));
        }

        let group = Group::from_indices(atom_indices, self.get_n_atoms());
        unsafe {
            match self.get_groups_as_ref_mut().insert(name.to_string(), group) {
                None => Ok(()),
                Some(_) => Err(GroupError::AlreadyExistsWarning(name.to_string())),
            }
        }
    }

    /// Create a group with a given name from the provided atom ranges.
    ///
    /// ## Returns
    /// - `Ok` if the group was successfully created without any warnings.
    /// - `GroupError::AlreadyExistsWarning` if the new group has overwritten a previous group.
    /// - `GroupError::InvalidName` if the name of the group is invalid (no group created).
    ///
    /// ## Example
    /// Creating a group "My Group" containing the atoms of the system with indices 0-25 and 50-75 (including).
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// match system.group_create_from_ranges("My Group", vec![(0, 25), (50, 75)]) {
    ///     Err(e) => eprintln!("{}", e),
    ///     Ok(_) => (),
    /// }
    /// ```
    ///
    /// ## Notes
    /// - In case a group with the given name already exists, it is replaced with the new group.
    /// - The following characters are not allowed in group names: '"&|!@()<>=
    /// - The created Group will be valid even if the input `atom_ranges` vector contains overlapping atom ranges.
    pub fn group_create_from_ranges(
        &mut self,
        name: &str,
        atom_ranges: Vec<(usize, usize)>,
    ) -> Result<(), GroupError> {
        if !Group::name_is_valid(name) {
            return Err(GroupError::InvalidName(name.to_string()));
        }

        let group = Group::from_ranges(atom_ranges, self.get_n_atoms());

        unsafe {
            match self.get_groups_as_ref_mut().insert(name.to_string(), group) {
                None => Ok(()),
                Some(_) => Err(GroupError::AlreadyExistsWarning(name.to_string())),
            }
        }
    }

    /// Split atoms by their residue numbers creating a group for each residue number.
    ///
    /// ## Returns
    /// A tuple of the following form:
    /// - First item: `GroupError::MultipleAlreadyExistWarning` if any new group has overwritten a previous group. Else `Ok`.
    /// - Second item: Vector containing all the names of the generated groups.
    ///
    /// ## Example
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// let (result, residues) = system.atoms_split_by_resid();
    /// if let Err(e) = result {
    ///     eprintln!("{}", e);
    /// }
    ///
    /// // iterate through residues of the system
    /// for resid in residues.iter() {
    ///     for atom in system.group_iter(resid).unwrap() {
    ///         println!("{:?}", atom);
    ///     }
    /// }
    /// ```
    ///
    /// ## Notes
    /// - The names of the groups have the following format: 'resid X' where 'X' is the residue number.
    /// - The atoms of the residues do not have positioned to be next to each other, i.e. groups can be broken.
    /// - The atoms will not be sorted, i.e. their position in the groups will be the same as in the System structure.
    /// - The order of the group names in the output vector will be the same as the order of residues in the system.
    /// - The properties of the atoms are not changed.
    pub fn atoms_split_by_resid(&mut self) -> (Result<(), GroupError>, Vec<String>) {
        match self.group_split_by_resid("all") {
            (Err(GroupError::NotFound(_)), _) => {
                panic!("Groan error. Default group `all` does not exist for the system.")
            }
            (result, vec) => (result, vec),
        }
    }

    /// Split atoms of a given group by their residue numbers creating a group for each residue number.
    ///
    /// ## Returns
    /// A tuple of the following form:
    /// - First item:
    ///     - `GroupError::NotFound` if the given group does not exist.
    ///     - `GroupError::MultipleAlreadyExistWarning` if any new group has overwritten a previous group.     
    ///     - `Ok` otherwise.
    /// - Second item: Vector containing all the names of the generated groups.
    ///
    /// ## Example
    /// Creating a group for each residue of the group 'Protein', i.e. creating a group for each amino acid.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// let (result, residues) = system.group_split_by_resid("Protein");
    /// if let Err(e) = result {
    ///     eprintln!("{}", e);
    /// }
    ///
    /// // iterate through residues of the group `Protein`
    /// for resid in residues.iter() {
    ///     for atom in system.group_iter(resid).unwrap() {
    ///         println!("{:?}", atom);
    ///     }
    /// }
    /// ```
    ///
    /// ## Notes
    /// - The names of the groups have the following format: 'resid X' where 'X' is the residue number.
    /// - The atoms of the residues do not have positioned to be next to each other, i.e. groups can be broken.
    /// - The atoms will not be sorted, i.e. their position in the groups will be the same as in the System structure.
    /// - The order of the group names in the output vector will be the same as the order of residues in the system.
    /// - The properties of the atoms are not changed.
    pub fn group_split_by_resid(&mut self, name: &str) -> (Result<(), GroupError>, Vec<String>) {
        let mut groups: IndexMap<String, Vec<usize>> = IndexMap::new();

        let iterator = match self.group_iter(name) {
            Ok(i) => i,
            Err(e) => return (Err(e), Vec::new()),
        };

        // collect atoms into groups
        for (atomid, atom) in iterator.enumerate() {
            let res = atom.get_residue_number();
            let group_name = format!("resid {}", res);

            match groups.get_mut(&group_name) {
                Some(group) => group.push(atomid),
                None => {
                    groups.insert(group_name, vec![atomid]);
                }
            }
        }

        let mut duplicate_names: HashSet<String> = HashSet::new();

        // create the groups
        for (group_name, atoms) in groups.iter() {
            if self
                .group_create_from_indices(group_name, atoms.clone())
                .is_err()
            // we know that the only possible error to be returned is `GroupError::AlreadyExistsWarning`
            {
                duplicate_names.insert(group_name.to_string());
            }
        }

        if duplicate_names.is_empty() {
            (Ok(()), groups.keys().cloned().collect())
        } else {
            (
                Err(GroupError::MultipleAlreadyExistWarning(Box::new(
                    duplicate_names,
                ))),
                groups.keys().cloned().collect(),
            )
        }
    }

    /// Split atoms by their residue names creating a group for each unique residue name.
    ///
    /// ## Returns
    /// A tuple of the following form:
    /// - First item: `GroupError::MultipleAlreadyExistWarning` if any new group has overwritten a previous group. Else `Ok`.
    /// - Second item: Vector containing all the names of the generated groups.
    ///
    /// ## Example
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// let (result, residues) = system.atoms_split_by_resname();
    /// if let Err(e) = result {
    ///     eprintln!("{}", e);
    /// }
    ///
    /// // iterate through residue types of the system
    /// for resid in residues.iter() {
    ///     for atom in system.group_iter(resid).unwrap() {
    ///         println!("{:?}", atom);
    ///     }
    /// }
    /// ```
    ///
    /// ## Notes
    /// - The names of the groups have the following format: 'resname ABC' where 'ABC' is the residue name.
    /// - The atoms will not be sorted, i.e. their position in the groups will be the same as in the System structure.
    /// - The order of the group names in the output vector will be the same as the order of residues in the system.
    /// - The properties of the atoms are not changed.
    pub fn atoms_split_by_resname(&mut self) -> (Result<(), GroupError>, Vec<String>) {
        match self.group_split_by_resname("all") {
            (Err(GroupError::NotFound(_)), _) => {
                panic!("Groan error. Default group `all` does not exist for the system.")
            }
            (result, vec) => (result, vec),
        }
    }

    /// Split atoms of a given group by their residue names creating a group for each unique residue name.
    ///
    /// ## Returns
    /// A tuple of the following form:
    /// - First item:
    ///     - `GroupError::NotFound` if the given group does not exist.
    ///     - `GroupError::MultipleAlreadyExistWarning` if any new group has overwritten a previous group.     
    ///     - `Ok` otherwise.
    /// - Second item: Vector containing all the names of the generated groups.
    ///
    /// ## Example
    /// Creating a group for each residue type of the group 'Protein', i.e. creating a group for each amino acid type.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// let (result, residues) = system.group_split_by_resname("Protein");
    /// if let Err(e) = result {
    ///     eprintln!("{}", e);
    /// }
    ///
    /// // iterate through residue types of the system
    /// for resid in residues.iter() {
    ///     for atom in system.group_iter(resid).unwrap() {
    ///         println!("{:?}", atom);
    ///     }
    /// }
    /// ```
    ///
    /// ## Notes
    /// - The names of the groups have the following format: 'resname ABC' where 'ABC' is the residue name.
    /// - The atoms will not be sorted, i.e. their position in the groups will be the same as in the System structure.
    /// - The order of the group names in the output vector will be the same as the order of residues in the system.
    /// - The properties of the atoms are not changed.
    pub fn group_split_by_resname(&mut self, name: &str) -> (Result<(), GroupError>, Vec<String>) {
        let mut groups: IndexMap<String, Vec<usize>> = IndexMap::new();

        let iterator = match self.group_iter(name) {
            Ok(i) => i,
            Err(e) => return (Err(e), Vec::new()),
        };

        for (atomid, atom) in iterator.enumerate() {
            let res = atom.get_residue_name();
            let group_name = format!("resname {}", res);

            match groups.get_mut(&group_name) {
                Some(group) => group.push(atomid),
                None => {
                    groups.insert(group_name, vec![atomid]);
                }
            }
        }

        let mut duplicate_names: HashSet<String> = HashSet::new();

        // create groups
        for (group_name, atoms) in groups.iter() {
            if self
                .group_create_from_indices(group_name, atoms.clone())
                .is_err()
            // we know that the only possible error to be returned is `GroupError::AlreadyExistsWarning`
            {
                duplicate_names.insert(group_name.to_string());
            }
        }

        if duplicate_names.is_empty() {
            (Ok(()), groups.keys().cloned().collect())
        } else {
            (
                Err(GroupError::MultipleAlreadyExistWarning(Box::new(
                    duplicate_names,
                ))),
                groups.keys().cloned().collect(),
            )
        }
    }

    /// Make target group ndx-writable, i.e. the group will be written into an ndx file when using `System::write_ndx`.
    ///
    /// ## Returns
    /// `Ok` if successful, `GroupError::NotFound` in case the group does not exist.
    pub fn group_make_writable(&mut self, name: &str) -> Result<(), GroupError> {
        unsafe {
            match self.get_groups_as_ref_mut().get_mut(name) {
                None => return Err(GroupError::NotFound(name.to_owned())),
                Some(group) => group.print_ndx = true,
            }
        }

        Ok(())
    }

    /// Make target group ndx-nonwritable, i.e. the group will NOT be written into an ndx file when using `System::write_ndx`.
    ///
    /// ## Returns
    /// `Ok` if successful, `GroupError::NotFound` in case the group does not exist.
    pub fn group_make_nonwritable(&mut self, name: &str) -> Result<(), GroupError> {
        unsafe {
            match self.get_groups_as_ref_mut().get_mut(name) {
                None => return Err(GroupError::NotFound(name.to_owned())),
                Some(group) => group.print_ndx = false,
            }
        }

        Ok(())
    }

    /// Remove target group from the system.
    ///
    /// ## Returns
    /// `Ok` if successful, `GroupError::NotFound` in case the group does not exist.
    ///
    /// ## Safety
    /// Do not use this function to remove any of the default groups ('all' or 'All').
    /// Doing so would make many functions associated with the `System` structure invalid.
    /// Removing other groups is generally safe to do.
    pub unsafe fn group_remove(&mut self, name: &str) -> Result<(), GroupError> {
        if self.get_groups_as_ref_mut().remove(name).is_none() {
            Err(GroupError::NotFound(name.to_owned()))
        } else {
            Ok(())
        }
    }

    /// Rename target group in the system.
    ///
    /// ## Returns
    /// `Ok` if successful, `GroupError::NotFound` in case the old name of the group does not exist,
    /// `GroupError::AlreadyExistsWarning` in case the new name already exists (it is overwritten).
    ///
    /// ## Safety
    /// Do not use this function to rename any of the default groups ('all' or 'All').
    /// Doing so would make many functions associated with the `System` structure invalid.
    /// Renaming other groups is generally safe to do.
    ///
    /// ## Notes
    /// - Note that if the `new` name already matches name of another group in the `System`,
    /// the previous group with this name is overwritten and `GroupError::AlreadyExistsWarning` is returned.
    pub unsafe fn group_rename(&mut self, old: &str, new: &str) -> Result<(), GroupError> {
        // get the old group
        let group = match self.get_groups_as_ref_mut().remove(old) {
            Some(x) => x,
            None => return Err(GroupError::NotFound(old.to_owned())),
        };

        // reinsert the group with the new name
        match self.get_groups_as_ref_mut().insert(new.to_owned(), group) {
            Some(_) => Err(GroupError::AlreadyExistsWarning(new.to_owned())),
            None => Ok(()),
        }
    }

    /// Check whether a group with a given name exists in the system.
    pub fn group_exists(&self, name: &str) -> bool {
        self.get_groups_as_ref().contains_key(name)
    }

    /// Check whether the target group contains the atom of target index.
    ///
    /// ## Returns
    /// `true` or `false` if the atom is part of the group or not, respectively.
    /// `GroupError` if the group does not exist.
    ///
    /// ## Notes
    /// - Note that `index` corresponds to the atom index in the `System` structure.
    /// The atoms in `System` are numbered starting from 0.
    /// - The time complexity of this operation is somewhere between `O(1)` and `O(n)`
    /// where `n` is the number of atoms in the group.
    pub fn group_isin(&self, name: &str, index: usize) -> Result<bool, GroupError> {
        let group = self
            .get_groups_as_ref()
            .get(name)
            .ok_or(GroupError::NotFound(name.to_string()))?;

        for (start, end) in group.atom_ranges.iter() {
            if index >= *start && index <= *end {
                return Ok(true);
            }
        }

        Ok(false)
    }

    /// Get the number of atoms in target group.
    ///
    /// ## Returns
    /// The number of atoms or `GroupError::NotFound` if the group does not exist.
    ///
    /// ## Notes
    /// - This is NOT an `O(1)` operation. In fact, the complexity of this operation
    /// depends on the complexity of the group. Length of a group composed
    /// of a continuous block of atoms will be computed more quickly
    /// than the length of a group composed of a large number of separate atoms.
    ///
    /// - The time complexity of this operation is thus somewhere between `O(1)` and `O(n)`.
    ///
    /// ## Example:
    /// Get the number of atoms in the group "Protein".
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let system = System::from_file("system.gro").unwrap();
    /// let n_atoms = system.group_get_n_atoms("Protein");
    /// ```
    pub fn group_get_n_atoms(&self, name: &str) -> Result<usize, GroupError> {
        let group = self
            .get_groups_as_ref()
            .get(name)
            .ok_or(GroupError::NotFound(name.to_string()))?;

        Ok(group.get_n_atoms())
    }

    /// Create a new group that is the union of the provided groups.
    ///
    /// ## Returns
    /// `GroupError::NotFound` if any of the input groups does not exist.
    /// `GroupError::AlreadyExistsWarning` if a group with the name `union` already exists and is overwritten.
    /// `Ok` otherwise.
    ///
    /// ## Notes
    /// - If a group with the name `union` already exists, it is overwritten.
    /// - `group1` and `group2` are unchanged.
    pub fn group_union(
        &mut self,
        group1: &str,
        group2: &str,
        union: &str,
    ) -> Result<(), GroupError> {
        let group1 = self
            .get_groups_as_ref()
            .get(group1)
            .ok_or(GroupError::NotFound(group1.to_string()))?;
        let group2 = self
            .get_groups_as_ref()
            .get(group2)
            .ok_or(GroupError::NotFound(group2.to_string()))?;

        let mut atom_ranges =
            Vec::with_capacity(group1.atom_ranges.len() + group2.atom_ranges.len());
        atom_ranges.extend(group1.atom_ranges.iter());
        atom_ranges.extend(group2.atom_ranges.iter());

        let group = Group::from_ranges(atom_ranges, self.get_n_atoms());

        unsafe {
            match self
                .get_groups_as_ref_mut()
                .insert(union.to_string(), group)
            {
                None => Ok(()),
                Some(_) => Err(GroupError::AlreadyExistsWarning(union.to_string())),
            }
        }
    }

    /// Extend the group `group` by the group `extend`.
    ///
    /// ## Returns
    /// `Ok` if the group was extended or `GroupError::NotFound` if any of the input groups does not exist.
    ///
    /// ## Notes
    /// - `extend` is unchanged.
    pub fn group_extend(&mut self, group: &str, extend: &str) -> Result<(), GroupError> {
        match self.group_union(group, extend, group) {
            Ok(_) | Err(GroupError::AlreadyExistsWarning(_)) => Ok(()),
            Err(e) => Err(e),
        }
    }

    /// Get all group names associated with the system including ndx-nonwritable (default) groups.
    ///
    /// ## Example
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("Custom Group", "resid 1 to 3").unwrap();
    ///
    /// let names = system.group_names();
    ///
    /// // the `names` vector contains the user-created `Custom Group` as well as the groups
    /// // `all` and `All` that are implicitly constructed when the `System` structure is created
    /// assert_eq!(names.len(), 3);
    /// assert!(names.contains(&"Custom Group".to_string()));
    /// assert!(names.contains(&"all".to_string()));
    /// assert!(names.contains(&"All".to_string()));
    /// ```
    pub fn group_names(&self) -> Vec<String> {
        self.get_groups_as_ref()
            .keys()
            .map(|key| key.to_owned())
            .collect()
    }

    /// Get all group names assocaited with the system excluding ndx-nonwritable (default) groups.
    ///
    /// ## Example
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("Custom Group", "resid 1 to 3").unwrap();
    ///
    /// let names = system.group_names_writable();
    ///
    /// // the `names` vector only contains the user-created `Custom Group`
    /// assert_eq!(names.len(), 1);
    /// assert!(names.contains(&"Custom Group".to_string()));
    /// assert!(!names.contains(&"all".to_string()));
    /// assert!(!names.contains(&"All".to_string()));
    /// ```
    pub fn group_names_writable(&self) -> Vec<String> {
        self.get_groups_as_ref()
            .iter()
            .filter_map(|(key, group)| {
                if group.print_ndx {
                    Some(key.to_owned())
                } else {
                    None
                }
            })
            .collect()
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structures::dimension::Dimension;
    use crate::structures::shape::*;

    #[test]
    fn group_create_basic() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.group_create("Membrane", "resname POPC").unwrap();

        assert!(system.group_exists("Membrane"));
        assert_eq!(system.group_get_n_atoms("Membrane").unwrap(), 6144);

        for i in 61..=6204 {
            assert!(system.group_isin("Membrane", i).unwrap());
        }

        system.group_create("Protein", "serial 1 to 61").unwrap();

        assert!(system.group_exists("Protein"));
        assert_eq!(system.group_get_n_atoms("Protein").unwrap(), 61);

        for i in 0..=60 {
            assert!(system.group_isin("Protein", i).unwrap());
        }
    }

    #[test]
    fn group_create_basic_fails() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        // invalid group name
        match system.group_create("Memb!rane", "resname POPC") {
            Ok(_) => panic!("Group should have an invalid name but it was accepted."),
            Err(e) => {
                assert!(e.to_string().contains("Memb!rane"));
            }
        }

        assert!(!system.group_exists("Memb!rane"));

        // invalid query 1
        match system.group_create("Membrane", "resname POPC &&") {
            Ok(_) => panic!("Parsing should have failed but it succeeded."),
            Err(e) => {
                assert!(e.to_string().contains("missing argument"));
            }
        }

        assert!(!system.group_exists("Membrane"));

        // invalid query 2
        match system.group_create("Membrane", "(resname POPC && resname POPE))") {
            Ok(_) => panic!("Parsing should have failed but it succeeded."),
            Err(e) => {
                assert!(e.to_string().contains("unmatching parentheses"));
            }
        }
    }

    #[test]
    fn group_create_basic_overwrite() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.group_create("Membrane", "serial 1").unwrap();

        // overwritting previous group
        match system.group_create("Membrane", "resname POPC") {
            Ok(_) => panic!("Warning was not raised."),
            Err(e) => {
                assert!(e.to_string().contains("already existed"));
            }
        }

        assert!(system.group_exists("Membrane"));
        assert_eq!(system.group_get_n_atoms("Membrane").unwrap(), 6144);
        for i in 61..=6204 {
            assert!(system.group_isin("Membrane", i).unwrap());
        }
    }

    #[test]
    fn group_create_chain() {
        let mut system = System::from_file("test_files/example.pdb").unwrap();

        system.group_create("Chains A+B", "chain A B").unwrap();

        assert!(system.group_exists("Chains A+B"));
        assert_eq!(system.group_get_n_atoms("Chains A+B").unwrap(), 31);

        for i in 0..=30 {
            assert!(system.group_isin("Chains A+B", i).unwrap());
        }

        system.group_create("Chain C", "chain C").unwrap();

        assert!(system.group_exists("Chain C"));
        assert_eq!(system.group_get_n_atoms("Chain C").unwrap(), 19);

        for i in 31..=49 {
            assert!(system.group_isin("Chain C", i).unwrap());
        }
    }

    #[test]
    fn group_create_chain_from_gro() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.group_create("Chains A+B", "chain A B").unwrap();

        assert!(system.group_exists("Chains A+B"));
        assert_eq!(system.group_get_n_atoms("Chains A+B").unwrap(), 0);
    }

    #[test]
    fn group_create_macro_protein() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.group_create("Protein", "@protein").unwrap();

        assert!(system.group_exists("Protein"));
        assert_eq!(system.group_get_n_atoms("Protein").unwrap(), 61);

        for i in 0..=60 {
            assert!(system.group_isin("Protein", i).unwrap());
        }
    }

    #[test]
    fn group_create_macro_membrane() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.group_create("Membrane", "@membrane").unwrap();

        assert!(system.group_exists("Membrane"));
        assert_eq!(system.group_get_n_atoms("Membrane").unwrap(), 6144);

        for i in 61..=6204 {
            assert!(system.group_isin("Membrane", i).unwrap());
        }
    }

    #[test]
    fn group_create_macro_water() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.group_create("Water", "@water").unwrap();
        assert!(system.group_exists("Water"));
        assert_eq!(system.group_get_n_atoms("Water").unwrap(), 10399);

        for i in 6205..=16603 {
            assert!(system.group_isin("Water", i).unwrap());
        }
    }

    #[test]
    fn group_create_macro_ion() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.group_create("Ion", "@ion").unwrap();
        assert!(system.group_exists("Ion"));
        assert_eq!(system.group_get_n_atoms("Ion").unwrap(), 240);

        for i in 16604..=16843 {
            assert!(system.group_isin("Ion", i).unwrap());
        }
    }

    #[test]
    fn group_create_macro_dna() {
        let mut system = System::from_file("test_files/protein_with_dna.pdb").unwrap();

        system.group_create("DNA", "@dna").unwrap();
        assert!(system.group_exists("DNA"));
        assert_eq!(system.group_get_n_atoms("DNA").unwrap(), 169);

        for i in 0..=168 {
            assert!(system.group_isin("DNA", i).unwrap());
        }
    }

    #[test]
    fn group_create_macro_rna() {
        let mut system = System::from_file("test_files/rna.pdb").unwrap();

        system.group_create("RNA", "@rna").unwrap();
        assert!(system.group_exists("RNA"));
        assert_eq!(system.group_get_n_atoms("RNA").unwrap(), 1108);

        for i in 0..=1107 {
            assert!(system.group_isin("RNA", i).unwrap());
        }
    }

    #[test]
    fn group_create_complex() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.read_ndx("test_files/index.ndx").unwrap();

        system.group_create(
            "Complex Group", 
            "((serial 1 - 15 or atomid 13 14 15 16 || atomnum 62 64 to 70) && Protein ION) or 
            (resid 11179 to 13000 or resnum 5400) and (resname W or (resname GLY LEU and (name BB or atomname SC1)))").unwrap();
        assert_eq!(system.group_get_n_atoms("Complex Group").unwrap(), 2);

        let mut iter = system.group_iter("Complex Group").unwrap();
        assert_eq!(iter.next().unwrap().get_atom_number(), 1);
        assert_eq!(iter.next().unwrap().get_atom_number(), 11064);
    }

    #[test]
    fn group_create_atom_numbers() {
        let mut system = System::from_file("test_files/example_renumbered.gro").unwrap();

        system.group_create("Serial 6", "serial 6").unwrap();
        system.group_create("Atomid 6", "atomid 6").unwrap();
        system.group_create("Atomnum 6", "atomnum 6").unwrap();

        assert_eq!(
            system.group_get_n_atoms("Serial 6"),
            system.group_get_n_atoms("Atomid 6")
        );
        assert_eq!(
            system.group_get_n_atoms("Serial 6"),
            system.group_get_n_atoms("Atomnum 6")
        );

        assert!(system.group_isin("Serial 6", 5).unwrap());
        assert!(system.group_isin("Atomid 6", 5).unwrap());
        assert!(system.group_isin("Atomnum 6", 5).unwrap());

        system.group_create("Serial 10", "serial 10").unwrap();
        for atom in system.group_iter("Serial 10").unwrap() {
            assert_eq!(atom.get_atom_number(), 49);
        }

        system.group_create("Atomid 49", "atomid 49").unwrap();
        assert_eq!(system.group_get_n_atoms("Atomid 49").unwrap(), 3);

        assert!(system.group_isin("Atomid 49", 9).unwrap());
        assert!(system.group_isin("Atomid 49", 24).unwrap());
        assert!(system.group_isin("Atomid 49", 48).unwrap());

        system.group_create("Atomnum 133", "atomnum 133").unwrap();
        assert_eq!(system.group_get_n_atoms("Atomnum 133").unwrap(), 1);

        system.group_create("Atomid 10", "atomid 10").unwrap();
        assert_eq!(system.group_get_n_atoms("Atomid 10").unwrap(), 0);
    }

    #[test]
    fn group_create_from_geometry_cylinder() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let cylinder = Cylinder::new([5.0, 8.0, 3.0].into(), 2.0, 6.0, Dimension::Y);

        system
            .group_create_from_geometry("Selected Membrane", "Membrane", cylinder.clone())
            .unwrap();

        assert!(system.group_exists("Selected Membrane"));
        assert_eq!(system.group_get_n_atoms("Selected Membrane").unwrap(), 206);

        for atom in system.group_iter("Selected Membrane").unwrap() {
            assert_eq!(atom.get_residue_name(), "POPC");
            assert!(cylinder.inside(atom.get_position(), system.get_box_as_ref()));
        }
    }

    #[test]
    fn group_create_from_geometry_sphere() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let sphere = Sphere::new([0.5, 4.5, 3.5].into(), 4.6);

        system
            .group_create_from_geometry("Selected Water", "resname W", sphere.clone())
            .unwrap();

        assert!(system.group_exists("Selected Water"));
        assert_eq!(system.group_get_n_atoms("Selected Water").unwrap(), 1881);

        for atom in system.group_iter("Selected Water").unwrap() {
            assert_eq!(atom.get_residue_name(), "W");
            assert!(sphere.inside(atom.get_position(), system.get_box_as_ref()));
        }
    }

    #[test]
    fn group_create_from_geometry_rectangular() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let rectangular = Rectangular::new([5.0, 0.0, 2.0].into(), 5.0, 4.0, 4.3);

        system
            .group_create_from_geometry("Selected Protein", "@protein", rectangular.clone())
            .unwrap();

        assert!(system.group_exists("Selected Protein"));
        assert_eq!(system.group_get_n_atoms("Selected Protein").unwrap(), 25);

        for atom in system.group_iter("Selected Protein").unwrap() {
            let resname = atom.get_residue_name();
            assert!(
                resname == "VAL"
                    || resname == "LEU"
                    || resname == "ALA"
                    || resname == "LYS"
                    || resname == "CYS"
            );
            assert!(rectangular.inside(atom.get_position(), system.get_box_as_ref()));
        }
    }

    #[test]
    fn group_create_open_ended_ranges() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.group_create("Group1", "resid < 380").unwrap();

        assert!(system.group_exists("Group1"));
        assert_eq!(system.group_get_n_atoms("Group1").unwrap(), 4261);

        for i in 0..4261 {
            assert!(system.group_isin("Group1", i).unwrap());
        }

        system.group_create("Group2", "resid <= 380").unwrap();
        assert!(system.group_exists("Group2"));
        assert_eq!(system.group_get_n_atoms("Group2").unwrap(), 4273);

        for i in 0..4273 {
            assert!(system.group_isin("Group2", i).unwrap());
        }

        system.group_create("Group3", "serial > 9143").unwrap();
        assert!(system.group_exists("Group3"));
        assert_eq!(system.group_get_n_atoms("Group3").unwrap(), 7701);

        for i in 9143..16844 {
            assert!(system.group_isin("Group3", i).unwrap());
        }

        system.group_create("Group4", "serial >= 9143").unwrap();
        assert!(system.group_exists("Group4"));
        assert_eq!(system.group_get_n_atoms("Group4").unwrap(), 7702);

        for i in 9142..16844 {
            assert!(system.group_isin("Group4", i).unwrap());
        }

        system
            .group_create("Group 5", "serial <= 10000 10005-10010")
            .unwrap();
        assert!(system.group_exists("Group 5"));
        assert_eq!(system.group_get_n_atoms("Group 5").unwrap(), 10006);

        for i in 0..10000 {
            assert!(system.group_isin("Group 5", i).unwrap());
        }

        for i in 10004..10010 {
            assert!(system.group_isin("Group 5", i).unwrap());
        }
    }

    #[test]
    fn group_create_invalid_names() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let names = vec![
            "Invalid&name",
            "Invalid| name",
            "Invalid(name",
            "Invalid)name",
            "Invalidn@me",
            "!nvalidname",
            "Inval'idname",
            "Inval\"idname",
            "    ",
            "Inval=idname",
            "Inval<idname",
            ">Invalidname",
        ];

        for name in names {
            assert_eq!(
                system.group_create(name, "serial 1 to 3"),
                Err(GroupError::InvalidName(name.to_string()))
            );
            assert!(!system.group_exists(name));
        }
    }

    #[test]
    fn group_create_from_indices() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system
            .group_create_from_indices("Test", vec![1, 3, 3, 6, 8])
            .unwrap();

        assert!(system.group_exists("Test"));
        assert_eq!(system.group_get_n_atoms("Test").unwrap(), 4);

        match system.group_create_from_indices("Test", vec![1, 10, 11, 12, 13]) {
            Err(GroupError::AlreadyExistsWarning(e)) => assert_eq!(e, "Test"),
            Ok(_) => panic!("Warning should have been returned, but it was not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }

        assert!(system.group_exists("Test"));
        assert_eq!(system.group_get_n_atoms("Test").unwrap(), 5);
    }

    #[test]
    fn group_create_from_indices_invalid_names() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let names = vec![
            "Invalid&name",
            "Invalid| name",
            "Invalid(name",
            "Invalid)name",
            "Invalidn@me",
            "!nvalidname",
            "Inval'idname",
            "Inval\"idname",
            "    ",
            "Inval=idname",
            "Inval<idname",
            ">Invalidname",
        ];

        for name in names {
            assert_eq!(
                system.group_create_from_indices(name, vec![1, 3, 3, 6, 8]),
                Err(GroupError::InvalidName(name.to_string()))
            );
            assert!(!system.group_exists(name));
        }
    }

    #[test]
    fn group_create_from_ranges() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system
            .group_create_from_ranges("Test", vec![(0, 15), (10, 20)])
            .unwrap();

        assert!(system.group_exists("Test"));
        assert_eq!(system.group_get_n_atoms("Test").unwrap(), 21);

        match system.group_create_from_ranges("Test", vec![(10, 15)]) {
            Err(GroupError::AlreadyExistsWarning(e)) => assert_eq!(e, "Test"),
            Ok(_) => panic!("Warning should have been returned, but it was not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }

        assert!(system.group_exists("Test"));
        assert_eq!(system.group_get_n_atoms("Test").unwrap(), 6);
    }

    #[test]
    fn group_create_from_ranges_invalid_names() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let names = vec![
            "Invalid&name",
            "Invalid| name",
            "Invalid(name",
            "Invalid)name",
            "Invalidn@me",
            "!nvalidname",
            "Inval'idname",
            "Inval\"idname",
            "    ",
            "Inval=idname",
            "Inval<idname",
            ">Invalidname",
        ];

        for name in names {
            assert_eq!(
                system.group_create_from_ranges(name, vec![(0, 15), (10, 20)]),
                Err(GroupError::InvalidName(name.to_string()))
            );
            assert!(!system.group_exists(name));
        }
    }

    #[test]
    fn group_isin_basic() {
        let system = System::from_file("test_files/example.gro").unwrap();

        for i in 0..system.get_n_atoms() {
            assert!(system
                .group_isin("all", i)
                .expect("Group `all` does not exist but it should."));
        }
    }

    #[test]
    fn group_isin_advanced() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.read_ndx("test_files/index.ndx").unwrap();

        for i in 0..=60 {
            assert!(system
                .group_isin("Transmembrane_all", i)
                .expect("Group `Transmembrane_all` does not exist but it should."));
        }

        for i in 61..=6204 {
            assert!(system
                .group_isin("Membrane", i)
                .expect("Group `Membrane` does not exist but it should."));
        }

        for i in 16604..=16843 {
            assert!(system
                .group_isin("ION", i)
                .expect("Group `ION` does not exist but it should."));
        }
    }

    #[test]
    fn group_isin_fails() {
        let system = System::from_file("test_files/example.gro").unwrap();

        assert_eq!(
            system.group_isin("Transmembrane_all", 100),
            Err(GroupError::NotFound("Transmembrane_all".to_string()))
        );
        assert_eq!(
            system.group_isin("Protein", 200),
            Err(GroupError::NotFound("Protein".to_string()))
        );
        assert_eq!(
            system.group_isin("ION", 250),
            Err(GroupError::NotFound("ION".to_string()))
        );
    }

    #[test]
    fn split_by_resid() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let (code, residues) = system.atoms_split_by_resid();
        if let Err(_) = code {
            panic!("Overlapping groups in atoms_split_by_resid.");
        }

        assert_eq!(residues.len(), 11180);
        for groupname in residues.iter() {
            assert!(system.group_exists(groupname));
            assert!(system.group_get_n_atoms(groupname).unwrap() > 0);
        }

        assert_eq!(system.group_get_n_atoms("resid 1").unwrap(), 1);
        assert_eq!(system.group_get_n_atoms("resid 6").unwrap(), 2);
        assert_eq!(system.group_get_n_atoms("resid 27").unwrap(), 3);
        assert_eq!(system.group_get_n_atoms("resid 11180").unwrap(), 1);

        // check order of the group names in the `residues` vector
        for i in 1..=11180 {
            let resname = format!("resid {}", i);
            assert_eq!(residues[i - 1], resname);
        }
    }

    #[test]
    fn split_by_resid_warnings() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system
            .read_ndx("test_files/index_group_by_res.ndx")
            .unwrap();

        let (code, residues) = system.atoms_split_by_resid();
        match code {
            Err(GroupError::MultipleAlreadyExistWarning(e)) => assert_eq!(
                e,
                Box::new(HashSet::from([
                    "resid 6".to_string(),
                    "resid 27".to_string(),
                    "resid 11180".to_string()
                ]))
            ),
            Ok(_) => panic!("Warning should have been returned, but it was not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }

        assert_eq!(residues.len(), 11180);
        for groupname in residues.iter() {
            assert!(system.group_exists(groupname));
            assert!(system.group_get_n_atoms(groupname).unwrap() > 0);
        }

        assert_eq!(system.group_get_n_atoms("resid 1").unwrap(), 1);
        assert_eq!(system.group_get_n_atoms("resid 6").unwrap(), 2);
        assert_eq!(system.group_get_n_atoms("resid 27").unwrap(), 3);
        assert_eq!(system.group_get_n_atoms("resid 11180").unwrap(), 1);
    }

    #[test]
    fn split_by_resid_broken() {
        let mut system = System::from_file("test_files/example_shuffled_residues.gro").unwrap();

        let (code, residues) = system.atoms_split_by_resid();
        if let Err(_) = code {
            panic!("Overlapping groups in atoms_split_by_resid.");
        }

        assert_eq!(residues.len(), 21);
        for groupname in residues.iter() {
            assert!(system.group_exists(groupname));
            assert!(system.group_get_n_atoms(groupname).unwrap() > 0);
        }

        let expected_n = vec![
            2, 3, 2, 2, 3, 2, 1, 2, 2, 4, 2, 2, 1, 3, 2, 4, 3, 2, 2, 3, 3,
        ];
        for i in 1..=21 {
            let resname = format!("resid {}", i);
            assert_eq!(
                expected_n[i - 1],
                system.group_get_n_atoms(&resname).unwrap()
            );
        }

        // check order of the group names in the `residues` vector
        assert_eq!(residues[0], "resid 20");
        assert_eq!(residues[1], "resid 1");
        assert_eq!(residues[20], "resid 21");
    }

    #[test]
    fn group_split_by_resid() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.read_ndx("test_files/index.ndx").unwrap();

        let (code, residues) = system.group_split_by_resid("Protein");
        if let Err(_) = code {
            panic!("Error occured when using `System::group_split_by_resid`");
        }

        assert_eq!(residues.len(), 29);
        for groupname in residues.iter() {
            assert!(system.group_exists(groupname));
            assert!(system.group_get_n_atoms(groupname).unwrap() > 0);
        }

        assert_eq!(system.group_get_n_atoms("resid 1").unwrap(), 1);
        assert_eq!(system.group_get_n_atoms("resid 2").unwrap(), 3);
        assert_eq!(system.group_get_n_atoms("resid 15").unwrap(), 2);
        assert_eq!(system.group_get_n_atoms("resid 29").unwrap(), 2);

        // check order of the group names in the `residues` vector
        for i in 1..=29 {
            let resname = format!("resid {}", i);
            assert_eq!(residues[i - 1], resname);
        }
    }

    #[test]
    fn group_split_by_resid_warnings() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.read_ndx("test_files/index.ndx").unwrap();
        system
            .read_ndx("test_files/index_group_by_res.ndx")
            .unwrap();

        let (code, residues) = system.group_split_by_resid("Protein");
        match code {
            Err(GroupError::MultipleAlreadyExistWarning(e)) => assert_eq!(
                e,
                Box::new(HashSet::from([
                    "resid 6".to_string(),
                    "resid 27".to_string(),
                ]))
            ),
            Ok(_) => panic!("Warning should have been returned, but it was not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }

        assert_eq!(residues.len(), 29);
        for groupname in residues.iter() {
            assert!(system.group_exists(groupname));
            assert!(system.group_get_n_atoms(groupname).unwrap() > 0);
        }

        assert_eq!(system.group_get_n_atoms("resid 1").unwrap(), 1);
        assert_eq!(system.group_get_n_atoms("resid 2").unwrap(), 3);
        assert_eq!(system.group_get_n_atoms("resid 15").unwrap(), 2);
        assert_eq!(system.group_get_n_atoms("resid 29").unwrap(), 2);

        // check order of the group names in the `residues` vector
        for i in 1..=29 {
            let resname = format!("resid {}", i);
            assert_eq!(residues[i - 1], resname);
        }
    }

    #[test]
    fn group_split_by_resid_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let (code, residues) = system.group_split_by_resid("Protein");
        match code {
            Err(GroupError::NotFound(e)) => assert_eq!(e, "Protein"),
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }

        assert_eq!(residues.len(), 0);
    }

    #[test]
    fn split_by_resname() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let (code, residues) = system.atoms_split_by_resname();
        if let Err(_) = code {
            panic!("Overlapping groups in atoms_split_by_resname.");
        }

        assert_eq!(residues.len(), 9);
        for groupname in residues.iter() {
            assert!(system.group_exists(groupname));
            assert!(system.group_get_n_atoms(groupname).unwrap() > 0);
        }

        assert_eq!(system.group_get_n_atoms("resname GLY").unwrap(), 1);
        assert_eq!(system.group_get_n_atoms("resname LYS").unwrap(), 12);
        assert_eq!(system.group_get_n_atoms("resname VAL").unwrap(), 22);
        assert_eq!(system.group_get_n_atoms("resname LEU").unwrap(), 2);
        assert_eq!(system.group_get_n_atoms("resname ALA").unwrap(), 22);
        assert_eq!(system.group_get_n_atoms("resname CYS").unwrap(), 2);
        assert_eq!(system.group_get_n_atoms("resname POPC").unwrap(), 6144);
        assert_eq!(system.group_get_n_atoms("resname W").unwrap(), 10399);
        assert_eq!(system.group_get_n_atoms("resname ION").unwrap(), 240);
    }

    #[test]
    fn split_by_resname_warnings() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system
            .read_ndx("test_files/index_group_by_resname.ndx")
            .unwrap();

        let (code, residues) = system.atoms_split_by_resname();
        match code {
            Err(GroupError::MultipleAlreadyExistWarning(e)) => assert_eq!(
                e,
                Box::new(HashSet::from([
                    "resname POPC".to_string(),
                    "resname LYS".to_string(),
                    "resname W".to_string()
                ]))
            ),
            Ok(_) => panic!("Warning should have been returned, but it was not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }

        assert_eq!(residues.len(), 9);
        for groupname in residues.iter() {
            assert!(system.group_exists(groupname));
            assert!(system.group_get_n_atoms(groupname).unwrap() > 0);
        }

        assert_eq!(system.group_get_n_atoms("resname GLY").unwrap(), 1);
        assert_eq!(system.group_get_n_atoms("resname LYS").unwrap(), 12);
        assert_eq!(system.group_get_n_atoms("resname VAL").unwrap(), 22);
        assert_eq!(system.group_get_n_atoms("resname LEU").unwrap(), 2);
        assert_eq!(system.group_get_n_atoms("resname ALA").unwrap(), 22);
        assert_eq!(system.group_get_n_atoms("resname CYS").unwrap(), 2);
        assert_eq!(system.group_get_n_atoms("resname POPC").unwrap(), 6144);
        assert_eq!(system.group_get_n_atoms("resname W").unwrap(), 10399);
        assert_eq!(system.group_get_n_atoms("resname ION").unwrap(), 240);
    }

    #[test]
    fn group_split_by_resname() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.read_ndx("test_files/index.ndx").unwrap();

        let (code, residues) = system.group_split_by_resname("Protein");
        if let Err(_) = code {
            panic!("Error occured when using `System::group_split_by_resname`");
        }

        assert_eq!(residues.len(), 6);
        for groupname in residues.iter() {
            assert!(system.group_exists(groupname));
            assert!(system.group_get_n_atoms(groupname).unwrap() > 0);
        }

        assert_eq!(system.group_get_n_atoms("resname GLY").unwrap(), 1);
        assert_eq!(system.group_get_n_atoms("resname LYS").unwrap(), 12);
        assert_eq!(system.group_get_n_atoms("resname VAL").unwrap(), 22);
        assert_eq!(system.group_get_n_atoms("resname LEU").unwrap(), 2);
        assert_eq!(system.group_get_n_atoms("resname ALA").unwrap(), 22);
        assert_eq!(system.group_get_n_atoms("resname CYS").unwrap(), 2);
    }

    #[test]
    fn group_split_by_resname_warnings() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.read_ndx("test_files/index.ndx").unwrap();
        system
            .read_ndx("test_files/index_group_by_resname.ndx")
            .unwrap();

        let (code, residues) = system.group_split_by_resname("Protein");
        match code {
            Err(GroupError::MultipleAlreadyExistWarning(e)) => {
                assert_eq!(e, Box::new(HashSet::from(["resname LYS".to_string(),])))
            }
            Ok(_) => panic!("Warning should have been returned, but it was not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }

        assert_eq!(residues.len(), 6);
        for groupname in residues.iter() {
            assert!(system.group_exists(groupname));
            assert!(system.group_get_n_atoms(groupname).unwrap() > 0);
        }

        assert_eq!(system.group_get_n_atoms("resname GLY").unwrap(), 1);
        assert_eq!(system.group_get_n_atoms("resname LYS").unwrap(), 12);
        assert_eq!(system.group_get_n_atoms("resname VAL").unwrap(), 22);
        assert_eq!(system.group_get_n_atoms("resname LEU").unwrap(), 2);
        assert_eq!(system.group_get_n_atoms("resname ALA").unwrap(), 22);
        assert_eq!(system.group_get_n_atoms("resname CYS").unwrap(), 2);
    }

    #[test]
    fn group_split_by_resname_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let (code, residues) = system.group_split_by_resname("Protein");
        match code {
            Err(GroupError::NotFound(e)) => assert_eq!(e, "Protein"),
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }

        assert_eq!(residues.len(), 0);
    }

    #[test]
    fn group_make_writable() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert!(!system.get_groups_as_ref().get("all").unwrap().print_ndx);
        system.group_make_writable("all").unwrap();
        assert!(system.get_groups_as_ref().get("all").unwrap().print_ndx);
    }

    #[test]
    fn group_make_writable_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.group_make_writable("Protein") {
            Ok(_) => panic!("Group should not exist."),
            Err(e) => assert_eq!(e, GroupError::NotFound(String::from("Protein"))),
        }
    }

    #[test]
    fn group_make_nonwritable() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.group_create("Group", "resid 1").unwrap();

        assert!(system.get_groups_as_ref().get("Group").unwrap().print_ndx);
        system.group_make_nonwritable("Group").unwrap();
        assert!(!system.get_groups_as_ref().get("Group").unwrap().print_ndx);
    }

    #[test]
    fn group_make_nonwritable_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.group_make_nonwritable("Protein") {
            Ok(_) => panic!("Group should not exist."),
            Err(e) => assert_eq!(e, GroupError::NotFound(String::from("Protein"))),
        }
    }

    #[test]
    fn group_union_simple() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.group_union("Protein", "W", "Protein_W").unwrap();

        assert!(system.group_exists("Protein_W"));
        assert_eq!(
            system.group_get_n_atoms("Protein_W").unwrap(),
            system.group_get_n_atoms("Protein").unwrap() + system.group_get_n_atoms("W").unwrap()
        );
    }

    #[test]
    fn group_union_overlap() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system
            .group_union("Membrane", "System", "Membrane_System")
            .unwrap();

        assert!(system.group_exists("Membrane_System"));
        assert_eq!(
            system.group_get_n_atoms("Membrane_System").unwrap(),
            system.group_get_n_atoms("System").unwrap()
        );
    }

    #[test]
    fn group_union_already_exists() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        match system.group_union("Protein", "ION", "W_ION") {
            Err(GroupError::AlreadyExistsWarning(e)) => assert_eq!(e, "W_ION"),
            Ok(_) => panic!("Warning should have been returned, but it was not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }

        assert_eq!(
            system.group_get_n_atoms("W_ION").unwrap(),
            system.group_get_n_atoms("Protein").unwrap() + system.group_get_n_atoms("ION").unwrap()
        );
    }

    #[test]
    fn group_union_invalid_group1() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        match system.group_union("Nonexistent", "ION", "Test") {
            Err(GroupError::NotFound(e)) => assert_eq!(e, "Nonexistent"),
            Ok(_) => panic!("Creating union should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }

        assert!(!system.group_exists("Test"));
    }

    #[test]
    fn group_union_invalid_group2() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        match system.group_union("ION", "Nonexistent", "Test") {
            Err(GroupError::NotFound(e)) => assert_eq!(e, "Nonexistent"),
            Ok(_) => panic!("Creating union should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }

        assert!(!system.group_exists("Test"));
    }

    #[test]
    fn group_extend_simple() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let old_n_atoms = system.group_get_n_atoms("Protein").unwrap();

        system.group_extend("Protein", "ION").unwrap();
        assert_eq!(
            system.group_get_n_atoms("Protein").unwrap(),
            old_n_atoms + system.group_get_n_atoms("ION").unwrap()
        );
    }

    #[test]
    fn group_extend_overlap() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let old_n_atoms = system.group_get_n_atoms("Protein_Membrane").unwrap();

        system.group_extend("Protein_Membrane", "Membrane").unwrap();
        assert_eq!(
            system.group_get_n_atoms("Protein_Membrane").unwrap(),
            old_n_atoms
        );
    }

    #[test]
    fn group_names() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.group_create("Custom Group", "resname LYS").unwrap();

        let names = system.group_names();

        assert_eq!(names.len(), 24);

        for name in names {
            assert!(system.group_exists(&name));
        }
    }

    #[test]
    fn group_names_writable() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.group_create("Custom Group", "resname LYS").unwrap();

        let names = system.group_names_writable();

        assert_eq!(names.len(), 22);

        for name in &names {
            assert!(system.group_exists(name));
        }

        assert!(!names.contains(&"All".to_string()));
        assert!(!names.contains(&"all".to_string()));

        // set "Custom Group" to nonwritable
        system.group_make_nonwritable("Custom Group").unwrap();

        let new_names = system.group_names_writable();

        assert_eq!(new_names.len(), 21);
        assert!(!new_names.contains(&"Custom Group".to_string()));
    }

    #[test]
    fn group_remove() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        unsafe {
            system.group_remove("Protein").unwrap();
        }

        assert!(!system.group_exists("Protein"));
    }

    #[test]
    fn group_remove_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        unsafe {
            match system.group_remove("Proin") {
                Err(GroupError::NotFound(e)) => assert_eq!("Proin", e),
                Ok(_) => panic!("Function should have failed, but it was successful."),
                Err(e) => panic!("Incorrect error type returned: {:?}", e),
            }
        }
    }

    #[test]
    fn group_rename() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        unsafe {
            system.group_rename("Protein", "My Protein Group").unwrap();
        }

        assert!(!system.group_exists("Protein"));
        assert!(system.group_exists("My Protein Group"));

        for index in 0..61 {
            assert!(system.group_isin("My Protein Group", index).unwrap());
        }
    }

    #[test]
    fn group_rename_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        unsafe {
            match system.group_rename("Proin", "Protein") {
                Err(GroupError::NotFound(e)) => assert_eq!("Proin", e),
                Ok(_) => panic!("Function should have failed, but it was successful."),
                Err(e) => panic!("Incorrect error type returned: {:?}", e),
            }
        }
    }

    #[test]
    fn group_rename_overwrite() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        unsafe {
            match system.group_rename("Protein", "Membrane") {
                Err(GroupError::AlreadyExistsWarning(e)) => assert_eq!("Membrane", e),
                Ok(_) => panic!("Function should have raised warning, but it did not."),
                Err(e) => panic!("Incorrect error type returned: {:?}", e),
            }
        }

        assert!(!system.group_exists("Protein"));
        assert!(system.group_exists("Membrane"));

        for index in 0..61 {
            assert!(system.group_isin("Membrane", index).unwrap());
        }
    }
}