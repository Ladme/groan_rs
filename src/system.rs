// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of the System structure and its basic methods.

use indexmap::IndexMap;
use std::collections::HashSet;
use std::error::Error;
use std::path::Path;

use crate::atom::Atom;
use crate::errors::{AtomError, GroupError, ParseFileError};
use crate::files::FileType;
use crate::gro_io;
use crate::group::Group;
use crate::pdb_io;
use crate::simbox::SimBox;
use crate::vector3d::Vector3D;

#[derive(Debug)]
pub struct System {
    /// Name of the molecular system.
    name: String,
    /// Vector of atoms in the system.
    atoms: Vec<Atom>,
    /// Size of the simulation box.
    simulation_box: SimBox,
    /// Groups of atoms associated with the system.
    groups: IndexMap<String, Group>,
    /// Current simulation step.
    simulation_step: u64,
    /// Current simulation time.
    simulation_time: f32,
    /// Precision of the coordinates.
    coordinates_precision: u64,
    /// Lambda
    lambda: f32,
}

impl System {
    /// Create new System structure with a given name from the provided vector of atoms and simulation box.
    ///
    /// ## Notes
    /// - The returned `System` structure will contain two default groups "all" and "All",
    /// each consisting of all the atoms in the system.
    ///
    /// ## Example 1: Manually creating a system
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let name = "My System";
    /// let atoms = Vec::new();
    ///
    /// // ... fill the `atoms` vector with Atom structures ...
    ///
    /// let simulation_box = SimBox::from([10.0, 10.0, 12.0]);
    ///
    /// // construct the molecular system
    /// let system = System::new(name, atoms, simulation_box);
    /// ```
    ///
    /// ## Example 2: Creating system from other system using `extract`
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// // load system from file
    /// let mut original_system = System::from_file("system.gro").unwrap();
    /// // create a group "Protein" consisting of atoms of residues 1 to 29
    /// original_system.group_create("Protein", "resid 1 to 29").unwrap();
    ///
    /// // extract atoms from group "Protein"
    /// let protein = original_system.group_extract("Protein").unwrap();
    /// // create a new system containing only "Protein" atoms
    /// let new_system = System::new(
    ///     "System containing protein atoms only",
    ///     protein,
    ///     original_system.get_box_copy());
    /// ```
    ///
    /// ## Example 3: Creating system from other system using iterators
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut original_system = System::from_file("system.gro").unwrap();
    ///
    /// // construct a sphere located at x = 1, y = 2, z = 3 with a radius of 2.5 nm
    /// let sphere = Sphere::new([1.0, 2.0, 3.0].into(), 2.5);
    ///
    /// // create iterator over the atoms of the system
    /// // only select atoms which are inside the above-defined sphere
    /// let iterator = original_system
    ///     .atoms_iter()
    ///     .filter_geometry(sphere);
    ///
    /// let new_system = System::new(
    ///     "System containing atoms located inside the sphere",
    ///     iterator.cloned().collect(),
    ///     original_system.get_box_copy());
    /// ```
    pub fn new(name: &str, atoms: Vec<Atom>, simulation_box: SimBox) -> Self {
        let mut system = System {
            name: name.to_string(),
            atoms,
            simulation_box,
            groups: IndexMap::new(),
            simulation_step: 0u64,
            simulation_time: 0.0f32,
            coordinates_precision: 100u64,
            lambda: 0.0,
        };

        match system.group_create_all() {
            Err(_) => {
                panic!("Groan error. Group `all` or `All` already exists as System is created.")
            }
            Ok(_) => system,
        }
    }

    /// Create a new System from gro file or pdb file.
    /// The method will attempt to automatically recognize gro or pdb file based on the file extension.
    ///
    /// ## Returns
    /// `System` structure if successful.
    /// `ParseFileError` if the file format is not supported.
    /// `ParseGroError` if parsing of the gro file fails.
    /// `ParsePdbError` if parsing of the pdb file fails.
    ///
    /// ## Example
    /// Reading gro file.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let system = match System::from_file("system.gro") {
    ///     Ok(x) => x,
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;
    ///     }
    /// };
    /// ```
    /// ## Notes
    /// - The returned System structure will contain two default groups "all" and "All"
    /// consisting of all the atoms in the system.
    pub fn from_file(filename: impl AsRef<Path>) -> Result<Self, Box<dyn Error>> {
        match FileType::from_name(&filename) {
            FileType::GRO => gro_io::read_gro(filename).map_err(Box::from),
            FileType::PDB => pdb_io::read_pdb(filename).map_err(Box::from),
            _ => Err(Box::from(ParseFileError::UnknownExtension(Box::from(
                filename.as_ref(),
            )))),
        }
    }

    /**************************/
    /*  ACCESSING PROPERTIES  */
    /**************************/

    /// Get the name of the molecular system.
    pub fn get_name(&self) -> &str {
        &self.name
    }

    /// Get immutable reference to the atoms in the system.
    pub fn get_atoms_as_ref(&self) -> &Vec<Atom> {
        &self.atoms
    }

    /// Get mutable reference to the atoms in the system.
    ///
    /// ## Safety
    /// - This function is unsafe as manually changing the `atoms` of the system
    /// can cause the system to become invalid. Other functions may then not work correctly.
    /// - Notably, no atoms can be added or removed from the `atoms` vector as such
    /// operation would make all the groups associated with the system invalid. The same goes
    /// for reordering the atoms.
    /// - The properties of the individual atoms can however be safely changed.
    pub unsafe fn get_atoms_as_ref_mut(&mut self) -> &mut Vec<Atom> {
        &mut self.atoms
    }

    /// Get copy of the atoms in the system.
    pub fn get_atoms_copy(&self) -> Vec<Atom> {
        self.atoms.clone()
    }

    /// Get immutable reference to the groups in the system.
    pub fn get_groups_as_ref(&self) -> &IndexMap<String, Group> {
        &self.groups
    }

    /// Get mutable reference to the groups in the system.
    ///
    /// ## Safety
    /// - This function is unsafe as manually changing the `groups` of the system
    /// can cause the system to become invalid.
    pub unsafe fn get_groups_as_ref_mut(&mut self) -> &mut IndexMap<String, Group> {
        &mut self.groups
    }

    /// Get copy of the groups in the system.
    pub fn get_groups_copy(&self) -> IndexMap<String, Group> {
        self.groups.clone()
    }

    /// Get immutable reference to the simulation box.
    pub fn get_box_as_ref(&self) -> &SimBox {
        &self.simulation_box
    }

    /// Get center of the simulation box.
    ///
    /// ## Warning
    /// Currently only implemented for orthogonal simulation boxes!
    pub fn get_box_center(&self) -> Vector3D {
        Vector3D {
            x: self.simulation_box.x / 2.0f32,
            y: self.simulation_box.y / 2.0f32,
            z: self.simulation_box.z / 2.0f32,
        }
    }

    /// Get mutable reference to the simulation box.
    pub fn get_box_as_ref_mut(&mut self) -> &mut SimBox {
        &mut self.simulation_box
    }

    /// Get copy of the simulation box.
    pub fn get_box_copy(&self) -> SimBox {
        self.simulation_box.clone()
    }

    /// Get the number of atoms in the system.
    pub fn get_n_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Get the number of groups in the system. This counts all groups, even the default ones.
    pub fn get_n_groups(&self) -> usize {
        self.groups.len()
    }

    /// Get the current simulation time.
    pub fn get_simulation_time(&self) -> f32 {
        self.simulation_time
    }

    /// Get the current simulation step.
    pub fn get_simulation_step(&self) -> u64 {
        self.simulation_step
    }

    /// Get the precision of the coordinates.
    pub fn get_precision(&self) -> u64 {
        self.coordinates_precision
    }

    /// Get the simulation lambda.
    pub fn get_lambda(&self) -> f32 {
        self.lambda
    }

    /// Set the simulation time.
    pub fn set_simulation_time(&mut self, time: f32) {
        self.simulation_time = time;
    }

    /// Set the simulation step.
    pub fn set_simulation_step(&mut self, step: u64) {
        self.simulation_step = step;
    }

    /// Set simulation box.
    pub fn set_box(&mut self, sim_box: SimBox) {
        self.simulation_box = sim_box;
    }

    /// Set precision of the coordinates.
    pub fn set_precision(&mut self, precision: u64) {
        self.coordinates_precision = precision;
    }

    /// Set the simulation lambda.
    pub fn set_lambda(&mut self, lambda: f32) {
        self.lambda = lambda;
    }

    /// Check whether velocities are present.
    ///
    /// ## Returns
    /// `true` if any of the atoms in the system has non-zero velocity. `false` otherwise.
    ///
    /// ## Notes
    /// - Complexity of this operation is O(n), where n is the number of atoms in the system.
    pub fn has_velocities(&self) -> bool {
        self.atoms.iter().any(|atom| atom.has_velocity())
    }

    /// Check whether forces are present.
    ///
    /// ## Returns
    /// `true` if any of the atoms in the system has non-zero force acting on it. `false` otherwise.
    ///
    /// ## Notes
    /// - Complexity of this operation is O(n), where n is the number of atoms in the system.
    pub fn has_forces(&self) -> bool {
        self.atoms.iter().any(|atom| atom.has_force())
    }

    /**************************/
    /*    CREATING GROUPS     */
    /**************************/

    /// Create two groups each containing all atoms in the system: "all" and "All".
    ///
    /// ## Returns
    /// - `Ok` if both groups were created or GroupError in case any group with the same name already exists.
    fn group_create_all(&mut self) -> Result<(), GroupError> {
        self.group_create_from_ranges("all", vec![(0, self.get_n_atoms())])?;
        self.group_create_from_ranges("All", vec![(0, self.get_n_atoms())])?;

        self.groups
            .get_mut("all")
            .expect("Groan error. Group `all` is not available after creating it.")
            .print_ndx = false;

        self.groups
            .get_mut("All")
            .expect("Groan error. Group `All` is not available after creating it.")
            .print_ndx = false;

        Ok(())
    }

    /// Make a group with a given name from the given Groan selection language query.
    ///
    /// ## Returns
    /// - `Ok` if the group was successfully created.
    /// - `GroupError::AlreadyExistsWarning` if the new group has overwritten a previous group.
    /// - `GroupError::InvalidName` if the name of the group is invalid (no group created).
    /// - `SelectError` if the query could not be parsed.
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
    /// - The following characters are not allowed in group names: '"&|!@()
    /// - The group will be created even if the query selects no atoms.
    pub fn group_create(&mut self, name: &str, query: &str) -> Result<(), Box<dyn Error>> {
        if !Group::name_is_valid(name) {
            return Err(Box::from(GroupError::InvalidName(name.to_string())));
        }

        let group = Group::from_query(query, self)?;

        match self.groups.insert(name.to_string(), group) {
            None => Ok(()),
            Some(_) => Err(Box::from(GroupError::AlreadyExistsWarning(
                name.to_string(),
            ))),
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
    /// - The following characters are not allowed in group names: '"&|!@()
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
        match self.groups.insert(name.to_string(), group) {
            None => Ok(()),
            Some(_) => Err(GroupError::AlreadyExistsWarning(name.to_string())),
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
    /// - The following characters are not allowed in group names: '"&|!@()
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
        match self.groups.insert(name.to_string(), group) {
            None => Ok(()),
            Some(_) => Err(GroupError::AlreadyExistsWarning(name.to_string())),
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
    /// let (result, residues) = system.group_by_resid();
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
    pub fn group_by_resid(&mut self) -> (Result<(), GroupError>, Vec<String>) {
        // split all atoms into residues by their number
        let mut groups: IndexMap<String, Vec<usize>> = IndexMap::new();

        for (atomid, atom) in self.atoms.iter().enumerate() {
            let res = atom.get_residue_number();
            let groupname = format!("resid {}", res);

            match groups.get_mut(&groupname) {
                Some(group) => group.push(atomid),
                None => {
                    groups.insert(groupname, vec![atomid]);
                }
            }
        }

        let mut duplicate_names: HashSet<String> = HashSet::new();

        // create groups
        for (groupname, atoms) in groups.iter() {
            if self
                .group_create_from_indices(groupname, atoms.clone())
                .is_err()
            {
                duplicate_names.insert(groupname.to_string());
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
    /// let (result, residues) = system.group_by_resname();
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
    pub fn group_by_resname(&mut self) -> (Result<(), GroupError>, Vec<String>) {
        // split all atoms into residues by their name
        let mut groups: IndexMap<String, Vec<usize>> = IndexMap::new();

        for (atomid, atom) in self.atoms.iter().enumerate() {
            let res = atom.get_residue_name();
            let groupname = format!("resname {}", res);

            match groups.get_mut(&groupname) {
                Some(group) => group.push(atomid),
                None => {
                    groups.insert(groupname, vec![atomid]);
                }
            }
        }

        let mut duplicate_names: HashSet<String> = HashSet::new();

        // create groups
        for (groupname, atoms) in groups.iter() {
            if self
                .group_create_from_indices(groupname, atoms.clone())
                .is_err()
            {
                duplicate_names.insert(groupname.to_string());
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

    /**************************/
    /* OPERATIONS WITH GROUPS */
    /**************************/

    /// Check whether a group with a given name exists in the system.
    pub fn group_exists(&self, name: &str) -> bool {
        self.groups.contains_key(name)
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
            .groups
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
            .groups
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
            .groups
            .get(group1)
            .ok_or(GroupError::NotFound(group1.to_string()))?;
        let group2 = self
            .groups
            .get(group2)
            .ok_or(GroupError::NotFound(group2.to_string()))?;

        let mut atom_ranges =
            Vec::with_capacity(group1.atom_ranges.len() + group2.atom_ranges.len());
        atom_ranges.extend(group1.atom_ranges.iter());
        atom_ranges.extend(group2.atom_ranges.iter());

        let group = Group::from_ranges(atom_ranges, self.get_n_atoms());
        match self.groups.insert(union.to_string(), group) {
            None => Ok(()),
            Some(_) => Err(GroupError::AlreadyExistsWarning(union.to_string())),
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

    /// Get all group names associated with the system including default groups.
    pub fn group_names(&self) -> Vec<String> {
        let mut names = Vec::with_capacity(self.groups.len());

        for key in self.groups.keys() {
            names.push(key.to_owned());
        }

        names
    }

    /**************************/
    /* ITERATING & EXTRACTING */
    /**************************/

    /// Copy the atoms in the system into an independent vector.
    /// Same as [`get_atoms_copy`].
    ///
    /// ## Example
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let system = System::from_file("system.gro").unwrap();
    /// let extracted: Vec<Atom> = system.atoms_extract();
    /// ```
    /// [`get_atoms_copy`]: System::get_atoms_copy
    pub fn atoms_extract(&self) -> Vec<Atom> {
        self.atoms.clone()
    }

    /// Copy the atoms in a group into an independent vector.
    ///
    /// ## Returns
    /// A vector containing copies of the atoms in the group.
    /// `GroupError::NotFound` if the group does not exist.
    ///
    /// ## Example
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// let extracted_group: Vec<Atom> = match system.group_extract("Protein") {
    ///     Ok(x) => x,
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;
    ///     }
    /// };
    /// ```
    pub fn group_extract(&self, name: &str) -> Result<Vec<Atom>, GroupError> {
        Ok(self.group_iter(name)?.cloned().collect())
    }

    /// Get immutable reference to an atom with target atom number.
    ///
    /// ## Returns
    /// Reference to `Atom` structure or `AtomError::OutOfRange` if `gmx_number` is out of range.
    pub fn get_atom_as_ref(&self, gmx_number: usize) -> Result<&Atom, AtomError> {
        if gmx_number == 0 || gmx_number > self.atoms.len() {
            return Err(AtomError::OutOfRange(gmx_number));
        }

        Ok(&self.atoms[gmx_number - 1])
    }

    /// Get mutable reference to an atom with target atom number.
    ///
    /// ## Returns
    /// Mutable reference to `Atom` structure or `AtomError::OutOfRange` if `gmx_number` is out of range.
    pub fn get_atom_as_ref_mut(&mut self, gmx_number: usize) -> Result<&mut Atom, AtomError> {
        if gmx_number == 0 || gmx_number > self.atoms.len() {
            return Err(AtomError::OutOfRange(gmx_number));
        }

        Ok(&mut self.atoms[gmx_number - 1])
    }

    /// Get copy of an atom with target atom number.
    ///
    /// ## Returns
    /// Copy of an `Atom` structure or `AtomError::OutOfRange` if `gmx_number` is out of range
    pub fn get_atom_copy(&self, gmx_number: usize) -> Result<Atom, AtomError> {
        if gmx_number == 0 || gmx_number > self.atoms.len() {
            return Err(AtomError::OutOfRange(gmx_number));
        }

        Ok(self.atoms[gmx_number - 1].clone())
    }

    /**************************/
    /*    MODIFYING SYSTEM    */
    /**************************/

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
            let simbox = &self.simulation_box as *const SimBox;

            for atom in self.group_iter_mut(name)? {
                atom.translate(
                    vector,
                    simbox
                        .as_ref()
                        .expect("Groan error. SimBox is NULL which is impossible."),
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
        for atom in self.atoms.iter_mut() {
            atom.translate(vector, &self.simulation_box);
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
        for (i, atom) in self.atoms.iter_mut().enumerate() {
            atom.set_atom_number(i + 1);
        }
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
    fn new() {
        let system = System::new(
            "System generated using the `groan_rs` library.",
            Vec::new(),
            [1.5, 3.3, 0.8].into(),
        );

        assert_eq!(
            system.get_name(),
            "System generated using the `groan_rs` library."
        );
        assert_eq!(system.get_atoms_as_ref().len(), 0);

        assert_approx_eq!(f32, system.get_box_as_ref().v1x, 1.5f32);
        assert_approx_eq!(f32, system.get_box_as_ref().v2y, 3.3f32);
        assert_approx_eq!(f32, system.get_box_as_ref().v3z, 0.8f32);
        assert_eq!(system.get_box_as_ref().v1y, 0.0f32);
        assert_eq!(system.get_box_as_ref().v1z, 0.0f32);
        assert_eq!(system.get_box_as_ref().v2x, 0.0f32);
        assert_eq!(system.get_box_as_ref().v2z, 0.0f32);
        assert_eq!(system.get_box_as_ref().v3x, 0.0f32);
        assert_eq!(system.get_box_as_ref().v3y, 0.0f32);

        assert!(system.group_exists("all"));
    }

    #[test]
    fn from_file() {
        let system_gro = System::from_file("test_files/example_novelocities.gro").unwrap();

        assert_eq!(system_gro.get_name(), "Buforin II peptide P11L");
        assert_eq!(system_gro.get_n_atoms(), 50);

        let simbox = system_gro.get_box_as_ref();
        assert_approx_eq!(f32, simbox.x, 6.08608);
        assert_approx_eq!(f32, simbox.y, 6.08608);
        assert_approx_eq!(f32, simbox.z, 6.08608);

        assert_eq!(simbox.v1y, 0.0f32);
        assert_eq!(simbox.v1z, 0.0f32);
        assert_eq!(simbox.v2x, 0.0f32);

        assert_eq!(simbox.v2z, 0.0f32);
        assert_eq!(simbox.v3x, 0.0f32);
        assert_eq!(simbox.v3y, 0.0f32);

        let system_pdb = System::from_file("test_files/example.pdb").unwrap();
        assert_eq!(system_pdb.get_name(), "Buforin II peptide P11L");
        assert_eq!(system_pdb.get_n_atoms(), 50);

        let simbox = system_pdb.get_box_as_ref();
        assert_approx_eq!(f32, simbox.x, 6.0861);
        assert_approx_eq!(f32, simbox.y, 6.0861);
        assert_approx_eq!(f32, simbox.z, 6.0861);

        assert_eq!(simbox.v1y, 0.0f32);
        assert_eq!(simbox.v1z, 0.0f32);
        assert_eq!(simbox.v2x, 0.0f32);

        assert_eq!(simbox.v2z, 0.0f32);
        assert_eq!(simbox.v3x, 0.0f32);
        assert_eq!(simbox.v3y, 0.0f32);

        // compare atoms from PDB an GRO file
        for (groa, pdba) in system_gro.atoms_iter().zip(system_pdb.atoms_iter()) {
            assert_eq!(groa.get_residue_number(), pdba.get_residue_number());
            assert_eq!(groa.get_residue_name(), pdba.get_residue_name());
            assert_eq!(groa.get_atom_number(), pdba.get_atom_number());
            assert_eq!(groa.get_atom_name(), pdba.get_atom_name());
            assert_approx_eq!(f32, groa.get_position().x, pdba.get_position().x);
            assert_approx_eq!(f32, groa.get_position().y, pdba.get_position().y);
            assert_approx_eq!(f32, groa.get_position().z, pdba.get_position().z);

            assert_eq!(groa.get_velocity(), pdba.get_velocity());
            assert_eq!(groa.get_force(), pdba.get_force());
        }
    }

    #[test]
    fn from_file_unknown() {
        match System::from_file("test_files/index.ndx") {
            Ok(_) => panic!("Parsing should have failed."),
            Err(e) => assert!(e.to_string().contains("test_files/index.ndx")),
        }
    }

    #[test]
    fn from_file_no_extension() {
        match System::from_file("LICENSE") {
            Ok(_) => panic!("Parsing should have failed."),
            Err(e) => assert!(e.to_string().contains("LICENSE")),
        }
    }

    #[test]
    fn get_n_atoms() {
        let system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(system.get_n_atoms(), 16844);
    }

    #[test]
    fn test_get_name() {
        let system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(
            system.get_name(),
            "INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1"
        );
    }

    #[test]
    fn get_box_copy() {
        let system = System::from_file("test_files/example_box9.gro").unwrap();

        let mut simbox = system.get_box_copy();

        assert_approx_eq!(f32, simbox.x, system.get_box_as_ref().x);
        assert_approx_eq!(f32, simbox.y, system.get_box_as_ref().y);
        assert_approx_eq!(f32, simbox.z, system.get_box_as_ref().z);

        assert_approx_eq!(f32, simbox.v1y, system.get_box_as_ref().v1y);
        assert_approx_eq!(f32, simbox.v1z, system.get_box_as_ref().v1z);
        assert_approx_eq!(f32, simbox.v2x, system.get_box_as_ref().v2x);

        assert_approx_eq!(f32, simbox.v2z, system.get_box_as_ref().v2z);
        assert_approx_eq!(f32, simbox.v3x, system.get_box_as_ref().v3x);
        assert_approx_eq!(f32, simbox.v3y, system.get_box_as_ref().v3y);

        simbox.v1x = 10.3;
        simbox.v2x = 8.4;
    }

    #[test]
    fn get_atoms_copy() {
        let system = System::from_file("test_files/example.gro").unwrap();

        let mut atoms = system.get_atoms_copy();

        for (extracted_atom, system_atom) in atoms.iter().zip(system.get_atoms_as_ref().iter()) {
            assert_eq!(
                system_atom.get_atom_number(),
                extracted_atom.get_atom_number()
            );
        }

        let _ = atoms.pop();
        assert_eq!(atoms.len(), 16843);
        assert_eq!(system.get_atoms_as_ref().len(), 16844);
    }

    #[test]
    fn get_groups_copy() {
        let system = System::from_file("test_files/example_box9.gro").unwrap();

        let mut groups = system.get_groups_copy();

        assert!(groups.contains_key("all"));

        let new_group = Group::from_indices(vec![1, 3, 6, 8], 1000);
        groups.insert("Test".to_string(), new_group);

        assert!(groups.contains_key("Test"));
        assert!(!system.get_groups_as_ref().contains_key("Test"));
    }

    #[test]
    fn has_velocities() {
        let system = System::from_file("test_files/example.gro").unwrap();
        assert!(system.has_velocities());

        let system = System::from_file("test_files/example_novelocities.gro").unwrap();
        assert!(!system.has_velocities());
    }

    #[test]
    fn has_forces() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert!(!system.has_forces());

        system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .next();
        assert!(system.has_forces());

        system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .nth(1);
        assert!(!system.has_forces());
    }

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
                assert!(e.to_string().contains("already exists"));
            }
        }

        assert!(system.group_exists("Membrane"));
        assert_eq!(system.group_get_n_atoms("Membrane").unwrap(), 6144);
        for i in 61..=6204 {
            assert!(system.group_isin("Membrane", i).unwrap());
        }
    }

    #[test]
    fn group_create_macro() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.group_create("Membrane", "@membrane").unwrap();

        assert!(system.group_exists("Membrane"));
        assert_eq!(system.group_get_n_atoms("Membrane").unwrap(), 6144);

        for i in 61..=6204 {
            assert!(system.group_isin("Membrane", i).unwrap());
        }

        system.group_create("Protein", "@protein").unwrap();

        assert!(system.group_exists("Protein"));
        assert_eq!(system.group_get_n_atoms("Protein").unwrap(), 61);

        for i in 0..=60 {
            assert!(system.group_isin("Protein", i).unwrap());
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
    fn atoms_extract() {
        let system = System::from_file("test_files/example.gro").unwrap();

        let mut atoms = system.atoms_extract();

        for (extracted_atom, system_atom) in atoms.iter().zip(system.get_atoms_as_ref().iter()) {
            assert_eq!(
                system_atom.get_atom_number(),
                extracted_atom.get_atom_number()
            );
        }

        let _ = atoms.pop();
        assert_eq!(atoms.len(), 16843);
        assert_eq!(system.get_atoms_as_ref().len(), 16844);
    }

    #[test]
    fn group_extract() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let mut ions = system.group_extract("ION").unwrap();

        for (extracted_atom, system_atom) in ions.iter().zip(system.group_iter("ION").unwrap()) {
            assert_eq!(
                system_atom.get_atom_number(),
                extracted_atom.get_atom_number()
            );
        }

        let _ = ions.pop();
        assert_eq!(ions.len(), 239);
        assert_eq!(system.group_get_n_atoms("ION").unwrap(), 240);
    }

    #[test]
    fn group_extract_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        match system.group_extract("Nonexistent") {
            Err(GroupError::NotFound(e)) => assert_eq!(e, "Nonexistent"),
            Ok(_) => panic!("Group extracting should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn get_atom_as_ref() {
        let system = System::from_file("test_files/example.gro").unwrap();

        assert!(system.get_atom_as_ref(0).is_err());
        assert!(system.get_atom_as_ref(16845).is_err());

        let atom = system.get_atom_as_ref(1).unwrap();
        assert_eq!(atom.get_atom_number(), 1);

        let atom = system.get_atom_as_ref(16844).unwrap();
        assert_eq!(atom.get_atom_number(), 16844);
    }

    #[test]
    fn get_atom_as_ref_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        assert!(system.get_atom_as_ref_mut(0).is_err());
        assert!(system.get_atom_as_ref_mut(16845).is_err());

        let atom = system.get_atom_as_ref_mut(1).unwrap();
        assert_eq!(atom.get_atom_number(), 1);

        let atom = system.get_atom_as_ref_mut(16844).unwrap();
        assert_eq!(atom.get_atom_number(), 16844);
    }

    #[test]
    fn get_atom_copy() {
        let system = System::from_file("test_files/example.gro").unwrap();

        assert!(system.get_atom_copy(0).is_err());
        assert!(system.get_atom_copy(16845).is_err());

        let atom = system.get_atom_copy(1).unwrap();
        assert_eq!(atom.get_atom_number(), 1);

        let atom = system.get_atom_copy(16844).unwrap();
        assert_eq!(atom.get_atom_number(), 16844);
    }

    #[test]
    fn group_by_resid() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let (code, residues) = system.group_by_resid();
        if let Err(_) = code {
            panic!("Overlapping groups in group_by_resid.");
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
    fn group_by_resid_warnings() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system
            .read_ndx("test_files/index_group_by_res.ndx")
            .unwrap();

        let (code, residues) = system.group_by_resid();
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
    fn group_by_resid_broken() {
        let mut system = System::from_file("test_files/example_shuffled_residues.gro").unwrap();

        let (code, residues) = system.group_by_resid();
        if let Err(_) = code {
            panic!("Overlapping groups in group_by_resid.");
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
    fn group_by_resname() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let (code, residues) = system.group_by_resname();
        if let Err(_) = code {
            panic!("Overlapping groups in group_by_resname.");
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
    fn group_by_resname_warnings() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system
            .read_ndx("test_files/index_group_by_resname.ndx")
            .unwrap();

        let (code, residues) = system.group_by_resname();
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
}
