// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of the `System` structure and methods for constructing the `System` and accessing its properties.

use indexmap::IndexMap;
use std::collections::HashSet;
use std::error::Error;
use std::path::Path;

use crate::errors::{AtomError, GroupError, ParseFileError};
use crate::files::FileType;
use crate::io::gro_io;
use crate::io::pdb_io;
use crate::structures::{atom::Atom, group::Group, simbox::SimBox, vector3d::Vector3D};

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
    /// Current simulation time in picoseconds.
    simulation_time: f32,
    /// Precision of the coordinates.
    coordinates_precision: u64,
    /// Lambda
    lambda: f32,
    /// Reference atoms for all polyatomic molecules.
    /// (Index of the first atom of each polyatomic molecule.)
    /// All functions changing the topology of the system, must set
    /// `mol_references` to `None`.
    mol_references: Option<Vec<usize>>,
}

/// ## Methods for creating `System` structures and accessing their properties.
impl System {
    /// Create new System structure with a given name from the provided vector of atoms and simulation box.
    ///
    /// ## Notes
    /// - The returned `System` structure will contain two default groups "all" and "All",
    /// each consisting of all the atoms in the system.
    ///
    /// ## Example 1: Manually creating a system
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
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
    /// # use groan_rs::prelude::*;
    /// #
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
    /// # use groan_rs::prelude::*;
    /// #
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
            mol_references: None,
        };

        match system.group_create_all() {
            Err(_) => {
                panic!("FATAL GROAN ERROR | System::new | Group 'all' or 'All' already exists as the System is created.");
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
    /// # use groan_rs::prelude::*;
    /// #
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
    pub fn from_file(filename: impl AsRef<Path>) -> Result<Self, Box<dyn Error + Send + Sync>> {
        match FileType::from_name(&filename) {
            FileType::GRO => gro_io::read_gro(filename).map_err(Box::from),
            FileType::PDB => pdb_io::read_pdb(filename).map_err(Box::from),
            _ => Err(Box::from(ParseFileError::UnknownExtension(Box::from(
                filename.as_ref(),
            )))),
        }
    }

    /// Create two groups each containing all atoms in the system: "all" and "All".
    ///
    /// ## Returns
    /// - `Ok` if both groups were created or GroupError in case any group with the same name already exists.
    fn group_create_all(&mut self) -> Result<(), GroupError> {
        self.group_create_from_ranges("all", vec![(0, self.get_n_atoms())])?;
        self.group_create_from_ranges("All", vec![(0, self.get_n_atoms())])?;

        unsafe {
            self.get_groups_as_ref_mut()
                .get_mut("all")
                .expect("FATAL GROAN ERROR | System::group_create_all | Group 'all' is not available immediately after its construction.")
                .print_ndx = false;

            self.get_groups_as_ref_mut()
                .get_mut("All")
                .expect("FATAL GROAN ERROR | System::group_create_all | Group 'All' is not available immediately after its construction.")
                .print_ndx = false;
        }

        Ok(())
    }

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
    /// - Notably, it is forbidden to modify the default groups 'all' and 'All' as changing
    /// these groups may cause the behavior of many other functions associated with `System`
    /// to become incorrect.
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

    /// Get reference atoms of all polyatomic molecules.
    /// This is mostly for internal use of the `groan_rs` library.
    pub fn get_mol_references(&self) -> Option<&Vec<usize>> {
        self.mol_references.as_ref()
    }

    /// Reset reference atoms of molecules.
    ///
    /// ## Notes
    /// - **This function must be called every time topology
    /// of the system is changed**.
    /// - (Safe native groan library functions handle this for you.)
    pub fn reset_mol_references(&mut self) {
        self.mol_references = None;
    }

    /// Set reference atoms of molecules.
    ///
    /// ## Safety
    /// Modifying `mol_references` may break the system.
    /// You should not set `mol_references` manually unless you know what you are doing.
    /// Which you do not. There really is no reason to use this method.
    pub unsafe fn set_mol_references(&mut self, indices: Vec<usize>) {
        self.mol_references = Some(indices);
    }

    /// Check whether positions are present.
    ///
    /// ## Returns
    /// `true` if all of the atoms in the system have information about their positions.
    /// `false` otherwise.
    ///
    /// ## Notes
    /// - Complexity of this operation is O(n), where n is the number of atoms in the system.
    pub fn has_positions(&self) -> bool {
        self.atoms.iter().all(|atom| atom.has_position())
    }

    /// Check whether velocities are present.
    ///
    /// ## Returns
    /// `true` if all of the atoms in the system have information about their velocities.
    /// `false` otherwise.
    ///
    /// ## Notes
    /// - Complexity of this operation is O(n), where n is the number of atoms in the system.
    pub fn has_velocities(&self) -> bool {
        self.atoms.iter().all(|atom| atom.has_velocity())
    }

    /// Check whether forces are present.
    ///
    /// ## Returns
    /// `true` if all of the atoms in the system have information about force acting on them.
    /// `false` otherwise.
    ///
    /// ## Notes
    /// - Complexity of this operation is O(n), where n is the number of atoms in the system.
    pub fn has_forces(&self) -> bool {
        self.atoms.iter().all(|atom| atom.has_force())
    }

    /// Check whether there are any atoms in the system which share atom number.
    ///
    /// ## Returns
    /// `true` if at least two atoms share the atom number. `false` otherwise.
    ///
    /// ## Notes
    /// - Complexity of this operation is O(n), where n is the number of atoms in the system.
    pub fn has_duplicate_atom_numbers(&self) -> bool {
        let mut set = HashSet::new();

        for atom in self.atoms.iter() {
            if !set.insert(atom.get_atom_number()) {
                return true;
            }
        }

        false
    }

    /// Check whether connectivity information is available for the system.
    ///
    /// ## Returns
    /// `true` if at least one atom in the system has more than 0 bonds.
    /// `false` otherwise.
    ///
    /// ## Notes
    /// - Complexity of this operation is O(n), where n is the number of atoms in the system.
    pub fn has_bonds(&self) -> bool {
        self.atoms.iter().any(|atom| atom.get_n_bonded() > 0)
    }

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

    /// Get immutable reference to an atom at target index. Atoms are indexed starting from 0.
    ///
    /// ## Returns
    /// Reference to `Atom` structure or `AtomError::OutOfRange` if `index` is out of range.
    pub fn get_atom_as_ref(&self, index: usize) -> Result<&Atom, AtomError> {
        if index >= self.atoms.len() {
            return Err(AtomError::OutOfRange(index));
        }

        Ok(&self.atoms[index])
    }

    /// Get mutable reference to an atom with target index. Atoms are indexed starting from 0.
    ///
    /// ## Returns
    /// Mutable reference to `Atom` structure or `AtomError::OutOfRange` if `index` is out of range.
    pub fn get_atom_as_ref_mut(&mut self, index: usize) -> Result<&mut Atom, AtomError> {
        if index >= self.atoms.len() {
            return Err(AtomError::OutOfRange(index));
        }

        Ok(&mut self.atoms[index])
    }

    /// Get copy of an atom with target index. Atoms are indexed starting from 0.
    ///
    /// ## Returns
    /// Copy of an `Atom` structure or `AtomError::OutOfRange` if `index` is out of range
    pub fn get_atom_copy(&self, index: usize) -> Result<Atom, AtomError> {
        if index >= self.atoms.len() {
            return Err(AtomError::OutOfRange(index));
        }

        Ok(self.atoms[index].clone())
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use crate::errors::ParsePdbConnectivityError;

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
            assert_approx_eq!(
                f32,
                groa.get_position().unwrap().x,
                pdba.get_position().unwrap().x
            );
            assert_approx_eq!(
                f32,
                groa.get_position().unwrap().y,
                pdba.get_position().unwrap().y
            );
            assert_approx_eq!(
                f32,
                groa.get_position().unwrap().z,
                pdba.get_position().unwrap().z
            );

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
    fn has_positions() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert!(system.has_positions());

        system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .nth(1);
        assert!(!system.has_positions());

        system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .nth(2);
        assert!(!system.has_positions());

        system
            .trr_iter("test_files/short_trajectory.trr")
            .unwrap()
            .nth(3);
        assert!(system.has_positions());
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
    fn has_duplicate_atom_numbers() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert!(!system.has_duplicate_atom_numbers());

        unsafe {
            system
                .get_atoms_as_ref_mut()
                .get_mut(10)
                .unwrap()
                .set_atom_number(44);
        }

        assert!(system.has_duplicate_atom_numbers());
    }

    #[test]
    fn has_bonds_1() {
        let mut system = System::from_file("test_files/example.pdb").unwrap();
        assert!(!system.has_bonds());

        match system.add_bonds_from_pdb("test_files/example.pdb") {
            Ok(_) => panic!("Should have returned NoBonds warning."),
            Err(ParsePdbConnectivityError::NoBondsWarning(_)) => assert!(!system.has_bonds()),
            Err(e) => panic!("Function failed with error type `{:?}`.", e),
        }
    }

    #[test]
    fn has_bonds_2() {
        let mut system = System::from_file("test_files/example.pdb").unwrap();
        assert!(!system.has_bonds());

        system
            .add_bonds_from_pdb("test_files/bonds_for_example.pdb")
            .unwrap();
        assert!(system.has_bonds());
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

        assert!(system.get_atom_as_ref(16844).is_err());

        let atom = system.get_atom_as_ref(0).unwrap();
        assert_eq!(atom.get_atom_number(), 1);

        let atom = system.get_atom_as_ref(16843).unwrap();
        assert_eq!(atom.get_atom_number(), 16844);
    }

    #[test]
    fn get_atom_as_ref_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        assert!(system.get_atom_as_ref_mut(16844).is_err());

        let atom = system.get_atom_as_ref_mut(0).unwrap();
        assert_eq!(atom.get_atom_number(), 1);

        let atom = system.get_atom_as_ref_mut(16843).unwrap();
        assert_eq!(atom.get_atom_number(), 16844);
    }

    #[test]
    fn get_atom_copy() {
        let system = System::from_file("test_files/example.gro").unwrap();

        assert!(system.get_atom_copy(16844).is_err());

        let atom = system.get_atom_copy(0).unwrap();
        assert_eq!(atom.get_atom_number(), 1);

        let atom = system.get_atom_copy(16843).unwrap();
        assert_eq!(atom.get_atom_number(), 16844);
    }
}
