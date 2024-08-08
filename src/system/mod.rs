// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of the System structure and its methods.

use hashbrown::HashMap;
use indexmap::IndexMap;
use std::collections::HashSet;
use std::error::Error;
use std::path::Path;

use crate::errors::{AtomError, GroupError, ParseFileError, SimBoxError};
use crate::files::FileType;
use crate::io::{gro_io, pqr_io};
use crate::io::{pdb_io, tpr_io};
use crate::structures::{atom::Atom, group::Group, simbox::SimBox, vector3d::Vector3D};

mod analysis;
mod groups;
pub mod guess;
pub(crate) mod iterating;
mod labeled_atoms;
mod modifying;
#[cfg(any(feature = "parallel", doc))]
mod parallel;
mod utility;

#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(deny_unknown_fields))]
pub struct System {
    /// Name of the molecular system.
    name: String,
    /// Vector of atoms in the system.
    atoms: Vec<Atom>,
    /// Size of the simulation box. (Optional.)
    simulation_box: Option<SimBox>,
    /// Groups of atoms associated with the system.
    groups: IndexMap<String, Group>,
    /// Atoms that have been specifically labeled with a string.
    /// Each atom can have multiple labels, but one label specifies a single atom.
    labeled_atoms: HashMap<String, usize>,
    /// Current simulation step.
    simulation_step: u64,
    /// Current simulation time in picoseconds.
    simulation_time: f32,
    /// Precision of the coordinates.
    coordinates_precision: u64,
    /// Lambda.
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
    /// let system = System::new(name, atoms, Some(simulation_box));
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
    pub fn new(name: &str, atoms: Vec<Atom>, simulation_box: Option<SimBox>) -> Self {
        let mut system = System {
            name: name.to_string(),
            atoms,
            simulation_box,
            groups: IndexMap::new(),
            labeled_atoms: HashMap::new(),
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

    /// Create a new System by reading a gro, pdb, pqr, or tpr file.
    /// The method will attempt to automatically recognize gro, pdb, tpr or a pqr file based on the file extension.
    ///
    /// ## Returns
    /// `System` structure if successful.
    /// `ParseFileError` if the file format is not supported.
    /// `ParseGroError` if parsing of the gro file fails.
    /// `ParsePdbError` if parsing of the pdb file fails.
    /// `ParseTprError` if parsing of the tpr file fails.
    /// `ParsePqrError` if parsing of the pqr file fails.
    ///
    /// ## Example
    /// Reading a gro file.
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
    ///
    /// ## Notes
    /// - The returned System structure will contain two default groups "all" and "All"
    /// consisting of all the atoms in the system.
    /// - When reading a pdb file, no connectivity information (bonds) is read, even if it is provided. You can add
    /// connectivity from a pdb file to your system using [`System::add_bonds_from_pdb`]. See more information
    /// about parsing PDB files in [`pdb_io::read_pdb`](`crate::io::pdb_io::read_pdb`).
    /// - Groups are not read from tpr files.
    #[inline(always)]
    pub fn from_file(filename: impl AsRef<Path>) -> Result<Self, Box<dyn Error + Send + Sync>> {
        let format = FileType::from_name(&filename);
        match format {
            FileType::GRO | FileType::PDB | FileType::TPR | FileType::PQR => {
                Self::from_file_with_format(filename, format)
            }
            _ => Err(Box::from(ParseFileError::UnknownExtension(Box::from(
                filename.as_ref(),
            )))),
        }
    }

    /// Create a new System by reading a file with the specified format.
    /// Same as [`System::from_file`](`crate::system::System::from_file`), but no automatic recognition of the file type is performed.
    ///
    /// ## Example
    /// Reading a file without an extension as a tpr file.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use groan_rs::files::FileType;
    /// #
    /// let system = match System::from_file_with_format("system", FileType::TPR) {
    ///     Ok(x) => x,
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;
    ///     }
    /// };
    /// ```
    pub fn from_file_with_format(
        filename: impl AsRef<Path>,
        filetype: FileType,
    ) -> Result<Self, Box<dyn Error + Send + Sync>> {
        match filetype {
            FileType::GRO => gro_io::read_gro(filename).map_err(Box::from),
            FileType::PDB => pdb_io::read_pdb(filename).map_err(Box::from),
            FileType::TPR => tpr_io::read_tpr(filename).map_err(Box::from),
            FileType::PQR => pqr_io::read_pqr(filename).map_err(Box::from),
            _ => Err(Box::from(ParseFileError::UnsupportedFileType(filetype))),
        }
    }

    /// Create two groups each containing all atoms in the system: "all" and "All".
    ///
    /// ## Returns
    /// - `Ok` if both groups were created or GroupError in case any group with the same name already exists.
    fn group_create_all(&mut self) -> Result<(), GroupError> {
        self.group_create_from_ranges("all", vec![(0, self.get_n_atoms())])?;
        self.group_create_from_ranges("All", vec![(0, self.get_n_atoms())])?;

        self.get_groups_as_mut()
            .get_mut("all")
            .expect("FATAL GROAN ERROR | System::group_create_all | Group 'all' is not available immediately after its construction.")
            .print_ndx = false;

        self.get_groups_as_mut()
            .get_mut("All")
            .expect("FATAL GROAN ERROR | System::group_create_all | Group 'All' is not available immediately after its construction.")
            .print_ndx = false;

        Ok(())
    }

    /// Get the name of the molecular system.
    #[inline(always)]
    pub fn get_name(&self) -> &str {
        &self.name
    }

    /// Get immutable reference to the atoms in the system.
    #[inline(always)]
    pub fn get_atoms_as_ref(&self) -> &Vec<Atom> {
        &self.atoms
    }

    /// Get mutable reference to the atoms in the system.
    ///
    /// ## Warning
    /// - Note that manually changing the `atoms` of the system
    /// can cause the system to become invalid. Other functions may then not work correctly.
    /// - Notably, no atoms can be added or removed from the `atoms` vector as such
    /// operation would make all the groups associated with the system invalid. The same goes
    /// for reordering the atoms.
    /// - The properties of the individual atoms can however be safely changed.
    #[inline(always)]
    pub(crate) fn get_atoms_as_mut(&mut self) -> &mut Vec<Atom> {
        &mut self.atoms
    }

    /// Get copy of the atoms in the system.
    #[inline(always)]
    pub fn get_atoms_copy(&self) -> Vec<Atom> {
        self.atoms.clone()
    }

    /// Get immutable reference to the groups in the system.
    #[inline(always)]
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
    #[inline(always)]
    pub(crate) fn get_groups_as_mut(&mut self) -> &mut IndexMap<String, Group> {
        &mut self.groups
    }

    /// Get copy of the groups in the system.
    #[inline(always)]
    pub fn get_groups_copy(&self) -> IndexMap<String, Group> {
        self.groups.clone()
    }

    /// Get immutable reference to the simulation box.
    #[inline(always)]
    pub fn get_box_as_ref(&self) -> Option<&SimBox> {
        self.simulation_box.as_ref()
    }

    /// Check whether the system has a simulation box.
    #[inline(always)]
    pub fn has_box(&self) -> bool {
        self.simulation_box.is_some()
    }

    /// Get center of the simulation box.
    ///
    /// ## Returns
    /// - `Vector3D` if successful.
    /// - `SimBoxError` if the system has no simulation box
    /// or if the simulation box is not orthogonal.
    #[inline(always)]
    pub fn get_box_center(&self) -> Result<Vector3D, SimBoxError> {
        match &self.simulation_box {
            Some(simbox) if simbox.is_orthogonal() => Ok(Vector3D::new(
                simbox.x / 2.0f32,
                simbox.y / 2.0f32,
                simbox.z / 2.0f32,
            )),
            Some(_) => Err(SimBoxError::NotOrthogonal),
            None => Err(SimBoxError::DoesNotExist),
        }
    }

    /// Get mutable reference to the simulation box.
    #[inline(always)]
    pub fn get_box_as_mut(&mut self) -> Option<&mut SimBox> {
        self.simulation_box.as_mut()
    }

    /// Get copy of the simulation box.
    #[inline(always)]
    pub fn get_box_copy(&self) -> Option<SimBox> {
        self.simulation_box.as_ref().cloned()
    }

    /// Get the number of atoms in the system.
    #[inline(always)]
    pub fn get_n_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Get the number of groups in the system. This counts all groups, even the default ones.
    #[inline(always)]
    pub fn get_n_groups(&self) -> usize {
        self.groups.len()
    }

    /// Get the current simulation time.
    #[inline(always)]
    pub fn get_simulation_time(&self) -> f32 {
        self.simulation_time
    }

    /// Get the current simulation step.
    #[inline(always)]
    pub fn get_simulation_step(&self) -> u64 {
        self.simulation_step
    }

    /// Get the precision of the coordinates.
    #[inline(always)]
    pub fn get_precision(&self) -> u64 {
        self.coordinates_precision
    }

    /// Get the simulation lambda.
    #[inline(always)]
    pub fn get_lambda(&self) -> f32 {
        self.lambda
    }

    /// Set the simulation time.
    #[inline(always)]
    pub fn set_simulation_time(&mut self, time: f32) {
        self.simulation_time = time;
    }

    /// Set the simulation step.
    #[inline(always)]
    pub fn set_simulation_step(&mut self, step: u64) {
        self.simulation_step = step;
    }

    /// Set simulation box.
    #[inline(always)]
    pub fn set_box(&mut self, sim_box: SimBox) {
        self.simulation_box = Some(sim_box);
    }

    /// Set simulation box to `None`.
    #[inline(always)]
    pub fn reset_box(&mut self) {
        self.simulation_box = None;
    }

    /// Set precision of the coordinates.
    #[inline(always)]
    pub fn set_precision(&mut self, precision: u64) {
        self.coordinates_precision = precision;
    }

    /// Set the simulation lambda.
    #[inline(always)]
    pub fn set_lambda(&mut self, lambda: f32) {
        self.lambda = lambda;
    }

    /// Get reference atoms of all polyatomic molecules.
    /// This is mostly for internal use of the `groan_rs` library.
    #[inline(always)]
    pub fn get_mol_references(&self) -> Option<&Vec<usize>> {
        self.mol_references.as_ref()
    }

    /// Reset reference atoms of molecules.
    ///
    /// ## Notes
    /// - **This function must be called every time topology
    /// of the system is changed**.
    /// - (Safe native groan library functions handle this for you.)
    #[inline(always)]
    pub fn reset_mol_references(&mut self) {
        self.mol_references = None;
    }

    /// Set reference atoms of molecules.
    #[inline(always)]
    pub(crate) fn set_mol_references(&mut self, indices: Vec<usize>) {
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
    #[inline(always)]
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
    #[inline(always)]
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
    #[inline(always)]
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
    #[inline(always)]
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
    #[inline(always)]
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
    #[inline(always)]
    pub fn group_extract(&self, name: &str) -> Result<Vec<Atom>, GroupError> {
        Ok(self.group_iter(name)?.cloned().collect())
    }

    /// Get immutable reference to an atom at target index. Atoms are indexed starting from 0.
    ///
    /// ## Returns
    /// Reference to `Atom` structure or `AtomError::OutOfRange` if `index` is out of range.
    #[inline(always)]
    pub fn get_atom_as_ref(&self, index: usize) -> Result<&Atom, AtomError> {
        self.atoms.get(index).ok_or(AtomError::OutOfRange(index))
    }

    /// Get mutable reference to an atom at target index. Atoms are indexed starting from 0.
    ///
    /// ## Returns
    /// Mutable reference to `Atom` structure or `AtomError::OutOfRange` if `index` is out of range.
    #[inline(always)]
    pub fn get_atom_as_mut(&mut self, index: usize) -> Result<&mut Atom, AtomError> {
        self.atoms
            .get_mut(index)
            .ok_or(AtomError::OutOfRange(index))
    }

    /// Get copy of an atom with target index. Atoms are indexed starting from 0.
    ///
    /// ## Returns
    /// Copy of an `Atom` structure or `AtomError::OutOfRange` if `index` is out of range
    #[inline(always)]
    pub fn get_atom_copy(&self, index: usize) -> Result<Atom, AtomError> {
        self.atoms
            .get(index)
            .cloned()
            .ok_or(AtomError::OutOfRange(index))
    }

    /// Get immutable reference to an atom at taget index WITHOUT performing boundary checks.
    /// Atoms are indexed starting from 0.
    ///
    /// ## Safety
    /// `index` must be lower than the number of atoms in the system.
    ///
    /// ## Notes
    /// - Always prefer to use [`System::get_atom_as_ref`], unless you are sure that the
    /// boundary checks measurably slow down your application.
    #[inline(always)]
    pub unsafe fn get_atom_unchecked_as_ref(&self, index: usize) -> &Atom {
        self.atoms.get_unchecked(index)
    }

    /// Get mutable reference to an atom at target index WITHOUT performing boundary checks.
    /// Atoms are indexed starting from 0.
    ///
    /// ## Safety
    /// `index` must be lower than the number of atoms in the system.
    ///
    /// ## Notes
    /// - Always prefer to use [`System::get_atom_as_mut`], unless you are sure that the
    /// boundary checks measurably slow down your application.
    #[inline(always)]
    pub unsafe fn get_atom_unchecked_as_mut(&mut self, index: usize) -> &mut Atom {
        self.atoms.get_unchecked_mut(index)
    }

    /// Get copy of an atom with target index WITHOUT performing boundary checks.
    /// Atoms are indexed starting from 0.
    ///
    /// ## Safety
    /// `index` must be lower than the number of atoms in the system.
    ///
    /// ## Notes
    /// - Always prefer to use [`System::get_atom_copy`], unless you are sure that the
    /// boundary checks measurably slow down your application.
    #[inline(always)]
    pub unsafe fn get_atom_unchecked_copy(&self, index: usize) -> Atom {
        self.atoms.get_unchecked(index).clone()
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use crate::{
        errors::ParsePdbConnectivityError,
        structures::element::Elements,
        test_utilities::utilities::{compare_atoms, compare_atoms_tpr_with_pdb, compare_box},
    };

    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn new() {
        let system = System::new(
            "System generated using the `groan_rs` library.",
            Vec::new(),
            Some([1.5, 3.3, 0.8].into()),
        );

        assert_eq!(
            system.get_name(),
            "System generated using the `groan_rs` library."
        );
        assert_eq!(system.get_atoms_as_ref().len(), 0);

        let simbox = system.get_box_as_ref().unwrap();

        assert_approx_eq!(f32, simbox.v1x, 1.5f32);
        assert_approx_eq!(f32, simbox.v2y, 3.3f32);
        assert_approx_eq!(f32, simbox.v3z, 0.8f32);
        assert_eq!(simbox.v1y, 0.0f32);
        assert_eq!(simbox.v1z, 0.0f32);
        assert_eq!(simbox.v2x, 0.0f32);
        assert_eq!(simbox.v2z, 0.0f32);
        assert_eq!(simbox.v3x, 0.0f32);
        assert_eq!(simbox.v3y, 0.0f32);

        assert!(system.group_exists("all"));
    }

    #[test]
    fn from_file() {
        let system_gro = System::from_file("test_files/example_novelocities.gro").unwrap();

        assert_eq!(system_gro.get_name(), "Buforin II peptide P11L");
        assert_eq!(system_gro.get_n_atoms(), 50);

        let simbox = system_gro.get_box_as_ref().unwrap();
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

        let simbox = system_pdb.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.x, 6.0861);
        assert_approx_eq!(f32, simbox.y, 6.0861);
        assert_approx_eq!(f32, simbox.z, 6.0861);

        assert_eq!(simbox.v1y, 0.0f32);
        assert_eq!(simbox.v1z, 0.0f32);
        assert_eq!(simbox.v2x, 0.0f32);

        assert_eq!(simbox.v2z, 0.0f32);
        assert_eq!(simbox.v3x, 0.0f32);
        assert_eq!(simbox.v3y, 0.0f32);

        let system_pqr = System::from_file("test_files/example.pqr").unwrap();
        assert_eq!(system_pqr.get_name(), "Buforin II peptide P11L");
        assert_eq!(system_pqr.get_n_atoms(), 50);

        let simbox = system_pqr.get_box_as_ref().unwrap();
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

        // compare atoms from PQR and PDB file
        for (pqra, pdba) in system_pqr.atoms_iter().zip(system_pdb.atoms_iter()) {
            assert_eq!(pqra.get_residue_number(), pdba.get_residue_number());
            assert_eq!(pqra.get_residue_name(), pdba.get_residue_name());
            assert_eq!(pqra.get_atom_number(), pdba.get_atom_number());
            assert_eq!(pqra.get_atom_name(), pdba.get_atom_name());
            assert_approx_eq!(
                f32,
                pqra.get_position().unwrap().x,
                pdba.get_position().unwrap().x
            );
            assert_approx_eq!(
                f32,
                pqra.get_position().unwrap().y,
                pdba.get_position().unwrap().y
            );
            assert_approx_eq!(
                f32,
                pqra.get_position().unwrap().z,
                pdba.get_position().unwrap().z
            );

            assert_eq!(pqra.get_velocity(), pdba.get_velocity());
            assert_eq!(pqra.get_force(), pdba.get_force());
            assert_eq!(pqra.get_chain(), pdba.get_chain());
        }
    }

    #[test]
    fn from_file_tpr() {
        let system_tpr = System::from_file("test_files/aa_for_testing_tpr.tpr").unwrap();
        let mut system_pdb = System::from_file("test_files/aa_for_testing_tpr.pdb").unwrap();
        system_pdb
            .add_bonds_from_pdb("test_files/aa_for_testing_tpr.pdb")
            .unwrap();
        system_pdb.guess_elements(Elements::default()).unwrap();

        assert_eq!(system_tpr.get_name(), system_pdb.get_name());
        compare_box(
            system_tpr.get_box_as_ref().unwrap(),
            system_pdb.get_box_as_ref().unwrap(),
        );

        // compare atoms (and bonds)
        for (atom1, atom2) in system_tpr.atoms_iter().zip(system_pdb.atoms_iter()) {
            compare_atoms_tpr_with_pdb(atom1, atom2);
        }
    }

    #[test]
    fn from_file_tpr_triclinic() {
        let system_tpr = System::from_file("test_files/triclinic.tpr").unwrap();
        let system_gro = System::from_file("test_files/triclinic.gro").unwrap();

        compare_box(
            system_tpr.get_box_as_ref().unwrap(),
            system_gro.get_box_as_ref().unwrap(),
        );
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
    fn from_file_with_format() {
        let system_with_format =
            System::from_file_with_format("test_files/example.gro", FileType::GRO).unwrap();
        let system_auto = System::from_file("test_files/example.gro").unwrap();

        assert_eq!(system_with_format.get_n_atoms(), system_auto.get_n_atoms());

        for (a1, a2) in system_with_format
            .atoms_iter()
            .zip(system_auto.atoms_iter())
        {
            crate::test_utilities::utilities::compare_atoms(a1, a2)
        }
    }

    #[test]
    fn from_file_with_format_unsupported() {
        match System::from_file_with_format("test_files/example.gro", FileType::XTC) {
            Ok(_) => panic!("Parsing should have failed."),
            Err(e) => assert!(e.to_string().contains("xtc")),
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
    fn reset_box() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert!(system.simulation_box.is_some());
        system.reset_box();
        assert!(system.simulation_box.is_none());
    }

    #[test]
    fn get_box_copy() {
        let system = System::from_file("test_files/example_box9.gro").unwrap();

        let simbox = system.get_box_copy().unwrap();
        let original_simbox = system.get_box_as_ref().unwrap();

        assert_approx_eq!(f32, simbox.x, original_simbox.x);
        assert_approx_eq!(f32, simbox.y, original_simbox.y);
        assert_approx_eq!(f32, simbox.z, original_simbox.z);

        assert_approx_eq!(f32, simbox.v1y, original_simbox.v1y);
        assert_approx_eq!(f32, simbox.v1z, original_simbox.v1z);
        assert_approx_eq!(f32, simbox.v2x, original_simbox.v2x);

        assert_approx_eq!(f32, simbox.v2z, original_simbox.v2z);
        assert_approx_eq!(f32, simbox.v3x, original_simbox.v3x);
        assert_approx_eq!(f32, simbox.v3y, original_simbox.v3y);
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

        system
            .get_atoms_as_mut()
            .get_mut(10)
            .unwrap()
            .set_atom_number(44);

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
    fn get_atom_unchecked_as_ref() {
        let system = System::from_file("test_files/example.gro").unwrap();

        let indices = [0, 329, 4938, 16843];
        for i in indices {
            let atom_safe = system.get_atom_as_ref(i).unwrap();
            let atom_unsafe = unsafe { system.get_atom_unchecked_as_ref(i) };

            compare_atoms(atom_safe, atom_unsafe);
        }
    }

    #[test]
    fn get_atom_as_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        assert!(system.get_atom_as_mut(16844).is_err());

        let atom = system.get_atom_as_mut(0).unwrap();
        assert_eq!(atom.get_atom_number(), 1);

        let atom = system.get_atom_as_mut(16843).unwrap();
        assert_eq!(atom.get_atom_number(), 16844);
    }

    #[test]
    fn get_atom_unchecked_as_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let indices = [0, 329, 4938, 16843];
        for i in indices {
            let atom_unsafe = unsafe { system.get_atom_unchecked_as_mut(i) as *mut Atom };
            let atom_safe = system.get_atom_as_mut(i).unwrap();
            compare_atoms(atom_safe, unsafe { &*atom_unsafe });
        }
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

    #[test]
    fn get_atom_unchecked_copy() {
        let system = System::from_file("test_files/example.gro").unwrap();

        let indices = [0, 329, 4938, 16843];
        for i in indices {
            let atom_safe = system.get_atom_copy(i).unwrap();
            let atom_unsafe = unsafe { system.get_atom_unchecked_copy(i) };

            compare_atoms(&atom_safe, &atom_unsafe);
        }
    }

    #[test]
    fn get_box_center() {
        let system = System::from_file("test_files/example.gro").unwrap();
        let center = system.get_box_center().unwrap();

        assert_approx_eq!(f32, center.x, 6.506655);
        assert_approx_eq!(f32, center.y, 6.506655);
        assert_approx_eq!(f32, center.z, 5.626735);
    }

    #[test]
    fn get_box_center_nosimbox() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.reset_box();

        match system.get_box_center() {
            Ok(_) => panic!("Function should have failed."),
            Err(SimBoxError::DoesNotExist) => (),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn get_box_center_notorthogonal() {
        let system = System::from_file("test_files/octahedron.gro").unwrap();

        match system.get_box_center() {
            Ok(_) => panic!("Function should have failed."),
            Err(SimBoxError::NotOrthogonal) => (),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }
    }
}

#[cfg(test)]
#[cfg(feature = "serde")]
mod serde_tests {
    use std::fs::read_to_string;

    use super::*;

    #[test]
    fn system_to_yaml() {
        let mut system = System::from_file("test_files/protein.gro").unwrap();
        system.group_create("Sidechains", "name r'^SC.*'").unwrap();

        let string = serde_yaml::to_string(&system).unwrap();
        let expected = read_to_string("test_files/serde_system.yaml").unwrap();

        assert_eq!(string, expected);
    }

    #[test]
    fn system_from_yaml() {
        let string = read_to_string("test_files/serde_system.yaml").unwrap();
        let system: System = serde_yaml::from_str(&string).unwrap();

        assert_eq!(system.get_n_atoms(), 61);
        assert_eq!(system.get_n_groups(), 3);
        assert!(system.get_box_as_ref().is_some());
    }
}
