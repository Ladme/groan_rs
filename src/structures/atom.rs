// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of the Atom structure and its methods.

use std::io::Write;

use crate::errors::{AtomError, PositionError, WriteGroError, WritePdbError, WritePqrError};
use crate::io::pqr_io::PqrPrecision;
use crate::structures::{
    container::AtomContainer, dimension::Dimension, simbox::SimBox, vector3d::Vector3D,
};

#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(deny_unknown_fields))]
pub struct Atom {
    /// Number (index) of the residue this atom is part of.
    residue_number: usize,
    /// Name of the residue this atom is part of.
    residue_name: String,
    /// Number (index) of the atom.
    atom_number: usize,
    /// Name of the atom.
    atom_name: String,
    /// Identifier of the chain this atom is part of. (Optional)
    chain: Option<char>,
    /// Charge of the atom. (Optional)
    charge: Option<f32>,
    /// Mass of the atom. (Optional)
    mass: Option<f32>,
    /// Van der Waals radius of the atom. (Optional)
    vdw: Option<f32>,
    /// Expected maximal number of bonded atoms. (Optional.)
    expected_max_bonds: Option<u8>,
    /// Expected minimal number of bonded atoms. (Optional.)
    expected_min_bonds: Option<u8>,
    /// Name of the element of the atom. (Optional.)
    element_name: Option<String>,
    /// Symbol of the element of the atom. (Optional.)
    element_symbol: Option<String>,
    /// Position of the atom in 3D space. (Optional)
    position: Option<Vector3D>,
    /// Velocity of the atom. (Optional)
    velocity: Option<Vector3D>,
    /// Force acting on the atom. (Optional)
    force: Option<Vector3D>,
    /// Indices of atoms that are bonded to this atom.
    bonded: AtomContainer,
}

impl Atom {
    /// Create new Atom structure with the specified properties.
    ///
    /// ## Notes
    /// - By default, atom is constructed with `position`, `velocity`, `force`, `chain`, `charge`,
    /// `mass`, `vdw`, `expected_max_bonds`, `expected_min_bonds`, `element_name`, and `element_symbol`
    /// set to `None`. You can provide this information using the `Atom::with_*` methods.
    pub fn new(
        residue_number: usize,
        residue_name: &str,
        atom_number: usize,
        atom_name: &str,
    ) -> Self {
        Atom {
            residue_number,
            residue_name: residue_name.to_string(),
            atom_number,
            atom_name: atom_name.to_string(),
            chain: None,
            charge: None,
            mass: None,
            vdw: None,
            expected_max_bonds: None,
            expected_min_bonds: None,
            element_name: None,
            element_symbol: None,
            position: None,
            velocity: None,
            force: None,
            bonded: AtomContainer::empty(),
        }
    }

    /// Add position information to target atom.
    pub fn with_position(mut self, position: Vector3D) -> Self {
        self.set_position(position);
        self
    }

    /// Add velocity information to target atom.
    pub fn with_velocity(mut self, velocity: Vector3D) -> Self {
        self.set_velocity(velocity);
        self
    }

    /// Add force information to target atom
    pub fn with_force(mut self, force: Vector3D) -> Self {
        self.set_force(force);
        self
    }

    /// Add chain information to target atom.
    pub fn with_chain(mut self, chain: char) -> Self {
        self.set_chain(chain);
        self
    }

    /// Add charge to target atom.
    pub fn with_charge(mut self, charge: f32) -> Self {
        self.set_charge(charge);
        self
    }

    /// Add mass to target atom.
    pub fn with_mass(mut self, mass: f32) -> Self {
        self.set_mass(mass);
        self
    }

    /// Add vdw to target atom.
    pub fn with_vdw(mut self, vdw: f32) -> Self {
        self.set_vdw(vdw);
        self
    }

    /// Add expected maximal number of bonds to target atom.
    pub fn with_expected_max_bonds(mut self, expected_max_bonds: u8) -> Self {
        self.set_expected_max_bonds(expected_max_bonds);
        self
    }

    /// Add expected minimal number of bonds to target atom.
    pub fn with_expected_min_bonds(mut self, expected_min_bonds: u8) -> Self {
        self.set_expected_min_bonds(expected_min_bonds);
        self
    }

    /// Add element name to target atom.
    pub fn with_element_name(mut self, name: &str) -> Self {
        self.set_element_name(name);
        self
    }

    /// Add element name to target atom.
    pub fn with_element_symbol(mut self, symbol: &str) -> Self {
        self.set_element_symbol(symbol);
        self
    }

    /// Get the number of the residue to which the atom belongs.
    pub fn get_residue_number(&self) -> usize {
        self.residue_number
    }

    /// Set the number of the residue to which the atom belongs.
    pub fn set_residue_number(&mut self, resnum: usize) {
        self.residue_number = resnum;
    }

    /// Get the name of the residue to which the atom belongs.
    pub fn get_residue_name(&self) -> &str {
        &self.residue_name
    }

    /// Set the name of the residue to which the atom belongs.
    pub fn set_residue_name(&mut self, resname: &str) {
        self.residue_name = resname.to_string();
    }

    /// Get the number of the atom as presented in gro or pdb file.
    pub fn get_atom_number(&self) -> usize {
        self.atom_number
    }

    /// Set the number of the atom as presented in gro or pdb file.
    pub fn set_atom_number(&mut self, atomnum: usize) {
        self.atom_number = atomnum;
    }

    /// Get the name of the atom.
    pub fn get_atom_name(&self) -> &str {
        &self.atom_name
    }

    /// Set the name of the atom.
    pub fn set_atom_name(&mut self, atomname: &str) {
        self.atom_name = atomname.to_string();
    }

    /// Get the chain this atom is part of.
    pub fn get_chain(&self) -> Option<char> {
        self.chain
    }

    /// Set the chain of the atom.
    pub fn set_chain(&mut self, chain: char) {
        self.chain = Some(chain);
    }

    /// Set the chain of the atom to `None`.
    pub fn reset_chain(&mut self) {
        self.chain = None;
    }

    /// Get the charge of the atom.
    pub fn get_charge(&self) -> Option<f32> {
        self.charge
    }

    /// Set the charge of the atom.
    pub fn set_charge(&mut self, charge: f32) {
        self.charge = Some(charge);
    }

    /// Set the charge of the atom to `None`.
    pub fn reset_charge(&mut self) {
        self.charge = None;
    }

    /// Get the mass of the atom.
    pub fn get_mass(&self) -> Option<f32> {
        self.mass
    }

    /// Set the mass of the atom.
    pub fn set_mass(&mut self, mass: f32) {
        self.mass = Some(mass);
    }

    /// Set mass of the atom to `None`.
    pub fn reset_mass(&mut self) {
        self.mass = None;
    }

    /// Get the vdW radius of the atom.
    pub fn get_vdw(&self) -> Option<f32> {
        self.vdw
    }

    /// Set the vdW radius of the atom.
    pub fn set_vdw(&mut self, vdw: f32) {
        self.vdw = Some(vdw);
    }

    /// Set the vdW radius of the atom to `None`.
    pub fn reset_vdw(&mut self) {
        self.vdw = None;
    }

    /// Get the expected maximal number of bonds of the atom.
    pub fn get_expected_max_bonds(&self) -> Option<u8> {
        self.expected_max_bonds
    }

    /// Set the expected maximal number of bonds of the atom.
    pub fn set_expected_max_bonds(&mut self, expected_max_bonds: u8) {
        self.expected_max_bonds = Some(expected_max_bonds);
    }

    /// Set the expected maximal number of bonds of the atom to `None`.
    pub fn reset_expected_max_bonds(&mut self) {
        self.expected_max_bonds = None;
    }

    /// Get the expected minimal number of bonds of the atom.
    pub fn get_expected_min_bonds(&self) -> Option<u8> {
        self.expected_min_bonds
    }

    /// Set the expected minimal number of bonds of the atom.
    pub fn set_expected_min_bonds(&mut self, expected_min_bonds: u8) {
        self.expected_min_bonds = Some(expected_min_bonds);
    }

    /// Set the expected minimal number of bonds of the atom to `None`.
    pub fn reset_expected_min_bonds(&mut self) {
        self.expected_min_bonds = None;
    }

    /// Get the element name of the atom.
    pub fn get_element_name(&self) -> Option<&str> {
        self.element_name.as_deref()
    }

    /// Set the element name of the atom.
    pub fn set_element_name(&mut self, name: &str) {
        self.element_name = Some(name.to_string());
    }

    /// Set the element name of the atom to `None`.
    pub fn reset_element_name(&mut self) {
        self.element_name = None;
    }

    /// Get the element symbol of the atom.
    pub fn get_element_symbol(&self) -> Option<&str> {
        self.element_symbol.as_deref()
    }

    /// Set the element symbol of the atom.
    pub fn set_element_symbol(&mut self, symbol: &str) {
        self.element_symbol = Some(symbol.to_string());
    }

    /// Set the element symbol of the atom to `None`.
    pub fn reset_element_symbol(&mut self) {
        self.element_symbol = None;
    }

    /// Get the coordinates of the atom.
    pub fn get_position(&self) -> Option<&Vector3D> {
        self.position.as_ref()
    }

    /// Get mutable reference to the coordinates of the atom.
    pub fn get_position_mut(&mut self) -> Option<&mut Vector3D> {
        self.position.as_mut()
    }

    /// Set the coordinates of the atom.
    pub fn set_position(&mut self, pos: Vector3D) {
        self.position = Some(pos);
    }

    /// Set the coordinates of the atom to `None`.
    pub fn reset_position(&mut self) {
        self.position = None;
    }

    /// Set the x-coordinate of the atom.
    pub fn set_position_x(&mut self, x: f32) {
        match self.position {
            None => self.position = Some(Vector3D::new(x, 0.0, 0.0)),
            Some(ref mut pos) => pos.x = x,
        }
    }

    /// Set the y-coordinate of the atom.
    pub fn set_position_y(&mut self, y: f32) {
        match self.position {
            None => self.position = Some(Vector3D::new(0.0, y, 0.0)),
            Some(ref mut pos) => pos.y = y,
        }
    }

    /// Set the z-coordinate of the atom.
    pub fn set_position_z(&mut self, z: f32) {
        match self.position {
            None => self.position = Some(Vector3D::new(0.0, 0.0, z)),
            Some(ref mut pos) => pos.z = z,
        }
    }

    /// Get the velocity vector of the atom.
    pub fn get_velocity(&self) -> Option<&Vector3D> {
        self.velocity.as_ref()
    }

    /// Get mutable reference to the velocity vector of the atom.
    pub fn get_velocity_mut(&mut self) -> Option<&mut Vector3D> {
        self.velocity.as_mut()
    }

    /// Set the velocity vector of the atom.
    pub fn set_velocity(&mut self, vel: Vector3D) {
        self.velocity = Some(vel);
    }

    /// Set the velocity vector of the atom to `None`.
    pub fn reset_velocity(&mut self) {
        self.velocity = None;
    }

    /// Get the vector of the total force acting on the atom.
    pub fn get_force(&self) -> Option<&Vector3D> {
        self.force.as_ref()
    }

    /// Get mutable reference to the total force acting on the atom.
    pub fn get_force_mut(&mut self) -> Option<&mut Vector3D> {
        self.force.as_mut()
    }

    /// Set the vector of the total force acting on the atom.
    pub fn set_force(&mut self, force: Vector3D) {
        self.force = Some(force);
    }

    /// Set the vector of the total force acting on the atom to `None`.
    pub fn reset_force(&mut self) {
        self.force = None;
    }

    /// Get the atoms bonded to this atom.
    pub fn get_bonded(&self) -> &AtomContainer {
        &self.bonded
    }

    /// Set the atoms bonded to this atom.
    ///
    /// ## Safety
    /// This method is only safe to use if all the following requirements are fulfilled:
    ///
    /// a) all indices are valid atom indices of the system,
    ///
    /// b) no index of `indices` corresponds to the index of the this atom,
    ///
    /// c) index of this atom is also added between bonded atoms of all atoms included in `indices`,
    ///
    /// d) `System::reset_mol_references` method is applied to the `System` this atom is part of.
    pub unsafe fn set_bonded(&mut self, indices: Vec<usize>) {
        self.bonded = AtomContainer::from_indices(indices, usize::MAX);
    }

    /// Add index of atom bonded to this atom.
    ///
    /// ## Safety
    /// This method is only safe to use if all the following requirements are fulfilled:
    ///
    /// a) the index is a valid atom index of the system,
    ///
    /// b) the index does not match the index of this atom,
    ///
    /// c) index of this atom is also added between bonded atoms of the atom `index`,
    ///
    /// d) `System::reset_mol_references` method is applied to the `System` this atom is part of.
    pub unsafe fn add_bonded(&mut self, index: usize) {
        self.bonded.add(index, usize::MAX);
    }

    /// Get the number of bonded atoms associated with this atom.
    ///
    /// ## Returns
    /// The number of atoms bonded to this atom. 0 if no connectivity information is available.
    pub fn get_n_bonded(&self) -> usize {
        self.bonded.get_n_atoms()
    }

    /// Check whether the atom has position information.
    ///
    /// ## Returns
    /// `true` if the atom has information about position. `false` otherwise.
    pub fn has_position(&self) -> bool {
        self.position.is_some()
    }

    /// Check whether the atom has velocity information.
    ///
    /// ## Returns
    /// `true` if the atom has information about velocity. `false` otherwise.
    pub fn has_velocity(&self) -> bool {
        self.velocity.is_some()
    }

    /// Check whether the atom has force information.
    ///
    /// ## Returns
    /// `true` if the atom has information about force acting on it. `false` otherwise.
    pub fn has_force(&self) -> bool {
        self.force.is_some()
    }

    /// Translates the position of the atom by the provided Vector3D.
    /// Wraps the atom to the simulation box.
    ///
    /// ## Returns
    /// `Ok` or `AtomError::InvalidPosition` if the atom has an undefined position.
    pub fn translate(&mut self, translate: &Vector3D, sbox: &SimBox) -> Result<(), AtomError> {
        if let Some(ref mut pos) = self.position {
            pos.x += translate.x;
            pos.y += translate.y;
            pos.z += translate.z;

            pos.wrap(sbox);
            Ok(())
        } else {
            Err(AtomError::InvalidPosition(PositionError::NoPosition(
                self.get_atom_number(),
            )))
        }
    }

    /// Translates the position of the atom by the provided Vector3D.
    /// Does **not** wrap the atom to the simulation box.
    ///
    /// ## Returns
    /// `Ok` of `AtomError::InvalidPosition` if the atom has an undefined position.
    pub fn translate_nopbc(&mut self, translate: &Vector3D) -> Result<(), AtomError> {
        if let Some(ref mut pos) = self.position {
            pos.0 += translate.0;
            Ok(())
        } else {
            Err(AtomError::InvalidPosition(PositionError::NoPosition(
                self.get_atom_number(),
            )))
        }
    }

    /// Wrap the atom into the simulation box.
    ///
    /// ## Returns
    /// `Ok` of `AtomError::InvalidPosition` if the atom has an undefined position.
    pub fn wrap(&mut self, sbox: &SimBox) -> Result<(), AtomError> {
        match self.position {
            None => Err(AtomError::InvalidPosition(PositionError::NoPosition(
                self.get_atom_number(),
            ))),
            Some(ref mut pos) => {
                pos.wrap(sbox);
                Ok(())
            }
        }
    }

    /// Write information about the atom in gro format.
    /// Only writes velocities if requested (if `write_velocities == true`).
    ///
    /// ## Notes
    /// - Allows for 0 to 5-letter atom names, 0 to 5-letter residue names, 1 to 5-digit atom numbers, and 1 to 5-digit residue numbers.
    /// - Longer names are shortened, longer numbers are wrapped to 0.
    /// - If atom has no position (or velocity, if requested), 0 is printed out for all dimensions.
    pub fn write_gro(
        &self,
        stream: &mut impl Write,
        write_velocities: bool,
    ) -> Result<(), WriteGroError> {
        let binding = Vector3D::default();
        let position = self.get_position().unwrap_or(&binding);

        let format_atomname = match self.get_atom_name().len() {
            0..=5 => format!("{:>5}", self.get_atom_name()),
            _ => format!(
                "{:>5}",
                self.get_atom_name().chars().take(5).collect::<String>()
            ),
        };

        let format_resname = match self.get_residue_name().len() {
            0..=5 => format!("{:<5}", self.get_residue_name()),
            _ => format!(
                "{:<5}",
                self.get_residue_name().chars().take(5).collect::<String>()
            ),
        };

        if write_velocities {
            let velocity = self.get_velocity().unwrap_or(&binding);

            writeln!(
                stream,
                "{:>5}{}{}{:>5}{:>8.3}{:>8.3}{:>8.3}{:>8.4}{:>8.4}{:>8.4}",
                self.get_residue_number() % 100_000,
                format_resname,
                format_atomname,
                self.get_atom_number() % 100_000,
                position.x,
                position.y,
                position.z,
                velocity.x,
                velocity.y,
                velocity.z
            )
            .map_err(|_| WriteGroError::CouldNotWrite)?;
        } else {
            writeln!(
                stream,
                "{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}",
                self.get_residue_number() % 100_000,
                format_resname,
                format_atomname,
                self.get_atom_number() % 100_000,
                position.x,
                position.y,
                position.z
            )
            .map_err(|_| WriteGroError::CouldNotWrite)?;
        }

        Ok(())
    }

    /// Write information about the atom in pdb format.
    ///
    /// ## Notes
    /// - All atoms are treated as 'ATOM'. 'HETATM' is not used at all.
    /// - Allows for 0 to 4-letter atom names, 0 to 4-letter residue names, 1 to 5-digit atom numbers and 1 to 4-digit residue numbers.
    /// - Longer names are shortened, longer numbers are wrapped to 0.
    /// - If atom has no position, 0 is printed out for all dimensions.
    pub fn write_pdb(&self, stream: &mut impl Write) -> Result<(), WritePdbError> {
        let binding = Vector3D::default();
        let position = self.get_position().unwrap_or(&binding);

        let format_resname = match self.get_residue_name().len() {
            0..=3 => format!("{:>3} ", self.get_residue_name()),
            4 => format!("{:>4}", self.get_residue_name()),
            _ => format!(
                "{:>4}",
                self.get_residue_name().chars().take(4).collect::<String>()
            ),
        };

        let format_atomname = match self.get_atom_name().len() {
            0..=3 => format!(" {:<3}", self.get_atom_name()),
            4 => format!("{:<4}", self.get_atom_name()),
            _ => format!(
                "{:<4}",
                self.get_atom_name().chars().take(4).collect::<String>()
            ),
        };

        let format_chain = self.get_chain().unwrap_or(' ');

        writeln!(
            stream,
            "ATOM  {:>5} {} {}{}{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00            ",
            self.get_atom_number() % 100000,
            format_atomname,
            format_resname,
            format_chain,
            self.get_residue_number() % 10000,
            position.x * 10.0,
            position.y * 10.0,
            position.z * 10.0,
        )
        .map_err(|_| WritePdbError::CouldNotWrite)?;

        Ok(())
    }

    /// Write information about the atom in whitespace-delimited PQR format.
    ///
    /// ## Parameters
    /// - `precision` parameters specify the number of decimal places to be printed for
    /// position, charge and radius.
    ///
    /// ## Notes
    /// - All atoms are treated as 'ATOM'. 'HETATM' is not used at all.
    /// - Allows for any atom and residue names of any length. No numbers are wrapped.
    /// - If atom has no position, 0 is printed out for all dimensions.
    /// - If atom has no charge or no vdw, 0 is used.
    pub fn write_pqr(
        &self,
        stream: &mut impl Write,
        precision: &PqrPrecision,
    ) -> Result<(), WritePqrError> {
        let format_resname = match self.get_residue_name().len() {
            0..=3 => format!("{:>3} ", self.get_residue_name()),
            _ => format!("{} ", self.get_residue_name()),
        };

        let format_atomname = match self.get_atom_name().len() {
            0..=3 => format!(" {:<3}", self.get_atom_name()),
            _ => self.get_atom_name().to_string(),
        };

        let resid = self.get_residue_number();
        let format_resid = match resid {
            0..=999 => format!("{:>4}    ", resid),
            1000..=9999 => format!("{:>5}   ", resid),
            10000..=99999 => format!("{:>6}  ", resid),
            100000..=999999 => format!("{:>7} ", resid),
            1000000..=9999999 => format!("{:>8}", resid),
            _ => format!(" {}", resid),
        };

        let format_atomid = match self.get_atom_number() {
            0..=99999 => format!(" {:>5}", self.get_atom_number()),
            _ => format!("{}", self.get_atom_number()),
        };

        let chain = self.get_chain().unwrap_or(' ');

        let binding = Vector3D::default();
        let position = self.get_position().unwrap_or(&binding);
        let charge = self.get_charge().unwrap_or(0.0);
        let vdw = self.get_vdw().unwrap_or(0.0);

        writeln!(
            stream,
            "ATOM {3} {4} {5}{6}{7} {8:>7.0$} {9:>7.0$} {10:>7.0$} {11:>7.1$} {12:>6.2$}",
            precision.position,
            precision.charge,
            precision.vdw,
            format_atomid,
            format_atomname,
            format_resname,
            chain,
            format_resid,
            position.x * 10.0,
            position.y * 10.0,
            position.z * 10.0,
            charge,
            vdw * 10.0,
        )
        .map_err(|_| WritePqrError::CouldNotWrite)?;

        Ok(())
    }

    /// Calculate distance between two atoms.
    /// Takes periodic boundary conditions into consideration.
    /// Returns oriented distance for 1D problems.
    ///
    /// ## Returns
    /// - `f32` if successful.
    /// - `AtomError::InvalidPosition` if any of the atoms has undefined position.
    ///
    /// ## Warning
    /// - Currently only works with orthogonal simulation boxes.
    ///
    /// ## Example
    /// Calculate distance between two atoms in the xy-plane.
    /// ```
    /// # use groan_rs::prelude::*;
    /// # use float_cmp::assert_approx_eq;
    /// #
    /// let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([1.0, 2.0, 3.0].into());
    /// let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([3.5, 1.0, 2.0].into());
    ///
    /// let simbox = SimBox::from([4.0, 4.0, 4.0]);
    ///
    /// let distance = atom1.distance(&atom2, Dimension::XY, &simbox).unwrap();
    /// assert_approx_eq!(f32, distance, 1.802776);
    /// ```
    ///
    /// ## Notes
    /// - If dimension is `Dimension::None` returns 0.
    pub fn distance(&self, atom: &Atom, dim: Dimension, sbox: &SimBox) -> Result<f32, AtomError> {
        match (&self.position, &atom.position) {
            (None, Some(_) | None) => Err(AtomError::InvalidPosition(PositionError::NoPosition(
                self.get_atom_number(),
            ))),
            (Some(_), None) => Err(AtomError::InvalidPosition(PositionError::NoPosition(
                atom.get_atom_number(),
            ))),
            (Some(ref pos1), Some(ref pos2)) => Ok(pos1.distance(pos2, dim, sbox)),
        }
    }

    /// Calculate distance between two atoms.
    /// **Ignores PBC.**
    /// Returns oriented distance for 1D problems.
    ///
    /// ## Returns
    /// - `f32` if successful.
    /// - `AtomError::InvalidPosition` if any of the atoms has undefined position.
    ///
    /// ## Example
    /// Calculate naive distance between two atoms in the xy-plane.
    /// ```
    /// # use groan_rs::prelude::*;
    /// # use float_cmp::assert_approx_eq;
    /// #
    /// let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([1.0, 2.0, 3.0].into());
    /// let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([3.5, 1.0, 2.0].into());
    ///
    /// let distance = atom1.distance_naive(&atom2, Dimension::XY).unwrap();
    /// assert_approx_eq!(f32, distance, 2.692582);
    /// ```
    ///
    /// ## Notes
    /// - If dimension is `Dimension::None` returns 0.
    pub fn distance_naive(&self, atom: &Atom, dim: Dimension) -> Result<f32, AtomError> {
        match (&self.position, &atom.position) {
            (None, Some(_) | None) => Err(AtomError::InvalidPosition(PositionError::NoPosition(
                self.get_atom_number(),
            ))),
            (Some(_), None) => Err(AtomError::InvalidPosition(PositionError::NoPosition(
                atom.get_atom_number(),
            ))),
            (Some(ref pos1), Some(ref pos2)) => Ok(pos1.distance_naive(pos2, dim)),
        }
    }

    /// Calculate distance between an atom and a point in space.
    /// Takes periodic boundary conditions into consideration.
    /// Returns oriented distance for 1D problems.
    ///
    /// ## Returns
    /// - `f32` if successful.
    /// - `AtomError::InvalidPosition` if the atom has undefined position.
    ///
    /// ## Warning
    /// - Currently only works with orthogonal simulation boxes.
    ///
    /// ## Panics
    /// Panics if the atom has no position.
    ///
    /// ## Example
    /// Calculate distance between an atom and a point in the xy-plane.
    /// ```
    /// # use groan_rs::prelude::*;
    /// # use float_cmp::assert_approx_eq;
    /// #
    /// let atom = Atom::new(1, "LYS", 1, "BB").with_position([1.0, 2.0, 3.0].into());
    /// let point = Vector3D::new(3.5, 1.0, 2.0);
    ///
    /// let simbox = SimBox::from([4.0, 4.0, 4.0]);
    ///
    /// let distance = atom.distance_from_point(&point, Dimension::XY, &simbox).unwrap();
    /// assert_approx_eq!(f32, distance, 1.802776);
    /// ```
    pub fn distance_from_point(
        &self,
        point: &Vector3D,
        dim: Dimension,
        sbox: &SimBox,
    ) -> Result<f32, AtomError> {
        match self.position {
            None => Err(AtomError::InvalidPosition(PositionError::NoPosition(
                self.get_atom_number(),
            ))),
            Some(ref pos) => Ok(pos.distance(point, dim, sbox)),
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

    fn make_default_atom() -> Atom {
        Atom::new(45, "GLY", 123, "BB")
            .with_position([15.123, 14.321, 9.834].into())
            .with_velocity([-3.432, 0.184, 1.234].into())
            .with_force([5.1235, 2.3451, -0.32145].into())
    }

    #[test]
    fn new() {
        let atom = make_default_atom();

        assert_eq!(atom.get_residue_number(), 45);
        assert_eq!(atom.get_residue_name(), "GLY");
        assert_eq!(atom.get_atom_number(), 123);
        assert_eq!(atom.get_atom_name(), "BB");

        let pos = atom.get_position().unwrap();

        assert_approx_eq!(f32, pos.x, 15.123);
        assert_approx_eq!(f32, pos.y, 14.321);
        assert_approx_eq!(f32, pos.z, 9.834);

        let vel = atom.get_velocity().unwrap();

        assert_approx_eq!(f32, vel.x, -3.432);
        assert_approx_eq!(f32, vel.y, 0.184);
        assert_approx_eq!(f32, vel.z, 1.234);

        let force = atom.get_force().unwrap();

        assert_approx_eq!(f32, force.x, 5.1235);
        assert_approx_eq!(f32, force.y, 2.3451);
        assert_approx_eq!(f32, force.z, -0.32145);
    }

    #[test]
    fn residue_number() {
        let mut atom = make_default_atom();
        assert_eq!(atom.get_residue_number(), 45);

        atom.set_residue_number(187);
        assert_eq!(atom.get_residue_number(), 187);
    }

    #[test]
    fn residue_name() {
        let mut atom = make_default_atom();
        assert_eq!(atom.get_residue_name(), "GLY");

        atom.set_residue_name("LYS");
        assert_eq!(atom.get_residue_name(), "LYS");
    }

    #[test]
    fn atom_number() {
        let mut atom = make_default_atom();
        assert_eq!(atom.get_atom_number(), 123);

        atom.set_atom_number(13);
        assert_eq!(atom.get_atom_number(), 13);
    }

    #[test]
    fn atom_name() {
        let mut atom = make_default_atom();
        assert_eq!(atom.get_atom_name(), "BB");

        atom.set_atom_name("SC1");
        assert_eq!(atom.get_atom_name(), "SC1");
    }

    #[test]
    fn chain() {
        let mut atom = make_default_atom();
        assert_eq!(atom.get_chain(), None);

        atom.set_chain('B');
        assert_eq!(atom.get_chain().unwrap(), 'B');

        atom.reset_chain();
        assert_eq!(atom.get_chain(), None);

        let atom2 = make_default_atom().with_chain('B');
        assert_eq!(atom2.get_chain().unwrap(), 'B');
    }

    #[test]
    fn charge() {
        let mut atom = make_default_atom();
        assert_eq!(atom.get_charge(), None);

        atom.set_charge(0.453);
        assert_approx_eq!(f32, atom.get_charge().unwrap(), 0.453);

        atom.reset_charge();
        assert_eq!(atom.get_charge(), None);

        let atom2 = make_default_atom().with_charge(0.453);
        assert_approx_eq!(f32, atom2.get_charge().unwrap(), 0.453);
    }

    #[test]
    fn mass() {
        let mut atom = make_default_atom();
        assert_eq!(atom.get_mass(), None);

        atom.set_mass(10.453);
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 10.453);

        atom.reset_mass();
        assert_eq!(atom.get_mass(), None);

        let atom2 = make_default_atom().with_mass(10.453);
        assert_approx_eq!(f32, atom2.get_mass().unwrap(), 10.453);
    }

    #[test]
    fn vdw() {
        let mut atom = make_default_atom();
        assert_eq!(atom.get_vdw(), None);

        atom.set_vdw(0.19);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.19);

        atom.reset_vdw();
        assert_eq!(atom.get_vdw(), None);

        let atom2 = make_default_atom().with_vdw(0.19);
        assert_approx_eq!(f32, atom2.get_vdw().unwrap(), 0.19);
    }

    #[test]
    fn expected_max_bonds() {
        let mut atom = make_default_atom();
        assert_eq!(atom.get_expected_max_bonds(), None);

        atom.set_expected_max_bonds(3);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 3);

        atom.reset_expected_max_bonds();
        assert_eq!(atom.get_expected_max_bonds(), None);

        let atom2 = make_default_atom().with_expected_max_bonds(3);
        assert_eq!(atom2.get_expected_max_bonds().unwrap(), 3);
    }

    #[test]
    fn expected_min_bonds() {
        let mut atom = make_default_atom();
        assert_eq!(atom.get_expected_min_bonds(), None);

        atom.set_expected_min_bonds(3);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 3);

        atom.reset_expected_min_bonds();
        assert_eq!(atom.get_expected_min_bonds(), None);

        let atom2 = make_default_atom().with_expected_min_bonds(3);
        assert_eq!(atom2.get_expected_min_bonds().unwrap(), 3);
    }

    #[test]
    fn element_name() {
        let mut atom = make_default_atom();
        assert_eq!(atom.get_element_name(), None);

        atom.set_element_name("carbon");
        assert_eq!(atom.get_element_name().unwrap(), String::from("carbon"));

        atom.reset_element_name();
        assert_eq!(atom.get_element_name(), None);

        let atom2 = make_default_atom().with_element_name("carbon");
        assert_eq!(atom2.get_element_name().unwrap(), String::from("carbon"));
    }

    #[test]
    fn element_symbol() {
        let mut atom = make_default_atom();
        assert_eq!(atom.get_element_symbol(), None);

        atom.set_element_symbol("C");
        assert_eq!(atom.get_element_symbol().unwrap(), String::from("C"));

        atom.reset_element_symbol();
        assert_eq!(atom.get_element_symbol(), None);

        let atom2 = make_default_atom().with_element_symbol("C");
        assert_eq!(atom2.get_element_symbol().unwrap(), String::from("C"));
    }

    #[test]
    fn get_position() {
        let atom = make_default_atom();
        assert_eq!(
            *atom.get_position().unwrap(),
            Vector3D::new(15.123, 14.321, 9.834)
        );
    }

    #[test]
    fn get_position_mut() {
        let mut atom = make_default_atom();

        atom.get_position_mut().unwrap().x += 3.4;
        let pos = atom.get_position().unwrap();
        assert_approx_eq!(f32, pos.x, 18.523);
        assert_approx_eq!(f32, pos.y, 14.321);
        assert_approx_eq!(f32, pos.z, 9.834);
    }

    #[test]
    fn set_position() {
        let mut atom = make_default_atom();
        let new_pos = Vector3D::new(1.764, 2.134, 19.129);
        atom.set_position(new_pos.clone());

        assert_eq!(*atom.get_position().unwrap(), new_pos);
    }

    #[test]
    fn set_position_x() {
        let mut atom = make_default_atom();
        let new_x = 19.129;
        let new_pos = Vector3D::new(19.129, 14.321, 9.834);
        atom.set_position_x(new_x);

        assert_eq!(*atom.get_position().unwrap(), new_pos);
    }

    #[test]
    fn set_position_y() {
        let mut atom = make_default_atom();
        let new_y = 19.129;
        let new_pos = Vector3D::new(15.123, 19.129, 9.834);
        atom.set_position_y(new_y);

        assert_eq!(*atom.get_position().unwrap(), new_pos);
    }

    #[test]
    fn set_position_z() {
        let mut atom = make_default_atom();
        let new_z = 19.129;
        let new_pos = Vector3D::new(15.123, 14.321, 19.129);
        atom.set_position_z(new_z);

        assert_eq!(*atom.get_position().unwrap(), new_pos);
    }

    #[test]
    fn get_velocity() {
        let atom = make_default_atom();
        assert_eq!(
            *atom.get_velocity().unwrap(),
            Vector3D::new(-3.432, 0.184, 1.234)
        );
    }

    #[test]
    fn get_velocity_mut() {
        let mut atom = make_default_atom();

        atom.get_velocity_mut().unwrap().y += 0.3;
        let vel = atom.get_velocity().unwrap();
        assert_approx_eq!(f32, vel.x, -3.432);
        assert_approx_eq!(f32, vel.y, 0.484);
        assert_approx_eq!(f32, vel.z, 1.234);
    }

    #[test]
    fn set_velocity() {
        let mut atom = make_default_atom();
        let new_vel = Vector3D::new(1.764, 2.134, 19.129);
        atom.set_velocity(new_vel.clone());

        assert_eq!(*atom.get_velocity().unwrap(), new_vel);
    }

    #[test]
    fn get_force() {
        let atom = make_default_atom();
        assert_eq!(
            *atom.get_force().unwrap(),
            Vector3D::new(5.1235, 2.3451, -0.32145)
        );
    }

    #[test]
    fn get_force_mut() {
        let mut atom = make_default_atom();

        atom.get_force_mut().unwrap().z -= 0.13;
        let force = atom.get_force().unwrap();
        assert_approx_eq!(f32, force.x, 5.1235);
        assert_approx_eq!(f32, force.y, 2.3451);
        assert_approx_eq!(f32, force.z, -0.45145);
    }

    #[test]
    fn set_force() {
        let mut atom = make_default_atom();
        let new_for = Vector3D::new(1.764, 2.134, 19.129);
        atom.set_force(new_for.clone());

        assert_eq!(*atom.get_force().unwrap(), new_for);
    }

    #[test]
    fn get_set_bonded() {
        let mut atom = make_default_atom();
        assert_eq!(atom.get_n_bonded(), 0);

        let bonded = atom.get_bonded();
        assert_eq!(bonded, &AtomContainer::empty());

        unsafe {
            atom.set_bonded(vec![1, 2, 5]);
        }

        assert_eq!(atom.get_n_bonded(), 3);
        let bonded = atom.get_bonded();
        assert_eq!(bonded, &AtomContainer::from_indices(vec![1, 2, 5], 10));
    }

    #[test]
    fn add_bonded() {
        let mut atom = make_default_atom();
        assert_eq!(atom.get_n_bonded(), 0);

        unsafe {
            atom.add_bonded(1);
            atom.add_bonded(2);
            atom.add_bonded(5);
        }

        assert_eq!(atom.get_n_bonded(), 3);
        let bonded = atom.get_bonded();
        assert_eq!(bonded, &AtomContainer::from_indices(vec![1, 2, 5], 10));
    }

    #[test]
    fn has_position() {
        let mut atom = make_default_atom();
        assert!(atom.has_position());

        atom.reset_position();
        assert!(!atom.has_position());
    }

    #[test]
    fn has_velocity() {
        let mut atom = make_default_atom();
        assert!(atom.has_velocity());

        atom.reset_velocity();
        assert!(!atom.has_velocity());
    }

    #[test]
    fn has_force() {
        let mut atom = make_default_atom();
        assert!(atom.has_force());

        atom.reset_force();
        assert!(!atom.has_force());
    }

    #[test]
    fn translate_nopbc() {
        let mut atom = make_default_atom();

        let shift = Vector3D::new(4.5, 2.3, -8.3);
        atom.translate_nopbc(&shift).unwrap();

        assert_approx_eq!(
            f32,
            atom.get_position().unwrap().x,
            19.623,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom.get_position().unwrap().y,
            16.621,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom.get_position().unwrap().z,
            1.534,
            epsilon = 0.00001
        );
    }

    #[test]
    fn translate_nopbc_fail() {
        let mut atom = make_default_atom();
        atom.reset_position();

        let shift = Vector3D::new(4.5, 2.3, -8.3);
        match atom.translate_nopbc(&shift) {
            Ok(_) => panic!("Function should have failed."),
            Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => {
                assert_eq!(x, atom.get_atom_number())
            }
            Err(e) => {
                panic!("Function failed successfully, but incorrect error type `{e}` was returned.")
            }
        }
    }

    #[test]
    fn translate() {
        let mut atom = make_default_atom();

        let shift = Vector3D::new(4.5, 2.3, -10.2);
        let simbox = SimBox::from([16.0, 16.0, 16.0]);

        atom.translate(&shift, &simbox).unwrap();

        assert_approx_eq!(
            f32,
            atom.get_position().unwrap().x,
            3.623,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom.get_position().unwrap().y,
            0.621,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom.get_position().unwrap().z,
            15.634,
            epsilon = 0.00001
        );
    }

    #[test]
    fn translate_fail() {
        let mut atom = make_default_atom();
        atom.reset_position();

        let shift = Vector3D::new(4.5, 2.3, -8.3);
        let simbox = SimBox::from([16.0, 16.0, 16.0]);
        match atom.translate(&shift, &simbox) {
            Ok(_) => panic!("Function should have failed."),
            Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => {
                assert_eq!(x, atom.get_atom_number())
            }
            Err(e) => {
                panic!("Function failed successfully, but incorrect error type `{e}` was returned.")
            }
        }
    }

    #[test]
    fn wrap() {
        let mut atom =
            Atom::new(45, "GLY", 123, "BB").with_position([15.123, 14.321, -1.743].into());

        let simbox = SimBox::from([15.0, 15.0, 15.0]);

        atom.wrap(&simbox).unwrap();

        assert_approx_eq!(
            f32,
            atom.get_position().unwrap().x,
            0.123,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom.get_position().unwrap().y,
            14.321,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom.get_position().unwrap().z,
            13.257,
            epsilon = 0.00001
        );
    }

    #[test]
    fn wrap_far() {
        let mut atom =
            Atom::new(45, "GLY", 123, "BB").with_position([60.123, 14.321, -31.743].into());

        let simbox = SimBox::from([15.0, 15.0, 15.0]);

        atom.wrap(&simbox).unwrap();

        assert_approx_eq!(
            f32,
            atom.get_position().unwrap().x,
            0.123,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom.get_position().unwrap().y,
            14.321,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom.get_position().unwrap().z,
            13.257,
            epsilon = 0.00001
        );
    }

    #[test]
    fn wrap_fail() {
        let mut atom = make_default_atom();
        atom.reset_position();

        let simbox = SimBox::from([15.0, 15.0, 15.0]);
        match atom.wrap(&simbox) {
            Ok(_) => panic!("Function should have failed."),
            Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => {
                assert_eq!(x, atom.get_atom_number())
            }
            Err(e) => {
                panic!("Function failed successfully, but incorrect error type `{e}` was returned.")
            }
        }
    }

    #[test]
    fn distance_x() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom1.distance(&atom2, Dimension::X, &simbox).unwrap(),
            -0.7,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::X, &simbox).unwrap(),
            0.7,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_x() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());

        assert_approx_eq!(
            f32,
            atom1.distance_naive(&atom2, Dimension::X).unwrap(),
            3.3,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance_naive(&atom1, Dimension::X).unwrap(),
            -3.3,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_y() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom1.distance(&atom2, Dimension::Y, &simbox).unwrap(),
            1.0,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::Y, &simbox).unwrap(),
            -1.0,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_y() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());

        assert_approx_eq!(
            f32,
            atom1.distance_naive(&atom2, Dimension::Y).unwrap(),
            1.0,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance_naive(&atom1, Dimension::Y).unwrap(),
            -1.0,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_z() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 1.0].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 3.5].into());
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom1.distance(&atom2, Dimension::Z, &simbox).unwrap(),
            1.5,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::Z, &simbox).unwrap(),
            -1.5,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_z() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());

        assert_approx_eq!(
            f32,
            atom1.distance_naive(&atom2, Dimension::Z).unwrap(),
            2.5,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance_naive(&atom1, Dimension::Z).unwrap(),
            -2.5,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_xy() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom1.distance(&atom2, Dimension::XY, &simbox).unwrap(),
            1.2206556,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::XY, &simbox).unwrap(),
            1.2206556,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_xy() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());

        assert_approx_eq!(
            f32,
            atom1.distance_naive(&atom2, Dimension::XY).unwrap(),
            3.448188,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance_naive(&atom1, Dimension::XY).unwrap(),
            3.448188,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_xz() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom1.distance(&atom2, Dimension::XZ, &simbox).unwrap(),
            1.6552945,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::XZ, &simbox).unwrap(),
            1.6552945,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_xz() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());

        assert_approx_eq!(
            f32,
            atom1.distance_naive(&atom2, Dimension::XZ).unwrap(),
            4.140048,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance_naive(&atom1, Dimension::XZ).unwrap(),
            4.140048,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_yz() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom1.distance(&atom2, Dimension::YZ, &simbox).unwrap(),
            1.8027756,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::YZ, &simbox).unwrap(),
            1.8027756,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_yz() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());

        assert_approx_eq!(
            f32,
            atom1.distance_naive(&atom2, Dimension::YZ).unwrap(),
            2.692582,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance_naive(&atom1, Dimension::YZ).unwrap(),
            2.692582,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_xyz() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom1.distance(&atom2, Dimension::XYZ, &simbox).unwrap(),
            1.933908,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::XYZ, &simbox).unwrap(),
            1.933908,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_xyz() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());

        assert_approx_eq!(
            f32,
            atom1.distance_naive(&atom2, Dimension::XYZ).unwrap(),
            4.259108,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance_naive(&atom1, Dimension::XYZ).unwrap(),
            4.259108,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_none() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom1.distance(&atom2, Dimension::None, &simbox).unwrap(),
            0.0,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::None, &simbox).unwrap(),
            0.0,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_none() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());
        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());

        assert_approx_eq!(
            f32,
            atom1.distance_naive(&atom2, Dimension::None).unwrap(),
            0.0,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance_naive(&atom1, Dimension::None).unwrap(),
            0.0,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_fail() {
        let mut atom1 = Atom::new(1, "LYS", 1, "BB");
        let mut atom2 = Atom::new(1, "LYS", 2, "SC1");
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        match atom1.distance(&atom2, Dimension::XYZ, &simbox) {
            Ok(_) => panic!("Function should have failed."),
            Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => {
                assert_eq!(x, atom1.get_atom_number())
            }
            Err(e) => {
                panic!("Function failed successfully, but incorrect error type `{e}` was returned.")
            }
        }

        atom2.set_position([1.0, 2.0, 3.0].into());

        match atom1.distance(&atom2, Dimension::XYZ, &simbox) {
            Ok(_) => panic!("Function should have failed."),
            Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => {
                assert_eq!(x, atom1.get_atom_number())
            }
            Err(e) => {
                panic!("Function failed successfully, but incorrect error type `{e}` was returned.")
            }
        }

        atom2.reset_position();
        atom1.set_position([1.0, 2.0, 3.0].into());

        match atom1.distance(&atom2, Dimension::XYZ, &simbox) {
            Ok(_) => panic!("Function should have failed."),
            Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => {
                assert_eq!(x, atom2.get_atom_number())
            }
            Err(e) => {
                panic!("Function failed successfully, but incorrect error type `{e}` was returned.")
            }
        }
    }

    #[test]
    fn distance_naive_fail() {
        let mut atom1 = Atom::new(1, "LYS", 1, "BB");
        let mut atom2 = Atom::new(1, "LYS", 2, "SC1");

        match atom1.distance_naive(&atom2, Dimension::XYZ) {
            Ok(_) => panic!("Function should have failed."),
            Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => {
                assert_eq!(x, atom1.get_atom_number())
            }
            Err(e) => {
                panic!("Function failed successfully, but incorrect error type `{e}` was returned.")
            }
        }

        atom2.set_position([1.0, 2.0, 3.0].into());

        match atom1.distance_naive(&atom2, Dimension::XYZ) {
            Ok(_) => panic!("Function should have failed."),
            Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => {
                assert_eq!(x, atom1.get_atom_number())
            }
            Err(e) => {
                panic!("Function failed successfully, but incorrect error type `{e}` was returned.")
            }
        }

        atom2.reset_position();
        atom1.set_position([1.0, 2.0, 3.0].into());

        match atom1.distance_naive(&atom2, Dimension::XYZ) {
            Ok(_) => panic!("Function should have failed."),
            Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => {
                assert_eq!(x, atom2.get_atom_number())
            }
            Err(e) => {
                panic!("Function failed successfully, but incorrect error type `{e}` was returned.")
            }
        }
    }

    #[test]
    fn distance_from_point_x() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::new(1.0, 0.5, 2.0);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::X, &simbox)
                .unwrap(),
            1.0,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_y() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::new(1.0, 0.5, 2.0);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::Y, &simbox)
                .unwrap(),
            -0.7,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_z() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::new(1.0, 0.5, 2.0);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::Z, &simbox)
                .unwrap(),
            1.5,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_xy() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::new(1.0, 0.5, 2.0);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::XY, &simbox)
                .unwrap(),
            1.220656,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_xz() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::new(1.0, 0.5, 2.0);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::XZ, &simbox)
                .unwrap(),
            1.802776,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_yz() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::new(1.0, 0.5, 2.0);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::YZ, &simbox)
                .unwrap(),
            1.6552945,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_xyz() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::new(1.0, 0.5, 2.0);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::XYZ, &simbox)
                .unwrap(),
            1.933908,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_none() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::new(1.0, 0.5, 2.0);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::None, &simbox)
                .unwrap(),
            0.0,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_fail() {
        let atom = Atom::new(1, "LYS", 1, "BB");
        let point = Vector3D::new(1.0, 0.5, 2.0);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        match atom.distance_from_point(&point, Dimension::XYZ, &simbox) {
            Ok(_) => panic!("Function should have failed."),
            Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => {
                assert_eq!(x, atom.get_atom_number())
            }
            Err(e) => {
                panic!("Function failed successfully, but incorrect error type `{e}` was returned.")
            }
        }
    }
}

#[cfg(test)]
#[cfg(feature = "serde")]
mod serde_tests {
    use std::fs::read_to_string;

    use float_cmp::assert_approx_eq;

    use super::*;

    #[test]
    fn atom_to_yaml() {
        let mut atom = Atom::new(10, "CYS", 24, "CA")
            .with_position(Vector3D::new(10.4, 7.5, 2.6))
            .with_velocity(Vector3D::new(-0.4332, 0.3221, 1.2314))
            .with_force(Vector3D::new(14.24898, -13.208472, 8.123864))
            .with_chain('B')
            .with_charge(-1.0)
            .with_mass(32.12)
            .with_vdw(0.64)
            .with_expected_max_bonds(3)
            .with_expected_min_bonds(1)
            .with_element_name("carbon")
            .with_element_symbol("C");

        unsafe {
            atom.add_bonded(4);
            atom.add_bonded(7);
            atom.add_bonded(8);
            atom.add_bonded(9);
        }

        let string = serde_yaml::to_string(&atom).unwrap();
        let expected = read_to_string("test_files/serde_atom.yaml").unwrap();

        assert_eq!(string, expected);
    }

    #[test]
    fn minimal_atom_to_yaml() {
        let atom = Atom::new(10, "CYS", 24, "CA");

        let string = serde_yaml::to_string(&atom).unwrap();
        let expected = read_to_string("test_files/serde_atom_minimal.yaml").unwrap();

        assert_eq!(string, expected);
    }

    #[test]
    fn atom_from_yaml() {
        let string = read_to_string("test_files/serde_atom.yaml").unwrap();
        let atom: Atom = serde_yaml::from_str(&string).unwrap();

        assert_eq!(atom.get_residue_number(), 10);
        assert_eq!(atom.get_residue_name(), "CYS");
        assert_eq!(atom.get_atom_number(), 24);
        assert_eq!(atom.get_atom_name(), "CA");

        let pos = atom.get_position().unwrap();

        assert_approx_eq!(f32, pos.x, 10.4);
        assert_approx_eq!(f32, pos.y, 7.5);
        assert_approx_eq!(f32, pos.z, 2.6);

        let vel = atom.get_velocity().unwrap();

        assert_approx_eq!(f32, vel.x, -0.4332);
        assert_approx_eq!(f32, vel.y, 0.3221);
        assert_approx_eq!(f32, vel.z, 1.2314);

        let force = atom.get_force().unwrap();

        assert_approx_eq!(f32, force.x, 14.24898);
        assert_approx_eq!(f32, force.y, -13.208472);
        assert_approx_eq!(f32, force.z, 8.123864);

        assert_eq!(atom.get_chain().unwrap(), 'B');
        assert_approx_eq!(f32, atom.get_charge().unwrap(), -1.0);
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 32.12);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.64);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 3);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);
        assert_eq!(atom.get_element_name().unwrap(), "carbon");
        assert_eq!(atom.get_element_symbol().unwrap(), "C");

        assert_eq!(
            atom.get_bonded(),
            &AtomContainer::from_indices(vec![4, 7, 8, 9], 100)
        );
    }

    #[test]
    fn minimal_atom_from_yaml() {
        let string = read_to_string("test_files/serde_atom_minimal.yaml").unwrap();
        let atom: Atom = serde_yaml::from_str(&string).unwrap();

        assert_eq!(atom.get_residue_number(), 10);
        assert_eq!(atom.get_residue_name(), "CYS");
        assert_eq!(atom.get_atom_number(), 24);
        assert_eq!(atom.get_atom_name(), "CA");

        assert!(atom.get_position().is_none());
        assert!(atom.get_velocity().is_none());
        assert!(atom.get_force().is_none());
        assert!(atom.get_chain().is_none());
        assert!(atom.get_charge().is_none());
        assert!(atom.get_mass().is_none());
        assert!(atom.get_vdw().is_none());
        assert!(atom.get_expected_max_bonds().is_none());
        assert!(atom.get_expected_min_bonds().is_none());
        assert!(atom.get_element_name().is_none());
        assert!(atom.get_element_symbol().is_none());

        assert_eq!(atom.get_bonded(), &AtomContainer::from_indices(vec![], 100));
    }

    #[test]
    fn atom_from_yaml_unknown_field() {
        let string = read_to_string("test_files/serde_atom_unknown_field.yaml").unwrap();
        let err = serde_yaml::from_str::<Atom>(&string).unwrap_err();
        assert!(err.to_string().contains("unknown field"));
    }
}
