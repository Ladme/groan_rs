// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of the Atom structure and its methods.

use std::io::Write;

use crate::errors::{WriteGroError, WritePdbError};
use crate::structures::{
    container::AtomContainer, dimension::Dimension, simbox::SimBox, vector3d::Vector3D,
};

#[derive(Debug, Clone)]
pub struct Atom {
    /// Number (index) of the residue this atom is part of.
    residue_number: usize,
    /// Name of the residue this atom is part of.
    residue_name: String,
    /// Number (index) of the atom.
    atom_number: usize,
    /// Name of the atom.
    atom_name: String,
    /// Identifier of the chain this atom is part of. (Optional.)
    chain: Option<char>,
    /// Charge of the atom. (Optional.)
    charge: Option<f32>,
    /// Position of the atom in 3D space. (Optional.)
    position: Option<Vector3D>,
    /// Velocity of the atom. (Optional.)
    velocity: Option<Vector3D>,
    /// Force acting on the atom. (Optional.)
    force: Option<Vector3D>,
    /// Indices of atoms that are bonded to this atom.
    bonded: AtomContainer,
}

impl Atom {
    /// Create new Atom structure with the specified properties.
    ///
    /// ## Notes
    /// - By default, atom is constructed with `position`, `velocity`, `force`, `chain`, and `charge`
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

    /// Get the charge of the atom.
    pub fn get_charge(&self) -> Option<f32> {
        self.charge
    }

    /// Set the charge of the atom.
    pub fn set_charge(&mut self, charge: f32) {
        self.charge = Some(charge);
    }

    /// Get the coordinates of the atom.
    pub fn get_position(&self) -> Option<&Vector3D> {
        self.position.as_ref()
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
            None => self.position = Some(Vector3D::from([x, 0.0, 0.0])),
            Some(ref mut pos) => pos.x = x,
        }
    }

    /// Set the y-coordinate of the atom.
    pub fn set_position_y(&mut self, y: f32) {
        match self.position {
            None => self.position = Some(Vector3D::from([0.0, y, 0.0])),
            Some(ref mut pos) => pos.y = y,
        }
    }

    /// Set the z-coordinate of the atom.
    pub fn set_position_z(&mut self, z: f32) {
        match self.position {
            None => self.position = Some(Vector3D::from([0.0, 0.0, z])),
            Some(ref mut pos) => pos.z = z,
        }
    }

    /// Get the velocity vector of the atom.
    pub fn get_velocity(&self) -> Option<&Vector3D> {
        self.velocity.as_ref()
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
    /// ## Panics
    /// Panics if atom has no position.
    pub fn translate(&mut self, translate: &Vector3D, sbox: &SimBox) {
        if let Some(ref mut pos) = self.position {
            pos.x += translate.x;
            pos.y += translate.y;
            pos.z += translate.z;

            pos.wrap(sbox);
        } else {
            panic!("FATAL GROAN ERROR | Atom::translate | Atom has no position.");
        }
    }

    /// Translates the position of the atom by the provided Vector3D.
    /// Does **not** wrap the atom to the simulation box.
    ///
    /// ## Panics
    /// Panics if atom has no position.
    pub fn translate_nopbc(&mut self, translate: &Vector3D) {
        if let Some(ref mut pos) = self.position {
            pos.x += translate.x;
            pos.y += translate.y;
            pos.z += translate.z;
        } else {
            panic!("FATAL GROAN ERROR | Atom::translate_nopbc | Atom has no position.");
        }
    }

    /// Wrap the atom into the simulation box.
    ///
    /// ## Panics
    /// Panics if atom has no position.
    pub fn wrap(&mut self, sbox: &SimBox) {
        match self.position {
            None => panic!("FATAL GROAN ERROR | Atom::wrap | Atom has no position."),
            Some(ref mut pos) => pos.wrap(sbox),
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

    /// Calculate distance between two atoms.
    /// Takes periodic boundary conditions into consideration.
    /// Returns oriented distance for 1D problems.
    ///
    /// ## Warning
    /// - Currently only works with orthogonal simulation boxes.
    ///
    /// ## Panics
    /// Panics if any of the atoms has no position.
    ///
    /// ## Example
    /// Calculate distance between two atoms in the xy-plane.
    /// ```
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([1.0, 2.0, 3.0].into());
    /// let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([3.5, 1.0, 2.0].into());
    ///
    /// let simbox = SimBox::from([4.0, 4.0, 4.0]);
    ///
    /// let distance = atom1.distance(&atom2, Dimension::XY, &simbox);
    /// assert_approx_eq!(f32, distance, 1.802776);
    /// ```
    pub fn distance(&self, atom: &Atom, dim: Dimension, sbox: &SimBox) -> f32 {
        match (&self.position, &atom.position) {
            (None | Some(_), None) | (None, Some(_)) => {
                panic!("FATAL GROAN ERROR | Atom::distance | Atom has no position.")
            }
            (Some(ref pos1), Some(ref pos2)) => pos1.distance(pos2, dim, sbox),
        }
    }

    /// Calculate distance between an atom and a point in space.
    /// Takes periodic boundary conditions into consideration.
    /// Returns oriented distance for 1D problems.
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
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let atom = Atom::new(1, "LYS", 1, "BB").with_position([1.0, 2.0, 3.0].into());
    /// let point = Vector3D::from([3.5, 1.0, 2.0]);
    ///
    /// let simbox = SimBox::from([4.0, 4.0, 4.0]);
    ///
    /// let distance = atom.distance_from_point(&point, Dimension::XY, &simbox);
    /// assert_approx_eq!(f32, distance, 1.802776);
    /// ```
    pub fn distance_from_point(&self, point: &Vector3D, dim: Dimension, sbox: &SimBox) -> f32 {
        match self.position {
            None => panic!("FATAL GROAN ERROR | Atom::distance_from_point | Atom has no position."),
            Some(ref pos) => pos.distance(point, dim, sbox),
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

        assert_approx_eq!(f32, atom.get_position().unwrap().x, 15.123);
        assert_approx_eq!(f32, atom.get_position().unwrap().y, 14.321);
        assert_approx_eq!(f32, atom.get_position().unwrap().z, 9.834);

        assert_approx_eq!(f32, atom.get_velocity().unwrap().x, -3.432);
        assert_approx_eq!(f32, atom.get_velocity().unwrap().y, 0.184);
        assert_approx_eq!(f32, atom.get_velocity().unwrap().z, 1.234);

        assert_approx_eq!(f32, atom.get_force().unwrap().x, 5.1235);
        assert_approx_eq!(f32, atom.get_force().unwrap().y, 2.3451);
        assert_approx_eq!(f32, atom.get_force().unwrap().z, -0.32145);
    }

    #[test]
    fn get_residue_number() {
        let atom = make_default_atom();
        assert_eq!(atom.get_residue_number(), 45);
    }

    #[test]
    fn set_residue_number() {
        let mut atom = make_default_atom();
        atom.set_residue_number(187);
        assert_eq!(atom.get_residue_number(), 187);
    }

    #[test]
    fn get_residue_name() {
        let atom = make_default_atom();
        assert_eq!(atom.get_residue_name(), "GLY");
    }

    #[test]
    fn set_residue_name() {
        let mut atom = make_default_atom();
        atom.set_residue_name("LYS");
        assert_eq!(atom.get_residue_name(), "LYS");
    }

    #[test]
    fn get_atom_number() {
        let atom = make_default_atom();
        assert_eq!(atom.get_atom_number(), 123);
    }

    #[test]
    fn set_atom_number() {
        let mut atom = make_default_atom();
        atom.set_atom_number(13);
        assert_eq!(atom.get_atom_number(), 13);
    }

    #[test]
    fn get_atom_name() {
        let atom = make_default_atom();
        assert_eq!(atom.get_atom_name(), "BB");
    }

    #[test]
    fn set_atom_name() {
        let mut atom = make_default_atom();
        atom.set_atom_name("SC1");
        assert_eq!(atom.get_atom_name(), "SC1");
    }

    #[test]
    fn get_position() {
        let atom = make_default_atom();
        assert_eq!(
            *atom.get_position().unwrap(),
            Vector3D::from([15.123, 14.321, 9.834])
        );
    }

    #[test]
    fn set_position() {
        let mut atom = make_default_atom();
        let new_pos = Vector3D::from([1.764, 2.134, 19.129]);
        atom.set_position(new_pos.clone());

        assert_eq!(*atom.get_position().unwrap(), new_pos);
    }

    #[test]
    fn set_position_x() {
        let mut atom = make_default_atom();
        let new_x = 19.129;
        let new_pos = Vector3D::from([19.129, 14.321, 9.834]);
        atom.set_position_x(new_x);

        assert_eq!(*atom.get_position().unwrap(), new_pos);
    }

    #[test]
    fn set_position_y() {
        let mut atom = make_default_atom();
        let new_y = 19.129;
        let new_pos = Vector3D::from([15.123, 19.129, 9.834]);
        atom.set_position_y(new_y);

        assert_eq!(*atom.get_position().unwrap(), new_pos);
    }

    #[test]
    fn set_position_z() {
        let mut atom = make_default_atom();
        let new_z = 19.129;
        let new_pos = Vector3D::from([15.123, 14.321, 19.129]);
        atom.set_position_z(new_z);

        assert_eq!(*atom.get_position().unwrap(), new_pos);
    }

    #[test]
    fn get_velocity() {
        let atom = make_default_atom();
        assert_eq!(
            *atom.get_velocity().unwrap(),
            Vector3D::from([-3.432, 0.184, 1.234])
        );
    }

    #[test]
    fn set_velocity() {
        let mut atom = make_default_atom();
        let new_vel = Vector3D::from([1.764, 2.134, 19.129]);
        atom.set_velocity(new_vel.clone());

        assert_eq!(*atom.get_velocity().unwrap(), new_vel);
    }

    #[test]
    fn get_force() {
        let atom = make_default_atom();
        assert_eq!(
            *atom.get_force().unwrap(),
            Vector3D::from([5.1235, 2.3451, -0.32145])
        );
    }

    #[test]
    fn set_force() {
        let mut atom = make_default_atom();
        let new_for = Vector3D::from([1.764, 2.134, 19.129]);
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

        let shift = Vector3D::from([4.5, 2.3, -8.3]);
        atom.translate_nopbc(&shift);

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
    fn translate() {
        let mut atom = make_default_atom();

        let shift = Vector3D::from([4.5, 2.3, -10.2]);
        let simbox = SimBox::from([16.0, 16.0, 16.0]);

        atom.translate(&shift, &simbox);

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
    fn wrap() {
        let mut atom =
            Atom::new(45, "GLY", 123, "BB").with_position([15.123, 14.321, -1.743].into());

        let simbox = SimBox::from([15.0, 15.0, 15.0]);

        atom.wrap(&simbox);

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

        atom.wrap(&simbox);

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
    fn distance_x() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([3.8, 2.0, 3.5].into());

        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([0.5, 1.0, 1.0].into());

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom1.distance(&atom2, Dimension::X, &simbox),
            -0.7,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::X, &simbox),
            0.7,
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
            atom1.distance(&atom2, Dimension::Y, &simbox),
            1.0,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::Y, &simbox),
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
            atom1.distance(&atom2, Dimension::Z, &simbox),
            1.5,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::Z, &simbox),
            -1.5,
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
            atom1.distance(&atom2, Dimension::XY, &simbox),
            1.2206556,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::XY, &simbox),
            1.2206556,
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
            atom1.distance(&atom2, Dimension::XZ, &simbox),
            1.6552945,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::XZ, &simbox),
            1.6552945,
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
            atom1.distance(&atom2, Dimension::YZ, &simbox),
            1.8027756,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::YZ, &simbox),
            1.8027756,
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
            atom1.distance(&atom2, Dimension::XYZ, &simbox),
            1.933908,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::XYZ, &simbox),
            1.933908,
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
            atom1.distance(&atom2, Dimension::None, &simbox),
            0.0,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            atom2.distance(&atom1, Dimension::None, &simbox),
            0.0,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_x() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::from([1.0, 0.5, 2.0]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::X, &simbox),
            1.0,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_y() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::from([1.0, 0.5, 2.0]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::Y, &simbox),
            -0.7,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_z() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::from([1.0, 0.5, 2.0]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::Z, &simbox),
            1.5,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_xy() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::from([1.0, 0.5, 2.0]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::XY, &simbox),
            1.220656,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_xz() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::from([1.0, 0.5, 2.0]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::XZ, &simbox),
            1.802776,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_yz() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::from([1.0, 0.5, 2.0]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::YZ, &simbox),
            1.6552945,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_xyz() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::from([1.0, 0.5, 2.0]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::XYZ, &simbox),
            1.933908,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_from_point_none() {
        let atom = Atom::new(1, "LYS", 1, "BB").with_position([2.0, 3.8, 3.5].into());

        let point = Vector3D::from([1.0, 0.5, 2.0]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            atom.distance_from_point(&point, Dimension::None, &simbox),
            0.0,
            epsilon = 0.00001
        );
    }
}
