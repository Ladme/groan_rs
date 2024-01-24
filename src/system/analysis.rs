// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of various methods for analysis of `System`.

use crate::errors::{AtomError, GroupError, MassError, PositionError};
use crate::structures::simbox::{simbox_check, SimBox};
use crate::structures::{dimension::Dimension, vector3d::Vector3D};
use crate::system::general::System;

use std::f32::consts;

// PI times 2.
const PI_X2: f32 = consts::PI * 2.0f32;

/// ## Methods for analyzing the properties of the system.
impl System {
    /// Calculate center of geometry of a group in `System`.
    /// Takes periodic boundary conditions into consideration.
    ///
    /// ## Returns
    /// - `Vector3D` corresponding to the geometric center of the group.
    /// - `GroupError::NotFound` if the group does not exist.
    /// - `GroupError::EmptyGroup` if the group contains no atoms.
    /// - `GroupError::InvalidSimBox` if the system has no simulation box
    /// or the simulation box is not orthogonal.
    /// - `GroupError::InvalidPosition` if any of the atoms in the group has no position.
    ///
    /// ## Notes
    /// - This calculation approach is adapted from Linge Bai & David Breen (2008).
    /// - It is able to calculate correct center of geometry for any distribution of atoms
    /// that is not completely homogeneous.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// // calculate center of geometry for group "Group"
    /// let center = match system.group_get_center("Group") {
    ///     Ok(x) => x,
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;    
    ///     }
    /// };
    /// ```
    pub fn group_get_center(&self, name: &str) -> Result<Vector3D, GroupError> {
        // we can't work with a nonexistent or empty group
        if self.group_isempty(name)? {
            return Err(GroupError::EmptyGroup(name.to_string()));
        }

        let simbox = simbox_check(self.get_box_as_ref()).map_err(GroupError::InvalidSimBox)?;

        let scaling = Vector3D::new(PI_X2 / simbox.x, PI_X2 / simbox.y, PI_X2 / simbox.z);

        let mut sum_xi = Vector3D::default();
        let mut sum_zeta = Vector3D::default();

        for atom in self.group_iter(name)
            .expect("FATAL GROAN ERROR | System::group_get_center | Group not found but this should have already been checked.") 
        {
            match atom.get_position() {
                Some(x) => center_atom_contribution(
                    x.clone(),
                    &scaling,
                    simbox,
                    1.0, // use 1 as mass since we are calculating center of geometry
                    &mut sum_xi,
                    &mut sum_zeta,
                ),
                None => {
                    return Err(GroupError::InvalidPosition(PositionError::NoPosition(
                        atom.get_atom_number(),
                    )));
                }
            }
        }

        // convert to real coordinates
        Ok(from_circle_to_line(sum_zeta, sum_xi, &scaling))
    }

    /// Calculate center of mass of a group in `System`.
    /// Takes periodic boundary conditions into consideration.
    ///
    /// ## Returns
    /// - `Vector3D` corresponding to the center of mass of the group.
    /// - `GroupError::NotFound` if the group does not exist.
    /// - `GroupError::EmptyGroup` if the group contains no atoms.
    /// - `GroupError::InvalidSimBox` if the system has no simulation box
    /// or the simulation box is not orthogonal.
    /// - `GroupError::InvalidPosition` if any of the atoms in the group has no position.
    /// - `GroupError::InvalidMass` if any of the atoms in the group has no mass.
    ///
    /// ## Notes
    /// - This calculation approach is adapted from Linge Bai & David Breen (2008).
    /// - It is able to calculate correct center of mass for any distribution of atoms
    /// that is not completely homogeneous.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// // ... assign masses to atoms ...
    ///
    /// // calculate center of mass for group "Group"
    /// let center = match system.group_get_com("Group") {
    ///     Ok(x) => x,
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;    
    ///     }
    /// };
    /// ```
    pub fn group_get_com(&self, name: &str) -> Result<Vector3D, GroupError> {
        // we can't work with a nonexistent or empty group
        if self.group_isempty(name)? {
            return Err(GroupError::EmptyGroup(name.to_string()));
        }

        let simbox = simbox_check(self.get_box_as_ref()).map_err(GroupError::InvalidSimBox)?;

        let scaling = Vector3D::new(PI_X2 / simbox.x, PI_X2 / simbox.y, PI_X2 / simbox.z);

        let mut sum_xi = Vector3D::default();
        let mut sum_zeta = Vector3D::default();

        for atom in self.group_iter(name)
            .expect("FATAL GROAN ERROR | System::group_get_com | Group not found but this should have already been checked.") 
        {
            let mass = atom
                .get_mass()
                .ok_or(GroupError::InvalidMass(MassError::NoMass(
                    atom.get_atom_number(),
                )))?;

            match atom.get_position() {
                Some(x) => center_atom_contribution(
                    x.clone(),
                    &scaling,
                    simbox,
                    mass,
                    &mut sum_xi,
                    &mut sum_zeta,
                ),
                None => {
                    return Err(GroupError::InvalidPosition(PositionError::NoPosition(
                        atom.get_atom_number(),
                    )));
                }
            }
        }

        Ok(from_circle_to_line(sum_zeta, sum_xi, &scaling))
    }

    /// Calculate distance between the centers of geometries of the specified groups.
    /// Oriented distance is used for 1D problems.
    ///
    /// ## Returns
    /// - `f32` corresponding to the distance between the groups.
    /// - `GroupError::NotFound` if any of the groups does not exist.
    /// - `GroupError::InvalidSimBox` if the system has no simulation box
    /// or the simulation box is not orthogonal.
    /// - `GroupError::InvalidPosition` if any of the atoms in the groups has no position.
    ///
    /// ## Example
    /// Calculating the distance between the group 'Protein' and 'Membrane' along the z-dimension.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// let distance = match system.group_distance("Protein", "Membrane", Dimension::Z) {
    ///     Ok(x) => x,
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;
    ///     }
    /// };
    /// ```
    pub fn group_distance(
        &self,
        group1: &str,
        group2: &str,
        dim: Dimension,
    ) -> Result<f32, GroupError> {
        let group1_center = self.group_get_center(group1)?;
        let group2_center = self.group_get_center(group2)?;

        let simbox = simbox_check(self.get_box_as_ref()).map_err(GroupError::InvalidSimBox)?;

        Ok(group1_center.distance(&group2_center, dim, simbox))
    }

    /// Calculate distances between each pair of atoms of the specified groups.
    /// Oriented distances are used for 1D problems.
    ///
    /// ## Returns
    /// - Two dimensional vector of distances.
    /// - `GroupError::NotFound` if any of the groups does not exist.
    /// - `GroupError::InvalidSimBox` if the system has no simulation box
    /// or the simulation box is not orthogonal.
    /// - `GroupError::InvalidPosition` if any of the atoms in the groups has no position.
    ///
    /// ## Example
    /// Calculate distances between atoms of group 'Protein' and 'Membrane'.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// // calculate the matrix of distances
    /// let distances: Vec<Vec<f32>> = match system.group_all_distances("Protein", "Membrane", Dimension::XYZ) {
    ///     Ok(x) => x,
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;
    ///     }
    /// };
    ///
    /// // get the maximal distance between the atoms
    /// let max = distances
    ///     .iter()
    ///     .flatten()
    ///     .fold(std::f32::NEG_INFINITY, |max, &current| max.max(current));
    /// // get the minimal distance between the atoms
    /// let min = distances
    ///     .iter()
    ///     .flatten()
    ///     .fold(std::f32::INFINITY, |min, &current| min.min(current));
    /// ```
    pub fn group_all_distances(
        &self,
        group1: &str,
        group2: &str,
        dim: Dimension,
    ) -> Result<Vec<Vec<f32>>, GroupError> {
        let n_atoms_group1 = self.group_get_n_atoms(group1)?;
        let n_atoms_group2 = self.group_get_n_atoms(group2)?;

        let simbox = simbox_check(self.get_box_as_ref()).map_err(GroupError::InvalidSimBox)?;

        let mut distances = Vec::with_capacity(n_atoms_group1);

        for (i, atom1) in self.group_iter(group1)?.enumerate() {
            distances.push(Vec::with_capacity(n_atoms_group2));
            for atom2 in self.group_iter(group2)? {
                let dist = match atom1.distance(atom2, dim, simbox) {
                    Ok(x) => x,
                    Err(AtomError::InvalidPosition(e)) => return Err(GroupError::InvalidPosition(e)),
                    _ => panic!("FATAL GROAN ERROR | System::group_all_distances | Invalid error type returned by Atom::distance."),
                };

                distances[i].push(dist);
            }
        }

        Ok(distances)
    }

    /// Calculate distance between two atoms in the system.
    /// Assumes orthogonal periodic boundary conditions.
    /// Atoms are indexed starting from 0.
    /// Oriented distance is used for 1D problems.
    ///
    /// ## Returns
    /// - `f32` corresponding to the distance between the atoms.
    /// - `AtomError::OutOfRange` if any of the atom indices is out of range.
    /// - `AtomError::InvalidSimBox` if the system has no simulation box
    /// or the simulation box is not orthogonal.
    /// - `AtomError::InvalidPosition` if any of the atoms has no position.
    ///
    /// ## Example
    /// Calculate distance between atoms with index 3 and 17 in the xy plane.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let system = System::from_file("system.gro").unwrap();
    ///
    /// let result = match system.atoms_distance(3, 17, Dimension::XY) {
    ///     Ok(x) => x,
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;    
    ///     }
    /// };
    ///
    /// // print the result
    /// println!("XY-distance between atoms 3 and 17 is {}.", result);
    /// ```
    pub fn atoms_distance(
        &self,
        index1: usize,
        index2: usize,
        dim: Dimension,
    ) -> Result<f32, AtomError> {
        let atom1 = self.get_atom_as_ref(index1)?;
        let atom2 = self.get_atom_as_ref(index2)?;

        let simbox = simbox_check(self.get_box_as_ref()).map_err(AtomError::InvalidSimBox)?;

        atom1.distance(atom2, dim, simbox)
    }
}

/// Calculate contribution of an atom to the center of mass.
#[inline(always)]
fn center_atom_contribution(
    mut position: Vector3D,
    scaling: &Vector3D,
    simbox: &SimBox,
    mass: f32,
    sum_xi: &mut Vector3D,
    sum_zeta: &mut Vector3D,
) {
    // wrap position into the box
    position.wrap(simbox);

    let theta = Vector3D::new(position.x * scaling.x, position.y * scaling.y, position.z * scaling.z);

    sum_xi.x += mass * theta.x.cos();
    sum_xi.y += mass * theta.y.cos();
    sum_xi.z += mass * theta.z.cos();

    sum_zeta.x += mass * theta.x.sin();
    sum_zeta.y += mass * theta.y.sin();
    sum_zeta.z += mass * theta.z.sin();
}

/// Convert coordinates from their representation on a circle to their representation on a line.
#[inline(always)]
fn from_circle_to_line(zeta: Vector3D, xi: Vector3D, scaling: &Vector3D) -> Vector3D {
    let theta = Vector3D::new(
        (-zeta.x).atan2(-xi.x) + consts::PI,
        (-zeta.y).atan2(-xi.y) + consts::PI,
        (-zeta.z).atan2(-xi.z) + consts::PI,
    );

    Vector3D::new(
        theta.x / scaling.x,
        theta.y / scaling.y,
        theta.z / scaling.z,
    )
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::assert_approx_eq;

    use crate::errors::SimBoxError;
    use crate::structures::atom::Atom;
    use crate::structures::element::Elements;

    #[test]
    fn center_single_atom() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([4.5, 3.2, 1.7].into());

        let atoms = vec![atom1];
        let system = System::new("Artificial system.", atoms, Some([10.0, 10.0, 10.0].into()));

        let center = system.group_get_center("all").unwrap();

        assert_approx_eq!(f32, center.x, 4.5);
        assert_approx_eq!(f32, center.y, 3.2);
        assert_approx_eq!(f32, center.z, 1.7);
    }

    #[test]
    fn center_two_atoms() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([4.5, 3.2, 1.7].into());

        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([4.0, 2.8, 3.0].into());

        let atoms = vec![atom1, atom2];
        let system = System::new("Artificial system.", atoms, Some([10.0, 10.0, 10.0].into()));

        let center = system.group_get_center("all").unwrap();

        assert_approx_eq!(f32, center.x, 4.25);
        assert_approx_eq!(f32, center.y, 3.0);
        assert_approx_eq!(f32, center.z, 2.35);
    }

    #[test]
    fn center_two_atoms_pbc() {
        let atom1 = Atom::new(1, "LYS", 1, "BB").with_position([4.5, 3.2, 1.7].into());

        let atom2 = Atom::new(1, "LYS", 2, "SC1").with_position([9.8, 9.5, 3.0].into());

        let atoms = vec![atom1, atom2];
        let system = System::new("Artificial system.", atoms, Some([10.0, 10.0, 10.0].into()));

        let center = system.group_get_center("all").unwrap();

        assert_approx_eq!(f32, center.x, 2.15);
        assert_approx_eq!(f32, center.y, 1.35);
        assert_approx_eq!(f32, center.z, 2.35);
    }

    #[test]
    fn center_several_atoms_pbc() {
        let atom_positions: Vec<[f32; 3]> = vec![
            [3.3, 0.3, 2.5],
            [4.3, 1.2, 9.8],
            [3.2, 5.6, 0.5],
            [0.2, 9.0, 6.6],
            [8.7, 5.0, 2.4],
        ];
        let mut atoms = Vec::new();
        for (i, position) in atom_positions.into_iter().enumerate() {
            let atom = Atom::new(i, "UNK", i, "BB").with_position(position.into());

            atoms.push(atom);
        }

        let system = System::new("Artificial system.", atoms, Some([10.0, 10.0, 10.0].into()));

        let center = system.group_get_center("all").unwrap();

        assert_approx_eq!(f32, center.x, 2.634386, epsilon = 0.0001);
        assert_approx_eq!(f32, center.y, 9.775156, epsilon = 0.0001);
        assert_approx_eq!(f32, center.z, 1.174800, epsilon = 0.0001);
    }

    #[test]
    fn center_several_atoms_outofbox() {
        let atom_positions: Vec<[f32; 3]> = vec![
            [3.3, 10.3, 2.5],
            [4.3, 1.2, -0.2],
            [13.2, 15.6, 0.5],
            [10.2, -1.0, 6.6],
            [-1.3, 5.0, 2.4],
        ];
        let mut atoms = Vec::new();
        for (i, position) in atom_positions.into_iter().enumerate() {
            let atom = Atom::new(i, "UNK", i, "BB").with_position(position.into());

            atoms.push(atom);
        }

        let system = System::new("Artificial system.", atoms, Some([10.0, 10.0, 10.0].into()));

        let center = system.group_get_center("all").unwrap();

        assert_approx_eq!(f32, center.x, 2.634386, epsilon = 0.0001);
        assert_approx_eq!(f32, center.y, 9.775156, epsilon = 0.0001);
        assert_approx_eq!(f32, center.z, 1.174800, epsilon = 0.0001);
    }

    #[test]
    fn center_real_system() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let center_mem = system.group_get_center("Membrane").unwrap();
        let center_prot = system.group_get_center("Protein").unwrap();

        assert_approx_eq!(f32, center_mem.x, 3.575004, epsilon = 0.0001);
        assert_approx_eq!(f32, center_mem.y, 8.009330, epsilon = 0.0001);
        assert_approx_eq!(f32, center_mem.z, 5.779888, epsilon = 0.0001);

        assert_approx_eq!(f32, center_prot.x, 9.857101, epsilon = 0.0001);
        assert_approx_eq!(f32, center_prot.y, 2.462601, epsilon = 0.0001);
        assert_approx_eq!(f32, center_prot.z, 5.461296, epsilon = 0.0001);
    }

    #[test]
    fn center_real_system_fail_invalid_group() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        match system.group_get_center("Nonexistent") {
            Err(GroupError::NotFound(e)) => assert_eq!(e, "Nonexistent"),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn center_real_system_fail_invalid_simbox() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();
        system.reset_box();

        match system.group_get_center("Protein") {
            Err(GroupError::InvalidSimBox(SimBoxError::DoesNotExist)) => (),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn center_real_system_fail_invalid_position() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.get_atom_as_ref_mut(15).unwrap().reset_position();

        match system.group_get_center("Protein") {
            Err(GroupError::InvalidPosition(PositionError::NoPosition(x))) => assert_eq!(x, 16),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn center_real_system_fail_empty_group() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.group_create("Empty", "resname NON").unwrap();

        match system.group_get_center("Empty") {
            Err(GroupError::EmptyGroup(x)) => assert_eq!(x, "Empty"),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn com_single_atom() {
        let atom1 = Atom::new(1, "LYS", 1, "BB")
            .with_position([4.5, 3.2, 1.7].into())
            .with_mass(12.8);

        let atoms = vec![atom1];
        let system = System::new("Artificial system.", atoms, Some([10.0, 10.0, 10.0].into()));

        let center = system.group_get_com("all").unwrap();

        assert_approx_eq!(f32, center.x, 4.5);
        assert_approx_eq!(f32, center.y, 3.2);
        assert_approx_eq!(f32, center.z, 1.7);
    }

    #[test]
    fn com_two_atoms() {
        let atom1 = Atom::new(1, "LYS", 1, "BB")
            .with_position([4.5, 3.2, 1.7].into())
            .with_mass(12.8);

        let atom2 = Atom::new(1, "LYS", 2, "SC1")
            .with_position([4.0, 2.8, 3.0].into())
            .with_mass(0.4);

        let atoms = vec![atom1, atom2];
        let system = System::new("Artificial system.", atoms, Some([10.0, 10.0, 10.0].into()));

        let center = system.group_get_com("all").unwrap();

        assert_approx_eq!(f32, center.x, 4.485, epsilon = 0.0001);
        assert_approx_eq!(f32, center.y, 3.188, epsilon = 0.0001);
        assert_approx_eq!(f32, center.z, 1.73549, epsilon = 0.0001);
    }

    #[test]
    fn com_two_atoms_pbc() {
        let atom1 = Atom::new(1, "LYS", 1, "BB")
            .with_position([4.5, 3.2, 1.7].into())
            .with_mass(12.8);

        let atom2 = Atom::new(1, "LYS", 2, "SC1")
            .with_position([9.8, 9.5, 3.0].into())
            .with_mass(0.4);

        let atoms = vec![atom1, atom2];
        let system = System::new("Artificial system.", atoms, Some([10.0, 10.0, 10.0].into()));

        let center = system.group_get_com("all").unwrap();

        assert_approx_eq!(f32, center.x, 4.4904, epsilon = 0.0001);
        assert_approx_eq!(f32, center.y, 3.1630, epsilon = 0.0001);
        assert_approx_eq!(f32, center.z, 1.7355, epsilon = 0.0001);
    }

    #[test]
    fn com_several_atoms_pbc() {
        let atom_positions: Vec<[f32; 3]> = vec![
            [3.3, 0.3, 2.5],
            [4.3, 1.2, 9.8],
            [3.2, 5.6, 0.5],
            [0.2, 9.0, 6.6],
            [8.7, 5.0, 2.4],
        ];

        let masses = vec![10.3, 5.4, 3.8, 10.1, 7.6];

        let mut atoms = Vec::new();
        for (i, position) in atom_positions.into_iter().enumerate() {
            let atom = Atom::new(i, "UNK", i, "BB")
                .with_position(position.into())
                .with_mass(masses[i]);

            atoms.push(atom);
        }

        let system = System::new("Artificial system.", atoms, Some([10.0, 10.0, 10.0].into()));

        let center = system.group_get_com("all").unwrap();

        assert_approx_eq!(f32, center.x, 1.9526, epsilon = 0.0001);
        assert_approx_eq!(f32, center.y, 9.7567, epsilon = 0.0001);
        assert_approx_eq!(f32, center.z, 1.8812, epsilon = 0.0001);
    }

    #[test]
    fn com_several_atoms_outofbox() {
        let atom_positions: Vec<[f32; 3]> = vec![
            [3.3, 10.3, 2.5],
            [4.3, 1.2, -0.2],
            [13.2, 15.6, 0.5],
            [10.2, -1.0, 6.6],
            [-1.3, 5.0, 2.4],
        ];

        let masses = vec![10.3, 5.4, 3.8, 10.1, 7.6];

        let mut atoms = Vec::new();
        for (i, position) in atom_positions.into_iter().enumerate() {
            let atom = Atom::new(i, "UNK", i, "BB")
                .with_position(position.into())
                .with_mass(masses[i]);

            atoms.push(atom);
        }

        let system = System::new("Artificial system.", atoms, Some([10.0, 10.0, 10.0].into()));

        let center = system.group_get_com("all").unwrap();

        assert_approx_eq!(f32, center.x, 1.9526, epsilon = 0.0001);
        assert_approx_eq!(f32, center.y, 9.7567, epsilon = 0.0001);
        assert_approx_eq!(f32, center.z, 1.8812, epsilon = 0.0001);
    }

    #[test]
    fn com_real_system_same_mass() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let center_mem = system.group_get_center("Membrane").unwrap();
        let center_prot = system.group_get_center("Protein").unwrap();

        for atom in system.atoms_iter_mut() {
            atom.set_mass(12.3);
        }

        let com_mem = system.group_get_com("Membrane").unwrap();
        let com_prot = system.group_get_com("Protein").unwrap();

        assert_approx_eq!(f32, center_prot.x, com_prot.x, epsilon = 0.0001);
        assert_approx_eq!(f32, center_prot.y, com_prot.y, epsilon = 0.0001);
        assert_approx_eq!(f32, center_prot.z, com_prot.z, epsilon = 0.0001);

        assert_approx_eq!(f32, center_mem.x, com_mem.x, epsilon = 0.0001);
        assert_approx_eq!(f32, center_mem.y, com_mem.y, epsilon = 0.0001);
        assert_approx_eq!(f32, center_mem.z, com_mem.z, epsilon = 0.0001);
    }

    #[test]
    fn com_real_system() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();

        system.group_create("Peptide", "@protein").unwrap();
        system.group_create("Membrane", "@membrane").unwrap();

        system.guess_elements(Elements::default()).unwrap();

        let com_prot = system.group_get_com("Peptide").unwrap();
        let com_mem = system.group_get_com("Membrane").unwrap();

        assert_approx_eq!(f32, com_prot.x, 4.047723, epsilon = 0.0001);
        assert_approx_eq!(f32, com_prot.y, 3.764632, epsilon = 0.0001);
        assert_approx_eq!(f32, com_prot.z, 3.2633042, epsilon = 0.0001);

        assert_approx_eq!(f32, com_mem.x, 1.44719, epsilon = 0.0001);
        assert_approx_eq!(f32, com_mem.y, 0.45375, epsilon = 0.0001);
        assert_approx_eq!(f32, com_mem.z, 3.74161, epsilon = 0.0001);
    }

    #[test]
    fn com_real_system_fail_invalid_group() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        match system.group_get_com("Nonexistent") {
            Err(GroupError::NotFound(e)) => assert_eq!(e, "Nonexistent"),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn com_real_system_fail_invalid_simbox() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();
        system.reset_box();

        match system.group_get_com("Protein") {
            Err(GroupError::InvalidSimBox(SimBoxError::DoesNotExist)) => (),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn com_real_system_fail_invalid_position() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        for atom in system.atoms_iter_mut() {
            atom.set_mass(10.3);
        }

        system.get_atom_as_ref_mut(15).unwrap().reset_position();

        match system.group_get_com("Protein") {
            Err(GroupError::InvalidPosition(PositionError::NoPosition(x))) => assert_eq!(x, 16),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn com_real_system_fail_invalid_mass() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        match system.group_get_com("Protein") {
            Err(GroupError::InvalidMass(MassError::NoMass(x))) => assert_eq!(x, 1),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn com_real_system_fail_empty_group() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.group_create("Empty", "resname NON").unwrap();

        match system.group_get_com("Empty") {
            Err(GroupError::EmptyGroup(x)) => assert_eq!(x, "Empty"),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn group_distance_x() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::X)
            .unwrap();
        assert_approx_eq!(f32, dist, 6.282097, epsilon = 0.0001);
    }

    #[test]
    fn group_distance_y() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::Y)
            .unwrap();
        assert_approx_eq!(f32, dist, -5.546729, epsilon = 0.0001);
    }

    #[test]
    fn group_distance_z() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::Z)
            .unwrap();
        assert_approx_eq!(f32, dist, -0.31859207, epsilon = 0.0001);
    }

    #[test]
    fn group_distance_xy() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::XY)
            .unwrap();
        assert_approx_eq!(f32, dist, 8.38039, epsilon = 0.0001);
    }

    #[test]
    fn group_distance_xz() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::XZ)
            .unwrap();
        assert_approx_eq!(f32, dist, 6.29017, epsilon = 0.0001);
    }

    #[test]
    fn group_distance_yz() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::YZ)
            .unwrap();
        assert_approx_eq!(f32, dist, 5.555871, epsilon = 0.0001);
    }

    #[test]
    fn group_distance_xyz() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::XYZ)
            .unwrap();
        assert_approx_eq!(f32, dist, 8.386444, epsilon = 0.0001);
    }

    #[test]
    fn group_distance_none() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::None)
            .unwrap();
        assert_approx_eq!(f32, dist, 0.0);
    }

    #[test]
    fn group_distance_fail_1() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        match system.group_distance("PRotein", "Membrane", Dimension::XYZ) {
            Err(GroupError::NotFound(x)) => assert_eq!(x, "PRotein"),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn group_distance_fail_2() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        match system.group_distance("Protein", "Nonexistent", Dimension::XYZ) {
            Err(GroupError::NotFound(x)) => assert_eq!(x, "Nonexistent"),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn group_distance_fail_simbox() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();
        system.reset_box();

        match system.group_distance("Protein", "Membrane", Dimension::XYZ) {
            Err(GroupError::InvalidSimBox(SimBoxError::DoesNotExist)) => (),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn group_distance_fail_position() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.get_atom_as_ref_mut(15).unwrap().reset_position();

        match system.group_distance("Protein", "Membrane", Dimension::XYZ) {
            Err(GroupError::InvalidPosition(PositionError::NoPosition(x))) => assert_eq!(x, 16),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn group_all_distances_xyz() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let n_atoms = system.group_get_n_atoms("Protein").unwrap() as usize;
        let distances = system
            .group_all_distances("Protein", "Protein", Dimension::XYZ)
            .unwrap();

        assert_eq!(distances.len(), n_atoms);
        assert_eq!(distances[0].len(), n_atoms);

        for i in 0..n_atoms {
            for j in 0..n_atoms {
                assert_approx_eq!(f32, distances[i][j], distances[j][i]);

                if i == j {
                    assert_eq!(distances[i][j], 0.0);
                }
            }
        }

        // get maximal value
        let max = distances
            .iter()
            .flatten()
            .fold(std::f32::NEG_INFINITY, |max, &current| max.max(current));
        assert_approx_eq!(f32, max, 4.597961);

        assert_approx_eq!(f32, distances[0][1], 0.31040135);
        assert_approx_eq!(f32, distances[n_atoms - 1][0], 4.266728);
        assert_approx_eq!(f32, distances[n_atoms - 1][n_atoms - 2], 0.31425142);
    }

    #[test]
    fn group_all_distances_z() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let n_atoms = system.group_get_n_atoms("Protein").unwrap() as usize;
        let distances = system
            .group_all_distances("Protein", "Protein", Dimension::Z)
            .unwrap();

        assert_eq!(distances.len(), n_atoms);
        assert_eq!(distances[0].len(), n_atoms);

        for i in 0..n_atoms {
            for j in 0..n_atoms {
                assert_approx_eq!(f32, distances[i][j], -distances[j][i]);

                if i == j {
                    assert_eq!(distances[i][j], 0.0);
                }
            }
        }

        // get maximal value
        let max = distances
            .iter()
            .flatten()
            .fold(std::f32::NEG_INFINITY, |max, &current| max.max(current));
        assert_approx_eq!(f32, max, 4.383, epsilon = 0.00001);

        // get the minimal value
        let min = distances
            .iter()
            .flatten()
            .fold(std::f32::INFINITY, |min, &current| min.min(current));
        assert_approx_eq!(f32, min, -4.383, epsilon = 0.00001);

        assert_approx_eq!(f32, distances[0][1], 0.0900, epsilon = 0.00001);
        assert_approx_eq!(f32, distances[n_atoms - 1][0], -4.213, epsilon = 0.00001);
        assert_approx_eq!(
            f32,
            distances[n_atoms - 1][n_atoms - 2],
            -0.101,
            epsilon = 0.00001
        );
    }

    #[test]
    fn group_all_distances_xy() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let n_atoms_membrane = system.group_get_n_atoms("Membrane").unwrap() as usize;
        let n_atoms_protein = system.group_get_n_atoms("Protein").unwrap() as usize;
        let distances = system
            .group_all_distances("Membrane", "Protein", Dimension::XY)
            .unwrap();

        assert_eq!(distances.len(), n_atoms_membrane);
        assert_eq!(distances[0].len(), n_atoms_protein);

        // get maximal value
        let max = distances
            .iter()
            .flatten()
            .fold(std::f32::NEG_INFINITY, |max, &current| max.max(current));
        assert_approx_eq!(f32, max, 9.190487, epsilon = 0.00001);

        // get the minimal value
        let min = distances
            .iter()
            .flatten()
            .fold(std::f32::INFINITY, |min, &current| min.min(current));
        assert_approx_eq!(f32, min, 0.02607, epsilon = 0.00001);

        assert_approx_eq!(f32, distances[0][0], 3.747651);
        assert_approx_eq!(f32, distances[1240][12], 3.7207017);
        assert_approx_eq!(f32, distances[12][34], 6.2494035);
        assert_approx_eq!(f32, distances[6143][60], 4.7850933);
    }

    #[test]
    fn group_all_distances_fail_1() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        match system.group_all_distances("Nonexistent", "Protein", Dimension::XYZ) {
            Err(GroupError::NotFound(x)) => assert_eq!(x, "Nonexistent"),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn group_all_distances_fail_2() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        match system.group_all_distances("Membrane", "Nonexistent", Dimension::XYZ) {
            Err(GroupError::NotFound(x)) => assert_eq!(x, "Nonexistent"),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn group_all_distances_fail_simbox() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();
        system.reset_box();

        match system.group_all_distances("Membrane", "Protein", Dimension::XYZ) {
            Err(GroupError::InvalidSimBox(SimBoxError::DoesNotExist)) => (),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn group_all_distances_fail_position() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();
        system.get_atom_as_ref_mut(15).unwrap().reset_position();

        match system.group_all_distances("Membrane", "Protein", Dimension::XYZ) {
            Err(GroupError::InvalidPosition(PositionError::NoPosition(x))) => assert_eq!(x, 16),
            Ok(_) => panic!("Calculating center should have failed, but it was successful."),
            Err(e) => panic!(
                "Failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn atoms_distance_xyz() {
        let system = System::from_file("test_files/example.gro").unwrap();

        let n_atoms = system.get_n_atoms();

        assert_approx_eq!(
            f32,
            system.atoms_distance(0, 1, Dimension::XYZ).unwrap(),
            0.31040135
        );
        assert_approx_eq!(
            f32,
            system
                .atoms_distance(n_atoms - 1, 0, Dimension::XYZ)
                .unwrap(),
            6.664787
        );
        assert_approx_eq!(
            f32,
            system
                .atoms_distance(n_atoms - 1, n_atoms - 2, Dimension::XYZ)
                .unwrap(),
            4.062491
        );
    }

    #[test]
    fn atoms_distance_fail_1() {
        let system = System::from_file("test_files/example.gro").unwrap();

        match system.atoms_distance(12, 16844, Dimension::XY) {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(AtomError::OutOfRange(e)) => assert_eq!(e, 16844),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }

        match system.atoms_distance(197_392, 12, Dimension::YZ) {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(AtomError::OutOfRange(e)) => assert_eq!(e, 197_392),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn atoms_distance_fail_simbox() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.reset_box();

        match system.atoms_distance(12, 15, Dimension::XYZ) {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(AtomError::InvalidSimBox(SimBoxError::DoesNotExist)) => (),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn atoms_distance_fail_position() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.get_atom_as_ref_mut(15).unwrap().reset_position();

        match system.atoms_distance(12, 15, Dimension::XYZ) {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(AtomError::InvalidPosition(PositionError::NoPosition(x))) => assert_eq!(x, 16),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }
}
