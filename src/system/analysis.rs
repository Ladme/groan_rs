// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of various methods for analysis of `System`.

use crate::errors::{AtomError, GroupError, PositionError};
use crate::structures::simbox::simbox_check;
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
        let simbox =
            simbox_check(self.get_box_as_ref()).map_err(|x| GroupError::InvalidSimBox(x))?;

        let reciprocal_box =
            Vector3D::from([1.0f32 / simbox.x, 1.0f32 / simbox.y, 1.0f32 / simbox.z]);

        let mut sum_xi = Vector3D::default();
        let mut sum_zeta = Vector3D::default();

        for atom in self.group_iter(name)? {
            // make sure that each coordinate is inside the box
            let mut coordinates = match atom.get_position() {
                Some(x) => x.clone(),
                None => return Err(GroupError::InvalidPosition(PositionError::NoPosition(atom.get_atom_number())))
            };
            coordinates.wrap(simbox);

            // calculate magic angles
            let theta = Vector3D {
                x: coordinates.x * reciprocal_box.x * PI_X2,
                y: coordinates.y * reciprocal_box.y * PI_X2,
                z: coordinates.z * reciprocal_box.z * PI_X2,
            };

            sum_xi.x += theta.x.cos();
            sum_xi.y += theta.y.cos();
            sum_xi.z += theta.z.cos();

            sum_zeta.x += theta.x.sin();
            sum_zeta.y += theta.y.sin();
            sum_zeta.z += theta.z.sin();
        }

        // transform magic angles into real coordinates
        let final_theta = Vector3D {
            x: (-sum_zeta.x).atan2(-sum_xi.x) + consts::PI,
            y: (-sum_zeta.y).atan2(-sum_xi.y) + consts::PI,
            z: (-sum_zeta.z).atan2(-sum_xi.z) + consts::PI,
        };

        let center = Vector3D {
            x: simbox.x * (final_theta.x / PI_X2),
            y: simbox.y * (final_theta.y / PI_X2),
            z: simbox.z * (final_theta.z / PI_X2),
        };

        Ok(center)
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

        let simbox =
            simbox_check(self.get_box_as_ref()).map_err(|x| GroupError::InvalidSimBox(x))?;

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

        let simbox =
            simbox_check(self.get_box_as_ref()).map_err(|x| GroupError::InvalidSimBox(x))?;

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

        let simbox =
            simbox_check(self.get_box_as_ref()).map_err(|x| AtomError::InvalidSimBox(x))?;

        atom1.distance(atom2, dim, simbox)
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::assert_approx_eq;

    use crate::structures::atom::Atom;

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

        assert_approx_eq!(f32, center.x, 2.634386);
        assert_approx_eq!(f32, center.y, 9.775156);
        assert_approx_eq!(f32, center.z, 1.174800);
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

        assert_approx_eq!(f32, center.x, 2.634386);
        assert_approx_eq!(f32, center.y, 9.775156);
        assert_approx_eq!(f32, center.z, 1.174800);
    }

    #[test]
    fn center_real_system() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let center_mem = system.group_get_center("Membrane").unwrap();
        let center_prot = system.group_get_center("Protein").unwrap();

        assert_approx_eq!(f32, center_mem.x, 3.575004);
        assert_approx_eq!(f32, center_mem.y, 8.009330);
        assert_approx_eq!(f32, center_mem.z, 5.779888);

        assert_approx_eq!(f32, center_prot.x, 9.857101);
        assert_approx_eq!(f32, center_prot.y, 2.462601);
        assert_approx_eq!(f32, center_prot.z, 5.461296);
    }

    #[test]
    fn center_real_system_fail() {
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
    fn group_distance_x() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::X)
            .unwrap();
        assert_approx_eq!(f32, dist, 6.282097);
    }

    #[test]
    fn group_distance_y() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::Y)
            .unwrap();
        assert_approx_eq!(f32, dist, -5.546729);
    }

    #[test]
    fn group_distance_z() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::Z)
            .unwrap();
        assert_approx_eq!(f32, dist, -0.31859207);
    }

    #[test]
    fn group_distance_xy() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::XY)
            .unwrap();
        assert_approx_eq!(f32, dist, 8.38039);
    }

    #[test]
    fn group_distance_xz() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::XZ)
            .unwrap();
        assert_approx_eq!(f32, dist, 6.29017);
    }

    #[test]
    fn group_distance_yz() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::YZ)
            .unwrap();
        assert_approx_eq!(f32, dist, 5.555871);
    }

    #[test]
    fn group_distance_xyz() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let dist = system
            .group_distance("Protein", "Membrane", Dimension::XYZ)
            .unwrap();
        assert_approx_eq!(f32, dist, 8.386444);
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

        if system
            .group_distance("PRotein", "Membrane", Dimension::XYZ)
            .is_ok()
        {
            panic!("`group_distance` should have failed but it was successful.");
        }
    }

    #[test]
    fn group_distance_fail_2() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        if system
            .group_distance("Protein", "Nonexistent", Dimension::XYZ)
            .is_ok()
        {
            panic!("`group_distance` should have failed but it was successful.");
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
    fn atoms_distance_fail() {
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
    fn group_all_distances_fail_1() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        if system
            .group_all_distances("Nonexistent", "Protein", Dimension::XYZ)
            .is_ok()
        {
            panic!("`group_all_distances` should have failed but it was successful.");
        }
    }

    #[test]
    fn group_all_distances_fail_2() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        if system
            .group_all_distances("Membrane", "Nonexistent", Dimension::XYZ)
            .is_ok()
        {
            panic!("`group_all_distances` should have failed but it was successful.");
        }
    }
}
