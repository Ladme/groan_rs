// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of methods for the calculation of RMSD.

use std::ops::Deref;

use crate::{
    errors::{GroupError, RMSDError},
    prelude::SimBox,
    structures::{atom::Atom, vector3d::Vector3D},
    system::System,
};
use nalgebra::Matrix3;

impl System {
    /// Calculate the Root Mean Square Deviation (RMSD) between the specified group of atoms
    /// in the current system and a reference system.
    ///
    /// Uses mass-weighted Kabsch algorithm.
    ///
    /// ## Returns
    /// - If successful, returns the minimum RMSD (in nanometers) between the atoms
    ///   of the specified group in `reference` and `self`.
    /// - Returns an `RMSDError` if the calculation fails.
    ///
    /// ## Example
    /// Calculating RMSD for every 10th frame of an XTC trajectory.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// fn calculate_rmsd() -> Result<Vec<f32>, Box<dyn std::error::Error + Send + Sync>> {
    ///     let mut system = System::from_file("system.gro")?;
    ///
    ///     // create a group `Protein` for which the RMSD will be calculated
    ///     system.group_create("Protein", "@protein")?;
    ///
    ///     // use the input structure file as a reference structure
    ///     let reference = system.clone();
    ///
    ///     // iterate over the xtc trajectory and calculate RMSD for every 10th frame
    ///     let mut rmsd = Vec::new();
    ///     for frame in system.xtc_iter("trajectory.xtc")?.with_step(10)? {
    ///         let frame = frame?;
    ///         rmsd.push(frame.calc_rmsd(&reference, "Protein")?);
    ///     }
    ///
    ///     // return the collected data
    ///     Ok(rmsd)
    /// }
    ///
    /// ```
    ///
    /// ## Notes
    /// - This method requires both systems to have a valid orthogonal simulation box;
    ///   otherwise, an error will be returned.
    /// - The method performs a rigid-body alignment of the atoms in the specified group using
    ///   the Kabsch algorithm before calculating the RMSD.
    /// - Mass weighting **is** performed during the alignment.
    /// - The RMSD is calculated in nanometers (nm).
    /// - Neither the current system (`self`) nor the reference system is modified by this method.
    /// - The group does not have to be centered or fully contained inside the simulation box.
    ///   The method itself performs the centering.
    #[inline]
    pub fn calc_rmsd(&self, reference: &System, group: &str) -> Result<f32, RMSDError> {
        let (_, _, rmsd) = self.calc_rmsd_rot_trans(reference, group)?;
        Ok(rmsd)
    }

    #[inline]
    pub fn calc_rmsd_and_fit(&mut self, reference: &System, group: &str) -> Result<f32, RMSDError> {
        let (rot, trans, rmsd) = self.calc_rmsd_rot_trans(reference, group)?;

        let rot_transposed = rot.transpose();
        //let simbox = self.get_box_copy().unwrap();
        let box_center = self.get_box_center().unwrap();

        self.atoms_iter_mut().for_each(|atom| {
            atom.translate_nopbc(&trans).unwrap();
            atom.rotate_nopbc(&rot_transposed).unwrap();
            atom.translate_nopbc(&box_center).unwrap();
        });

        Ok(rmsd)
    }

    /// Calculate RMSD, rotation matrix and translation vector.
    fn calc_rmsd_rot_trans(
        &self,
        reference: &System,
        group: &str,
    ) -> Result<(Matrix3<f32>, Vector3D, f32), RMSDError> {
        // check that the group exists and has a consistent number of atoms
        let n_atoms_target = extract_n_atoms(group, self.group_get_n_atoms(group))?;
        let n_atoms_reference = extract_n_atoms(group, reference.group_get_n_atoms(group))?;

        if n_atoms_target != n_atoms_reference {
            return Err(RMSDError::InconsistentGroup(
                group.to_owned(),
                n_atoms_reference,
                n_atoms_target,
            ));
        }

        if n_atoms_target == 0 {
            return Err(RMSDError::EmptyGroup(group.to_owned()));
        }

        // check that the simulation box is defined and orthogonal & extract box center
        let reference_box_center = reference
            .get_box_center()
            .map_err(RMSDError::InvalidSimBox)?;
        let target_box_center = self.get_box_center().map_err(RMSDError::InvalidSimBox)?;

        // get the center of mass of the reference
        let reference_center = get_com(reference, group)?;
        let target_center = get_com(self, group)?;

        // extract the coordinates of the atoms of the group and shift them
        let reference_shift = Vector3D(reference_box_center.deref() - reference_center.deref());
        let target_shift = Vector3D(target_box_center.deref() - target_center.deref());

        let reference_coordinates = shift_and_wrap_coordinates(
            extract_coordinates(reference.group_iter(group)),
            &reference_shift,
            reference.get_box().expect(
                "FATAL GROAN ERROR | System::calc_rmsd_rot_trans | Reference SimBox should exist.",
            ),
        );

        let target_coordinates = shift_and_wrap_coordinates(
            extract_coordinates(self.group_iter(group)),
            &target_shift,
            self.get_box().expect(
                "FATAL GROAN ERROR | System::calc_rmsd_rot_trans | Current SimBox should exist.",
            ),
        );

        // exctract masses
        let masses = self
            .group_iter(group)
            .unwrap()
            .map(|atom| {
                atom.get_mass().expect(
                    "FATAL GROAN ERROR | System::calc_rmsd_rot_trans | Mass should be defined.",
                )
            })
            .collect::<Vec<f32>>();

        let sum_masses = masses.iter().sum::<f32>();

        // calculate RMSD
        Ok(kabsch_rmsd(
            &target_coordinates,
            &reference_coordinates,
            &masses,
            &target_box_center,
            &reference_box_center,
            sum_masses,
        ))
    }
}

/// Auxiliary function for extracting the number of atoms in a group from system.
#[inline(always)]
fn extract_n_atoms(group: &str, result: Result<usize, GroupError>) -> Result<usize, RMSDError> {
    match result {
            Ok(0) => Err(RMSDError::EmptyGroup(group.to_owned())),
            Ok(x) => Ok(x),
            Err(GroupError::NotFound(_)) => Err(RMSDError::NonexistentGroup(group.to_owned())),
            Err(_) => panic!("FATAL GROAN ERROR | rmsd::extract_n_atoms | System::group_get_n_atoms returned an unexpected error."),
        }
}

/// Auxiliary function for extracting the coordinates of relevant atoms from the system.
#[inline(always)]
fn extract_coordinates<'a, I>(iter: Result<I, GroupError>) -> Vec<Vector3D>
where
    I: Iterator<Item = &'a Atom>,
{
    iter.expect("FATAL GROAN ERROR | rmsd::extract_coordinates | Group should exist.")
        .map(|x| {
            x.get_position()
                .cloned()
                .expect("FATAL GROAN ERROR | rmsd::extract_coordinates | All atoms should have defined positions.")
        })
        .collect::<Vec<Vector3D>>()
}

/// Auxiliary function for shifting and wrapping coordinates.
#[inline(always)]
fn shift_and_wrap_coordinates(
    coordinates: Vec<Vector3D>,
    shift: &Vector3D,
    simbox: &SimBox,
) -> Vec<Vector3D> {
    coordinates
        .iter()
        .map(|pos| {
            let mut new = Vector3D(pos.deref() + shift.deref());
            new.wrap(simbox);
            new
        })
        .collect::<Vec<Vector3D>>()
}

/// Auxiliary function for obtaining the center of mass of the specified group.
#[inline(always)]
fn get_com(system: &System, group: &str) -> Result<Vector3D, RMSDError> {
    match system.group_get_com(group) {
        Ok(x) => Ok(x),
        Err(GroupError::InvalidPosition(x)) => Err(RMSDError::InvalidPosition(x)),
        Err(GroupError::InvalidMass(x)) => Err(RMSDError::InvalidMass(x)),
        Err(_) => panic!("FATAL GROAN ERROR | rmsd::get_center | Unexpected error type returned from System::group_get_center."),
    }
}

/// Calculate the optimal rotation matrix, translation vector, and RMSD
/// to align two sets of points (P -> Q) using the Kabsch algorithm.
///
/// ## Parameters
/// - `p`: target vector of points
/// - `q`: reference vector of points
/// - `w`: weight (mass) of the individual points
/// - `centroid_p`: center of geometry (or mass) of points `p`
/// - `centroid_q`: center of geometry (or mass) of points `q`
/// - `sum_w`: sum of weights of the individual points
///
/// ## Returns
/// A tuple containing the optimal rotation matrix, the optimal translation vector, and the RMSD.
///
/// ## Panics
/// - Panics if the number of points `p` does not match the number of points `q`.
/// - Panics if the number of items in `w` does not match the number of points.
fn kabsch_rmsd(
    p: &[Vector3D],
    q: &[Vector3D],
    w: &[f32],
    centroid_p: &Vector3D,
    centroid_q: &Vector3D,
    sum_w: f32,
) -> (Matrix3<f32>, Vector3D, f32) {
    assert_eq!(
        p.len(),
        q.len(),
        "FATAL GROAN ERROR | rmsd::kabsch_rmsd | Number of points `p` and `q` does not match."
    );
    assert_eq!(p.len(), w.len(), "FATAL GROAN ERROR | rmsd::kabsch_rmsd | Number of points does not match the number of weights.");

    // center the points
    let p_centered: Vec<Vector3D> = p.iter().map(|point| point - centroid_p).collect();
    let q_centered: Vec<Vector3D> = q.iter().map(|point| point - centroid_q).collect();

    // compute the covariance matrix
    let mut h = Matrix3::zeros();
    for (p_c, q_c) in p_centered.iter().zip(q_centered.iter()) {
        h += p_c.deref() * q_c.transpose();
    }

    // perform Singular Value Decomposition (SVD)
    let svd = h.svd(true, true);
    let u = svd.u.unwrap();
    let v_t = svd.v_t.unwrap();

    let mut d = Matrix3::identity();
    if (u * v_t).determinant() < 0.0 {
        d[(2, 2)] = -1.0;
    }

    // calculate the rotation
    let r = u * d * v_t;

    // apply the rotation to the centered points
    let p_rotated: Vec<Vector3D> = p_centered
        .iter()
        .map(|point| Vector3D(r.transpose() * point.deref()))
        .collect();

    // calculate RMSD
    let rmsd = (p_rotated
        .iter()
        .zip(q_centered.iter())
        .zip(w.iter())
        .map(|((p_r, q_c), w)| w * (p_r.deref() - q_c.deref()).norm_squared())
        .sum::<f32>()
        / sum_w)
        .sqrt();

    // return the rotation matrix, translation vector, and RMSD
    (r, centroid_q - centroid_p, rmsd)
}

#[cfg(test)]
mod tests {
    use crate::errors::{MassError, PositionError, SimBoxError};

    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn test_kabsch_no_rotation_no_translation() {
        // no rotation, no translation, RMSD should be zero
        let p = [
            Vector3D::new(1.0, 0.0, 0.0),
            Vector3D::new(0.0, 1.0, 0.0),
            Vector3D::new(0.0, 0.0, 1.0),
        ];
        let q = [
            Vector3D::new(1.0, 0.0, 0.0),
            Vector3D::new(0.0, 1.0, 0.0),
            Vector3D::new(0.0, 0.0, 1.0),
        ];

        let masses = [1.0, 1.0, 1.0];

        let (rotation_matrix, translation_vector, rmsd) = kabsch_rmsd(
            &p,
            &q,
            &masses,
            &Vector3D::new(0.3333333, 0.3333333, 0.3333333),
            &Vector3D::new(0.3333333, 0.3333333, 0.3333333),
            3.0,
        );

        assert!(rotation_matrix.is_identity(1e-6));
        assert_approx_eq!(f32, translation_vector.norm(), 0.0, epsilon = 1e-6);
        assert_approx_eq!(f32, rmsd, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_kabsch_rotation_no_translation() {
        // rotation, no translation, RMSD should be zero
        let p = [
            Vector3D::new(1.0, 0.0, 0.0),
            Vector3D::new(0.0, 1.0, 0.0),
            Vector3D::new(0.0, 0.0, 1.0),
        ];

        let q = [
            Vector3D::new(0.6666667, 1.0, 0.0),
            Vector3D::new(-0.3333333, 0.0, 0.0),
            Vector3D::new(0.6666667, 0.0, 1.0),
        ];

        let masses = [1.0, 1.0, 1.0];

        let (rotation_matrix, translation_vector, rmsd) = kabsch_rmsd(
            &p,
            &q,
            &masses,
            &Vector3D::new(0.3333333, 0.3333333, 0.3333333),
            &Vector3D::new(0.3333333, 0.3333333, 0.3333333),
            3.0,
        );

        let expected_rotation_matrix =
            Matrix3::from([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        assert!(rotation_matrix.relative_eq(&expected_rotation_matrix, 1e-6, 1e-6));
        assert_approx_eq!(f32, translation_vector.norm(), 0.0, epsilon = 1e-6);
        assert_approx_eq!(f32, rmsd, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_kabsch_translation_no_rotation() {
        // translation, no rotation, RMSD should be zero
        let p = [
            Vector3D::new(1.0, 0.0, 0.0),
            Vector3D::new(0.0, 1.0, 0.0),
            Vector3D::new(0.0, 0.0, 1.0),
        ];

        let q = [
            Vector3D::new(2.0, 1.0, 1.0),
            Vector3D::new(1.0, 2.0, 1.0),
            Vector3D::new(1.0, 1.0, 2.0),
        ];

        let masses = [1.0, 1.0, 1.0];

        let (rotation_matrix, translation_vector, rmsd) = kabsch_rmsd(
            &p,
            &q,
            &masses,
            &Vector3D::new(0.3333333, 0.3333333, 0.3333333),
            &Vector3D::new(1.3333333, 1.3333333, 1.3333333),
            3.0,
        );

        assert!(rotation_matrix.is_identity(1e-6));
        let expected_translation = Vector3D::new(1.0, 1.0, 1.0);
        assert!(translation_vector.relative_eq(&expected_translation, 1e-6, 1e-6));
        assert_approx_eq!(f32, rmsd, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_kabsch_rotation_and_translation() {
        // translation, rotation, RMSD should be zero
        let p = [
            Vector3D::new(1.0, 0.0, 0.0),
            Vector3D::new(0.0, 1.0, 0.0),
            Vector3D::new(0.0, 0.0, 1.0),
        ];

        let q = [
            Vector3D::new(1.6666666, 2.0, 1.0),
            Vector3D::new(0.6666666, 1.0, 1.0),
            Vector3D::new(1.6666666, 1.0, 2.0),
        ];

        let masses = [1.0, 1.0, 1.0];

        let (rotation_matrix, translation_vector, rmsd) = kabsch_rmsd(
            &p,
            &q,
            &masses,
            &Vector3D::new(0.3333333, 0.3333333, 0.3333333),
            &Vector3D::new(1.3333333, 1.3333333, 1.3333333),
            3.0,
        );

        let expected_rotation_matrix =
            Matrix3::from([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        assert!(rotation_matrix.relative_eq(&expected_rotation_matrix, 1e-6, 1e-6));
        let expected_translation = Vector3D::new(1.0, 1.0, 1.0);
        assert!(translation_vector.relative_eq(&expected_translation, 1e-6, 1e-6));
        assert_approx_eq!(f32, rmsd, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn test_kabsch_nonzero_rmsd() {
        // translation, rotation, RMSD should be nonzero
        let p = [
            Vector3D::new(4.3, 2.1, -5.2),
            Vector3D::new(1.4, 2.1, 3.9),
            Vector3D::new(2.4, -3.3, 1.8),
        ];

        let q = [
            Vector3D::new(2.2, 0.0, 4.6),
            Vector3D::new(-1.4, 0.2, 0.3),
            Vector3D::new(1.3, 9.9, 11.3),
        ];

        let masses = [1.0, 1.0, 1.0];

        let center_p = Vector3D::new(2.7, 0.3, 0.16666667);
        let center_q = Vector3D::new(0.7, 3.3666667, 5.4);

        let (rotation_matrix, translation_vector, rmsd) =
            kabsch_rmsd(&p, &q, &masses, &center_p, &center_q, 3.0);

        let expected_rotation_matrix = Matrix3::from([
            [0.8842437, -0.10340805, -0.45543456],
            [0.2840647, -0.65496445, 0.70023507],
            [-0.37070346, -0.7485511, -0.5497733],
        ]);

        let expected_translation = Vector3D::new(-2.0, 3.066666, 5.233333);

        assert!(rotation_matrix.relative_eq(&expected_rotation_matrix, 1e-6, 1e-6));
        assert!(translation_vector.relative_eq(&expected_translation, 1e-6, 1e-6));
        assert_approx_eq!(f32, rmsd, 4.471225, epsilon = 1e-6);
    }

    #[test]
    fn test_calc_rmsd_same_structure() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        system.group_create("Protein", "@protein").unwrap();
        assert_approx_eq!(
            f32,
            system.calc_rmsd(&system, "Protein").unwrap(),
            0.0,
            epsilon = 1e-4
        )
    }

    #[test]
    fn test_calc_rmsd_trajectory() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();

        system.group_create("Protein", "@protein").unwrap();

        let mut reference = system.clone();

        // should work even if we remove position of some atom that is not in the group
        reference.get_atom_mut(176).unwrap().reset_position();

        let rmsd = system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .map(|frame| frame.unwrap().calc_rmsd(&reference, "Protein").unwrap())
            .collect::<Vec<f32>>();

        let expected = [
            0.23680355, 0.26356277, 0.26030675, 0.21396181, 0.22212411, 0.19429775, 0.26472768,
            0.27031693, 0.26426846, 0.23497732, 0.24261881,
        ];

        assert_eq!(rmsd.len(), expected.len());
        for (x, y) in rmsd.into_iter().zip(expected.into_iter()) {
            assert_approx_eq!(f32, x, y);
        }
    }

    #[test]
    fn test_calc_rmsd_fail_group_does_not_exist_in_reference() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        let reference = system.clone();

        system.group_create("Protein", "@protein").unwrap();

        match system.calc_rmsd(&reference, "Protein") {
            Ok(_) => panic!("Function should have failed."),
            Err(RMSDError::NonexistentGroup(x)) => assert_eq!(x, "Protein"),
            Err(_) => panic!("Function failed but incorrect error type returned."),
        }
    }

    #[test]
    fn test_calc_rmsd_fail_group_does_not_exist_in_self() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        system.group_create("Protein", "@protein").unwrap();
        let reference = system.clone();

        system.group_remove("Protein").unwrap();

        match system.calc_rmsd(&reference, "Protein") {
            Ok(_) => panic!("Function should have failed."),
            Err(RMSDError::NonexistentGroup(x)) => assert_eq!(x, "Protein"),
            Err(_) => panic!("Function failed but incorrect error type returned."),
        }
    }

    #[test]
    fn test_calc_rmsd_fail_inconsistent_group() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        system.group_create("Protein", "@protein").unwrap();
        let mut reference = system.clone();
        match reference.group_create("Protein", "serial 1 to 4") {
            Err(GroupError::AlreadyExistsWarning(_)) => (),
            _ => panic!("Group could not be created or was created without warning."),
        }

        match system.calc_rmsd(&reference, "Protein") {
            Ok(_) => panic!("Function should have failed."),
            Err(RMSDError::InconsistentGroup(x, i, j)) => {
                assert_eq!(x, "Protein");
                assert_eq!(i, 4);
                assert_eq!(j, 61);
            }
            Err(_) => panic!("Function failed but incorrect error type returned."),
        }
    }

    #[test]
    fn test_calc_rmsd_fail_empty_group() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        system.group_create("Protein", "not all").unwrap();

        let reference = system.clone();

        match system.calc_rmsd(&reference, "Protein") {
            Ok(_) => panic!("Function should have failed."),
            Err(RMSDError::EmptyGroup(x)) => assert_eq!(x, "Protein"),
            Err(_) => panic!("Function failed but incorrect error type returned."),
        }
    }

    #[test]
    fn test_calc_rmsd_fail_no_position() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let reference = system.clone();

        system.get_atom_mut(14).unwrap().reset_position();

        match system.calc_rmsd(&reference, "Protein") {
            Ok(_) => panic!("Function should have failed."),
            Err(RMSDError::InvalidPosition(PositionError::NoPosition(x))) => assert_eq!(x, 14),
            Err(_) => panic!("Function failed but incorrect error type returned."),
        }
    }

    #[test]
    fn test_calc_rmsd_fail_no_box() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let mut reference = system.clone();

        reference.reset_box();

        match system.calc_rmsd(&reference, "Protein") {
            Ok(_) => panic!("Function should have failed."),
            Err(RMSDError::InvalidSimBox(SimBoxError::DoesNotExist)) => (),
            Err(_) => panic!("Function failed but incorrect error type returned."),
        }
    }

    #[test]
    fn test_calc_rmsd_fail_no_mass() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let reference = system.clone();

        match system.calc_rmsd(&reference, "Protein") {
            Ok(_) => panic!("Function should have failed."),
            Err(RMSDError::InvalidMass(MassError::NoMass(_))) => (),
            Err(_) => panic!("Function failed but incorrect error type returned."),
        }
    }
}
