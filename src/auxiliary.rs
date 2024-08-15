// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Simple, auxiliary functions and constants used through the `groan_rs` library.

use std::{f32::consts, ops::Deref};

use crate::structures::{simbox::SimBox, vector3d::Vector3D};
use nalgebra::Matrix3;

/******************************/
/*         CONSTANTS          */
/******************************/

/// PI times 2.
pub(crate) const PI_X2: f32 = consts::PI * 2.0f32;

/// Smallest coordinate supported by PDB. The actual minimal supported coordinate is
/// -99.9999 nm but due to floating point shenanigans, we are slightly more restrictive to be safe.
pub(crate) const PDB_MIN_COORDINATE: f32 = -99.0;
/// Largest coordinate supported by PDB. The actual maximal supported coordinate is
/// 999.9999 nm but due to floating point shenanigans, we are slightly more restrictive to be safe.
pub(crate) const PDB_MAX_COORDINATE: f32 = 999.0;

/// Smallest coordinate supported by GRO. The actual minimal supported coordinate is
/// -999.999 nm but due to floating point shenanigans, we are slightly more restrictive to be safe.
pub(crate) const GRO_MIN_COORDINATE: f32 = -999.0;
/// Largest coordinate supported by GRO. The actual maximal supported coordinate is
/// 9999.999 nm but due to floating point shenanigans, we are slightly more restrictive to be safe.
pub(crate) const GRO_MAX_COORDINATE: f32 = 9999.0;

/******************************/
/*      GROUP/LABEL NAMES     */
/******************************/

/// Check whether the name for the group is a valid group name or a label name.
/// Characters '"&|!@()<>= are not allowed. Names containing whitespace only are also not allowed.
pub(crate) fn name_is_valid(string: &str) -> bool {
    if string.trim().is_empty() {
        return false;
    }

    let forbidden_chars = "'\"&|!@()<>=";

    for c in string.chars() {
        if forbidden_chars.contains(c) {
            return false;
        }
    }

    true
}

/******************************/
/*   CENTER OF MASS/GEOMETRY  */
/******************************/

/// Calculate contribution of an atom to the center of mass.
#[inline(always)]
pub(crate) fn center_atom_contribution(
    mut position: Vector3D,
    scaling: &Vector3D,
    simbox: &SimBox,
    mass: f32,
    sum_xi: &mut Vector3D,
    sum_zeta: &mut Vector3D,
) {
    // wrap position into the box
    position.wrap(simbox);

    let theta = Vector3D::new(
        position.x * scaling.x,
        position.y * scaling.y,
        position.z * scaling.z,
    );

    sum_xi.x += mass * theta.x.cos();
    sum_xi.y += mass * theta.y.cos();
    sum_xi.z += mass * theta.z.cos();

    sum_zeta.x += mass * theta.x.sin();
    sum_zeta.y += mass * theta.y.sin();
    sum_zeta.z += mass * theta.z.sin();
}

/// Convert coordinates from their representation on a circle to their representation on a line.
#[inline(always)]
pub(crate) fn from_circle_to_line(zeta: Vector3D, xi: Vector3D, scaling: &Vector3D) -> Vector3D {
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
/*      RMSD CALCULATION      */
/******************************/

/// Calculate the optimal rotation matrix, translation vector, and RMSD
/// to align two sets of points (P -> Q) using the Kabsch algorithm.
///
/// ## Parameters
/// - `p`: vector of points
/// - `q`: another vector of points
/// - `centroid_p`: center of geometry (or mass) of points `p`
/// - `centroid_q`: center of geometry (or mass) of points `q`
///
/// ## Returns
/// A tuple containing the optimal rotation matrix, the optimal translation vector, and the RMSD.
///
/// ## Panics
/// Panics if the number of points `p` does not match the number of points `q`.
pub(crate) fn kabsch_rmsd(
    p: &[Vector3D],
    q: &[Vector3D],
    centroid_p: &Vector3D,
    centroid_q: &Vector3D,
) -> (Matrix3<f32>, Vector3D, f32) {
    assert_eq!(
        p.len(),
        q.len(),
        "FATAL GROAN ERROR | auxiliary::kabsch_rmsd | Number of points `p` and `q` does not match."
    );

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
        .map(|(p_r, q_c)| (p_r.deref() - q_c.deref()).norm_squared())
        .sum::<f32>()
        / p.len() as f32)
        .sqrt();

    // return the rotation matrix, translation vector, and RMSD
    (r, centroid_q - centroid_p, rmsd)
}

#[cfg(test)]
mod tests_rmsd {
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

        let (rotation_matrix, translation_vector, rmsd) = kabsch_rmsd(
            &p,
            &q,
            &Vector3D::new(0.3333333, 0.3333333, 0.3333333),
            &Vector3D::new(0.3333333, 0.3333333, 0.3333333),
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
            Vector3D::new(0.66666666, 1.0, 0.0),
            Vector3D::new(-0.3333333, 0.0, 0.0),
            Vector3D::new(0.6666666, 0.0, 1.0),
        ];

        let (rotation_matrix, translation_vector, rmsd) = kabsch_rmsd(
            &p,
            &q,
            &Vector3D::new(0.3333333, 0.3333333, 0.3333333),
            &Vector3D::new(0.3333333, 0.3333333, 0.3333333),
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

        let (rotation_matrix, translation_vector, rmsd) = kabsch_rmsd(
            &p,
            &q,
            &Vector3D::new(0.3333333, 0.3333333, 0.3333333),
            &Vector3D::new(1.3333333, 1.3333333, 1.3333333),
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
            Vector3D::new(1.66666666, 2.0, 1.0),
            Vector3D::new(0.6666666, 1.0, 1.0),
            Vector3D::new(1.6666666, 1.0, 2.0),
        ];

        let (rotation_matrix, translation_vector, rmsd) = kabsch_rmsd(
            &p,
            &q,
            &Vector3D::new(0.3333333, 0.3333333, 0.3333333),
            &Vector3D::new(1.3333333, 1.3333333, 1.3333333),
        );

        let expected_rotation_matrix =
            Matrix3::from([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        assert!(rotation_matrix.relative_eq(&expected_rotation_matrix, 1e-6, 1e-6));
        let expected_translation = Vector3D::new(1.0, 1.0, 1.0);
        assert!(translation_vector.relative_eq(&expected_translation, 1e-6, 1e-6));
        assert_approx_eq!(f32, rmsd, 0.0, epsilon = 1e-6);
    }
}
