// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of methods for the calculation of RMSD.

use std::ops::Deref;

use crate::{
    errors::{GroupError, RMSDError},
    prelude::{SimBox, TrajRead, TrajReader},
    structures::{atom::Atom, vector3d::Vector3D},
    system::System,
};
use nalgebra::Matrix3;

impl System {
    /// Calculate the Root Mean Square Deviation (RMSD) between the specified group of atoms
    /// in the current system and the reference system.
    ///
    /// Uses the mass-weighted Kabsch algorithm.
    ///
    /// ## Returns
    /// - If successful, returns the minimum RMSD (in nanometers) between the atoms
    ///   of the specified group in `reference` and `self`.
    /// - Returns an `RMSDError` if the calculation fails.
    ///
    /// ## Example
    /// Calculating RMSD for every 10th frame of an XTC trajectory.
    /// Note that it is much faster to use [`RMSDIterator`] for this kind of calculation.
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
    /// - Atoms for which RMSD is calculated must be assigned masses in both systems.
    /// - The method performs a rigid-body alignment of the atoms in the specified group using
    ///   the Kabsch algorithm before calculating the RMSD.
    /// - Mass weighting **is** performed during the alignment. The reference structure must be assigned masses.
    /// - The RMSD is calculated in nanometers (nm).
    /// - Neither the current system (`self`) nor the reference system is modified by this method.
    /// - The group does not have to be centered or fully contained inside the simulation box.
    #[inline]
    pub fn calc_rmsd(&self, reference: &System, group: &str) -> Result<f32, RMSDError> {
        let (_, _, rmsd) = self.calc_rmsd_rot_trans(reference, group)?;
        Ok(rmsd)
    }

    /// Fits the selected `group` from the system to the `reference` structure.
    /// Also returns the Root Mean Square Deviation (RMSD) between the specified group of atoms
    /// in the current system and the reference system.
    ///
    /// Uses the mass-weighted Kabsch algorithm.
    ///
    /// ## Returns
    /// - If successful, returns the minimum RMSD (in nanometers) between the atoms
    ///   of the specified group in `reference` and `self`.
    /// - Returns an `RMSDError` if the calculation fails. In this case, `self` is not modified.
    ///
    /// ## Example
    /// Fitting a protein structure to a reference:
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// fn rmsd_fit() -> Result<System, Box<dyn std::error::Error + Send + Sync>> {
    ///     // load the reference structure, which contains only the atoms of the protein
    ///     // (in general, the reference may of course also contain other atoms)
    ///     let mut reference = System::from_file("protein.gro")?;
    ///     reference.group_create("Protein", "all")?;
    ///
    ///     // load the target structure, which may also contain other atoms besides the protein
    ///     let mut system = System::from_file("protein_in_membrane.gro")?;
    ///     system.group_create("Protein", "@protein")?;
    ///
    ///     // fit the system to the reference
    ///     let rmsd = system.calc_rmsd_and_fit(&reference, "Protein")?;
    ///
    ///     println!("{}", rmsd);
    ///
    ///     Ok(system)
    /// }
    /// ```
    /// ## Notes
    /// - This method requires both systems to have a valid orthogonal simulation box;
    ///   otherwise, an error will be returned.
    /// - Atoms for which RMSD is calculated must be assigned masses in both systems.
    /// - The method performs a rigid-body alignment of the atoms in the specified group using
    ///   the Kabsch algorithm before calculating the RMSD.
    /// - Mass weighting **is** performed during the alignment.
    /// - The RMSD is calculated in nanometers (nm).
    /// - The reference system is not modified by this function.
    /// - The group does not need to be centered or fully contained within the simulation box.
    #[inline]
    pub fn calc_rmsd_and_fit(&mut self, reference: &System, group: &str) -> Result<f32, RMSDError> {
        let (rotation, _, rmsd) = self.calc_rmsd_rot_trans(reference, group)?;

        let system_box = self.get_box_copy().unwrap();
        let box_center = self.get_box_center().unwrap();
        let group_com = self.group_get_com(group).unwrap();
        let ref_group_com = reference.group_get_com(group).unwrap();
        let inverse_box_center = Vector3D::new(-box_center.x, -box_center.y, -box_center.z);
        let shift_to_box_center = box_center - group_com;

        self.atoms_iter_mut().for_each(|atom| {
            atom.translate(&shift_to_box_center, &system_box).unwrap();

            atom.translate_nopbc(&inverse_box_center).unwrap();
            atom.rotate_nopbc(&rotation).unwrap();
            atom.translate_nopbc(&ref_group_com).unwrap();
        });

        Ok(rmsd)
    }

    /// Calculate RMSD, rotation matrix and translation vector.
    fn calc_rmsd_rot_trans(
        &self,
        reference: &System,
        group: &str,
    ) -> Result<(Matrix3<f32>, Vector3D, f32), RMSDError> {
        let (reference_coordinates, reference_box_center) =
            extract_data_from_system(reference, group)?;
        let (target_coordinates, target_box_center) = extract_data_from_system(self, group)?;

        // check that the group has a consistent number of atoms
        number_of_positions_consistent(&reference_coordinates, &target_coordinates, group)?;

        // extract masses from reference
        let masses = extract_masses(reference.group_iter(group));
        let sum_masses = masses.iter().sum::<f32>();

        // calculate RMSD
        Ok(kabsch_rmsd(
            &reference_coordinates,
            &target_coordinates,
            &masses,
            &reference_box_center,
            &target_box_center,
            sum_masses,
        ))
    }
}

/// Structure for efficient calculation of RMSD from trajectories.
/// `RMSDIterator` is a trajectory converter that returns both the
/// current state of the system and the current RMSD value.
pub struct RMSDIterator<'a, R>
where
    R: TrajRead<'a>,
{
    // precomputed, extracted and shifted reference coordinates
    reference_coordinates: Vec<Vector3D>,
    // center of the reference simulation box
    reference_center: Vector3D,
    // extracted masses of the relevant atoms
    masses: Vec<f32>,
    // sum of the extracted masses
    sum_masses: f32,
    // name of the group to work with
    group: String,
    // trajectory reader to work with
    reader: TrajReader<'a, R>,
}

impl<'a, R> TrajReader<'a, R>
where
    R: TrajRead<'a>,
{
    pub fn calc_rmsd(
        self,
        reference: &System,
        group: &str,
    ) -> Result<RMSDIterator<'a, R>, RMSDError> {
        let (coordinates, center) = extract_data_from_system(reference, group)?;
        let masses = extract_masses(reference.group_iter(group));
        let sum_masses = masses.iter().sum::<f32>();

        Ok(RMSDIterator {
            reference_coordinates: coordinates,
            reference_center: center,
            masses,
            sum_masses,
            group: group.to_owned(),
            reader: self,
        })
    }
}

impl<'a, R: TrajRead<'a>> Iterator for RMSDIterator<'a, R> {
    type Item = Result<(&'a mut System, f32), RMSDError>;

    fn next(&mut self) -> Option<Self::Item> {
        let frame = match self.reader.next()? {
            Ok(x) => x,
            Err(e) => return Some(Err(RMSDError::TrajectoryError(e))),
        };

        let (target_coordinates, target_center) = match extract_data_from_system(frame, &self.group)
        {
            Ok(x) => x,
            Err(e) => return Some(Err(e)),
        };

        match number_of_positions_consistent(
            &self.reference_coordinates,
            &target_coordinates,
            &self.group,
        ) {
            Ok(x) => x,
            Err(e) => return Some(Err(e)),
        };

        let (_, _, rmsd) = kabsch_rmsd(
            &self.reference_coordinates,
            &target_coordinates,
            &self.masses,
            &self.reference_center,
            &target_center,
            self.sum_masses,
        );

        Some(Ok((frame, rmsd)))
    }
}

/// Auxiliary function for checking consistency in the number of positions.
#[inline(always)]
fn number_of_positions_consistent(
    reference: &Vec<Vector3D>,
    target: &Vec<Vector3D>,
    group: &str,
) -> Result<(), RMSDError> {
    let n_atoms_ref = reference.len();
    let n_atoms_target = target.len();

    if n_atoms_ref != n_atoms_target {
        Err(RMSDError::InconsistentGroup(
            group.to_owned(),
            n_atoms_ref,
            n_atoms_target,
        ))
    } else {
        Ok(())
    }
}

/// Auxiliary function for extracting data required for the RMSD calculation from system.
fn extract_data_from_system(
    system: &System,
    group: &str,
) -> Result<(Vec<Vector3D>, Vector3D), RMSDError> {
    // check that the simulation box is defined and orthogonal & extract box center
    let box_center = system.get_box_center().map_err(RMSDError::InvalidSimBox)?;

    // get the center of mass of the group
    let group_center = get_com(system, group)?;

    // extract the coordinates of the atoms of the group and shift them
    let shift = Vector3D(box_center.deref() - group_center.deref());
    let coordinates = shift_and_wrap_coordinates(
        extract_coordinates(system.group_iter(group)),
        &shift,
        system.get_box().expect(
            "FATAL GROAN ERROR | rmsd::extract_data_from_system | Simulation box should exist.",
        ),
    );

    Ok((coordinates, box_center))
}

/// Auxiliary function for extracting masses of relevant atoms from system.
#[inline(always)]
fn extract_masses<'a, I>(iter: Result<I, GroupError>) -> Vec<f32>
where
    I: Iterator<Item = &'a Atom>,
{
    iter.expect("FATAL GROAN ERROR | rmsd::extract_masses | Group should exist.")
        .map(|atom| {
            atom.get_mass()
                .expect("FATAL GROAN ERROR | System::calc_rmsd_rot_trans | Mass should be defined.")
        })
        .collect::<Vec<f32>>()
}

/// Auxiliary function for extracting the coordinates of relevant atoms from system.
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
        Err(GroupError::NotFound(x)) => Err(RMSDError::NonexistentGroup(x)),
        Err(GroupError::InvalidPosition(x)) => Err(RMSDError::InvalidPosition(x)),
        Err(GroupError::InvalidMass(x)) => Err(RMSDError::InvalidMass(x)),
        Err(GroupError::EmptyGroup(x)) => Err(RMSDError::EmptyGroup(x)),
        Err(_) => panic!("FATAL GROAN ERROR | rmsd::get_center | Unexpected error type returned from System::group_get_com."),
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
    use std::fs::File;

    use crate::errors::{MassError, PositionError, SimBoxError};

    use super::*;
    use float_cmp::{approx_eq, assert_approx_eq};
    use tempfile::NamedTempFile;

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
            0.23680364, 0.26356384, 0.26030704, 0.2139618, 0.2221243, 0.19429797, 0.26472777,
            0.27031714, 0.26426855, 0.23497728, 0.2426188,
        ];

        assert_eq!(rmsd.len(), expected.len());
        for (x, y) in rmsd.into_iter().zip(expected.into_iter()) {
            assert_approx_eq!(f32, x, y);
        }
    }

    #[test]
    fn test_rmsd_iterator() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let reference = system.clone();

        let rmsd = system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .map(|frame| frame.unwrap().calc_rmsd(&reference, "Protein").unwrap())
            .collect::<Vec<f32>>();

        let rmsd2 = system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .calc_rmsd(&reference, "Protein")
            .unwrap()
            .map(|result| result.unwrap().1)
            .collect::<Vec<f32>>();

        assert_eq!(rmsd.len(), rmsd2.len());
        for (x, y) in rmsd.into_iter().zip(rmsd2.into_iter()) {
            assert_approx_eq!(f32, x, y);
        }
    }

    #[test]
    fn test_calc_rmsd_partial() {
        // the reference system does not have to contain all the atoms of the target system
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let mut reference = System::new(
            "Reference system",
            system.group_extract("Protein").unwrap(),
            system.get_box_copy(),
        );
        reference.group_create("Protein", "all").unwrap();

        assert_approx_eq!(
            f32,
            system.calc_rmsd(&reference, "Protein").unwrap(),
            0.0,
            epsilon = 1e-4
        );
    }

    #[test]
    fn test_calc_rmsd_broken_at_pbc() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let mut reference = system.clone();
        // break peptide at PBC
        reference
            .atoms_translate(&Vector3D::new(3.2, -2.1, -4.6))
            .unwrap();

        assert_approx_eq!(
            f32,
            system.calc_rmsd(&reference, "Protein").unwrap(),
            0.0,
            epsilon = 1e-4
        );
        assert_approx_eq!(
            f32,
            reference.calc_rmsd(&system, "Protein").unwrap(),
            0.0,
            epsilon = 1e-4
        );
    }

    fn compare_positions_with_wrapping(pos1: &Vector3D, pos2: &Vector3D, simbox: &SimBox) {
        if !approx_eq!(f32, pos1.x, pos2.x, epsilon = 1e-3)
            && !approx_eq!(f32, pos1.x + simbox.x, pos2.x, epsilon = 1e-3)
            && !approx_eq!(f32, pos1.x - simbox.x, pos2.x, epsilon = 1e-3)
        {
            panic!(
                "Assertion failed: {} != {} even after wrapping",
                pos1.x, pos2.x
            );
        }

        if !approx_eq!(f32, pos1.y, pos2.y, epsilon = 1e-3)
            && !approx_eq!(f32, pos1.y + simbox.y, pos2.y, epsilon = 1e-3)
            && !approx_eq!(f32, pos1.y - simbox.y, pos2.y, epsilon = 1e-3)
        {
            panic!(
                "Assertion failed: {} != {} even after wrapping",
                pos1.y, pos2.y
            );
        }

        if !approx_eq!(f32, pos1.z, pos2.z, epsilon = 1e-3)
            && !approx_eq!(f32, pos1.z + simbox.z, pos2.z, epsilon = 1e-3)
            && !approx_eq!(f32, pos1.z - simbox.z, pos2.z, epsilon = 1e-3)
        {
            panic!(
                "Assertion failed: {} != {} even after wrapping",
                pos1.z, pos2.z
            );
        }
    }

    #[test]
    fn test_rmsd_fit_same_structure() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let reference = system.clone();

        let rmsd = system.calc_rmsd_and_fit(&reference, "Protein").unwrap();
        assert_approx_eq!(f32, rmsd, 0.0, epsilon = 1e-4);

        let simbox = system.get_box().unwrap();

        for (pos1, pos2) in reference
            .atoms_iter()
            .map(|atom| atom.get_position().unwrap())
            .zip(system.atoms_iter().map(|atom| atom.get_position().unwrap()))
        {
            compare_positions_with_wrapping(pos1, pos2, simbox);
        }
    }

    #[test]
    fn test_rmsd_fit_shifted_and_rotated_copy() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let reference = system.clone();

        system
            .atoms_translate(&Vector3D::new(-1.1, 3.4, 2.7))
            .unwrap();
        let rotation_matrix = Matrix3::new(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        system
            .atoms_iter_mut()
            .for_each(|atom| atom.rotate_nopbc(&rotation_matrix).unwrap());

        let rmsd = system.calc_rmsd_and_fit(&reference, "Protein").unwrap();
        assert_approx_eq!(f32, rmsd, 0.0, epsilon = 1e-4);

        let simbox = system.get_box().unwrap();

        for (pos1, pos2) in reference
            .atoms_iter()
            .map(|atom| atom.get_position().unwrap())
            .zip(system.atoms_iter().map(|atom| atom.get_position().unwrap()))
        {
            compare_positions_with_wrapping(pos1, pos2, simbox);
        }
    }

    #[test]
    fn test_rmsd_fit_trajectory() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let mut reference = system.clone();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();
        system.xtc_writer_init(path_to_output).unwrap();

        // should work even if we remove position of some atom that is not in the group
        reference.get_atom_mut(176).unwrap().reset_position();

        let rmsd = system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .map(|frame| {
                let frame = frame.unwrap();
                let rmsd = frame.calc_rmsd_and_fit(&reference, "Protein").unwrap();
                frame.traj_write_frame().unwrap();

                rmsd
            })
            .collect::<Vec<f32>>();

        // close the output xtc file
        system.traj_close();

        let expected = [
            0.23680364, 0.26356384, 0.26030704, 0.2139618, 0.2221243, 0.19429797, 0.26472777,
            0.27031714, 0.26426855, 0.23497728, 0.2426188,
        ];

        assert_eq!(rmsd.len(), expected.len());
        for (x, y) in rmsd.into_iter().zip(expected.into_iter()) {
            assert_approx_eq!(f32, x, y);
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory_fit.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn test_rmsd_fit_trajectory_partial() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let mut reference = System::new(
            "Reference system",
            system.group_extract("Protein").unwrap(),
            system.get_box_copy(),
        );
        reference.group_create("Protein", "@protein").unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();
        system.xtc_writer_init(path_to_output).unwrap();

        let rmsd = system
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .map(|frame| {
                let frame = frame.unwrap();
                let rmsd = frame.calc_rmsd_and_fit(&reference, "Protein").unwrap();
                frame.traj_write_frame().unwrap();

                rmsd
            })
            .collect::<Vec<f32>>();

        // close the output xtc file
        system.traj_close();

        let expected = [
            0.23680364, 0.26356384, 0.26030704, 0.2139618, 0.2221243, 0.19429797, 0.26472777,
            0.27031714, 0.26426855, 0.23497728, 0.2426188,
        ];

        assert_eq!(rmsd.len(), expected.len());
        for (x, y) in rmsd.into_iter().zip(expected.into_iter()) {
            assert_approx_eq!(f32, x, y);
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory_fit.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn test_rmsd_fit_trajectory_broken_at_pbc() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let mut reference = system.clone();
        // break peptide at PBC
        reference
            .atoms_translate(&Vector3D::new(3.2, -2.1, -4.6))
            .unwrap();

        let xtc_output = NamedTempFile::new().unwrap();
        let path_to_output = xtc_output.path();
        system.xtc_writer_init(path_to_output).unwrap();

        for frame in system.xtc_iter("test_files/short_trajectory.xtc").unwrap() {
            let frame = frame.unwrap();
            frame.calc_rmsd_and_fit(&reference, "Protein").unwrap();
            frame.traj_write_frame().unwrap();
        }

        // close the output xtc file
        system.traj_close();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/short_trajectory_broken_fit.xtc").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn test_calc_rmsd_fail_group_does_not_exist_in_reference() {
        let mut system = System::from_file("test_files/example.tpr").unwrap();
        let reference = system.clone();

        system.group_create("Protein", "@protein").unwrap();

        match system.calc_rmsd_rot_trans(&reference, "Protein") {
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

        match system.calc_rmsd_rot_trans(&reference, "Protein") {
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

        match system.calc_rmsd_rot_trans(&reference, "Protein") {
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

        match system.calc_rmsd_rot_trans(&reference, "Protein") {
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

        match system.calc_rmsd_rot_trans(&reference, "Protein") {
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

        match system.calc_rmsd_rot_trans(&reference, "Protein") {
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

        match system.calc_rmsd_rot_trans(&reference, "Protein") {
            Ok(_) => panic!("Function should have failed."),
            Err(RMSDError::InvalidMass(MassError::NoMass(_))) => (),
            Err(_) => panic!("Function failed but incorrect error type returned."),
        }
    }
}
