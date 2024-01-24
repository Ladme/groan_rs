// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of the SimBox structure and its methods.

use crate::{errors::SimBoxError, structures::vector3d::Vector3D};
use std::ops::Deref;

/// Structure defining simulation box shape and dimensions.
#[derive(Debug, Clone)]
pub struct SimBox {
    /// You can also use `.x` to reach this value.
    pub v1x: f32,
    /// You can also use `.y` to reach this value.
    pub v2y: f32,
    /// You can also use `.z` to reach this value.
    pub v3z: f32,
    pub v1y: f32,
    pub v1z: f32,
    pub v2x: f32,
    pub v2z: f32,
    pub v3x: f32,
    pub v3y: f32,
}

impl From<[f32; 9]> for SimBox {
    /// Convert 9-member array to SimBox structure.
    /// The order of the members of the array should be the same as in a gro file.
    ///
    /// ## Panics
    /// Panics if the `SimBox` is not a simulation box supported by Gromacs,
    /// i.e. if `v1y`, `v1z`, and `v2z` are not zero.
    fn from(arr: [f32; 9]) -> Self {
        if arr[3] != 0.0 || arr[4] != 0.0 || arr[6] != 0.0 {
            panic!("FATAL GROAN ERROR | SimBox::from | Unsupported Gromacs simulation box.");
        }

        SimBox {
            v1x: arr[0],
            v2y: arr[1],
            v3z: arr[2],
            v1y: arr[3],
            v1z: arr[4],
            v2x: arr[5],
            v2z: arr[6],
            v3x: arr[7],
            v3y: arr[8],
        }
    }
}

impl From<[f32; 3]> for SimBox {
    /// Convert 3-member array to SimBox structure. Last 6 values of SimBox are set to 0.
    /// The order of the members of the array should be the same as in a gro file.
    fn from(arr: [f32; 3]) -> Self {
        SimBox {
            v1x: arr[0],
            v2y: arr[1],
            v3z: arr[2],
            v1y: 0.0f32,
            v1z: 0.0f32,
            v2x: 0.0f32,
            v2z: 0.0f32,
            v3x: 0.0f32,
            v3y: 0.0f32,
        }
    }
}

impl SimBox {
    /// Create new simulation box from lengths and angles (in degrees).
    ///
    /// ## Example
    /// ```
    /// # use groan_rs::prelude::*;
    /// # use float_cmp::assert_approx_eq;
    /// #
    /// let simbox = SimBox::from_lengths_angles([5.0, 4.0, 3.0].into(), [80.0, 70.0, 120.0].into());
    ///
    /// assert_approx_eq!(f32, simbox.v1x,  5.000000, epsilon = 0.0001);
    /// assert_approx_eq!(f32, simbox.v2y,  3.464102, epsilon = 0.0001);
    /// assert_approx_eq!(f32, simbox.v3z,  2.553768, epsilon = 0.0001);
    /// assert_approx_eq!(f32, simbox.v1y,  0.000000, epsilon = 0.0001);
    /// assert_approx_eq!(f32, simbox.v1z,  0.000000, epsilon = 0.0001);
    /// assert_approx_eq!(f32, simbox.v2x, -2.000000, epsilon = 0.0001);
    /// assert_approx_eq!(f32, simbox.v2z,  0.000000, epsilon = 0.0001);
    /// assert_approx_eq!(f32, simbox.v3x,  1.026060, epsilon = 0.0001);
    /// assert_approx_eq!(f32, simbox.v3y,  1.193930, epsilon = 0.0001);
    /// ```
    ///
    /// ## Notes
    /// - Adapted from Tsjerk Wassenaar's `triclinic` function:
    /// <https://www.mail-archive.com/gmx-users@gromacs.org/msg28032.html>
    pub fn from_lengths_angles(lengths: Vector3D, angles: Vector3D) -> Self {
        let mut simbox = SimBox {
            v1x: lengths.x,
            ..Default::default()
        };

        if angles.x == 90.0 && angles.y == 90.0 && angles.z == 90.0 {
            simbox.v2y = lengths.y;
            simbox.v3z = lengths.z;
        } else {
            // convert to radians
            let alpha = angles.x * std::f32::consts::PI / 180.0;
            let beta = angles.y * std::f32::consts::PI / 180.0;
            let gamma = angles.z * std::f32::consts::PI / 180.0;

            simbox.v2x = lengths.y * gamma.cos();
            simbox.v2y = lengths.y * gamma.sin();

            simbox.v3x = lengths.z * beta.cos();
            simbox.v3y = lengths.z * (alpha.cos() - beta.cos() * gamma.cos()) / gamma.sin();
            simbox.v3z =
                (lengths.z * lengths.z - simbox.v3x * simbox.v3x - simbox.v3y * simbox.v3y).sqrt();

            // v1y, v1z, v2z are guaranteed to be zero as required by Gromacs
        }

        simbox
    }

    /// Calculate box lengths and angles from simulation box.
    ///
    /// ## Returns
    /// (`lengths`, `angles`) of the simulation box.
    ///
    /// ## Example
    /// ```
    /// # use groan_rs::prelude::*;
    /// # use float_cmp::assert_approx_eq;
    /// #
    /// let simbox = SimBox::from([
    ///     5.0, 3.464102, 2.553768,
    ///     0.0, 0.0, -2.000000,
    ///     0.000000, 1.026060, 1.193930]);
    ///
    /// let (lengths, angles) = simbox.to_lengths_angles();
    ///
    /// assert_approx_eq!(f32, angles.x,  80.0, epsilon = 0.0001);
    /// assert_approx_eq!(f32, angles.y,  70.0, epsilon = 0.0001);
    /// assert_approx_eq!(f32, angles.z, 120.0, epsilon = 0.0001);
    ///
    /// assert_approx_eq!(f32, lengths.x, 5.0, epsilon = 0.0001);
    /// assert_approx_eq!(f32, lengths.y, 4.0, epsilon = 0.0001);
    /// assert_approx_eq!(f32, lengths.z, 3.0, epsilon = 0.0001);
    /// ```
    pub fn to_lengths_angles(&self) -> (Vector3D, Vector3D) {
        if self.is_orthogonal() {
            // orthogonal simulation box
            let lengths = Vector3D::new(
                self.v1x,
                self.v2y,
                self.v3z,
            );

            let angles = Vector3D::new(
                90.0,
                90.0,
                90.0,
            );

            (lengths, angles)
        } else {
            // triclinic simulation box
            let gamma = self.v2y.atan2(self.v2x);

            let lengths = Vector3D::new(
                self.v1x,
                (self.v2x * self.v2x + self.v2y * self.v2y).sqrt(),
                (self.v3x * self.v3x + self.v3y * self.v3y + self.v3z * self.v3z).sqrt(),
            );

            let beta = (self.v3x / lengths.z).acos();
            let alpha = ((self.v3y * gamma.sin()) / lengths.z + beta.cos() * gamma.cos()).acos();

            let angles = Vector3D::new(
                180.0 * alpha / std::f32::consts::PI,
                180.0 * beta / std::f32::consts::PI,
                180.0 * gamma / std::f32::consts::PI,
            );

            (lengths, angles)
        }
    }

    /// Check that the simulation box is orthogonal.
    ///
    /// ## Returns
    /// `true` if the simulation box is orthogonal. Otherwise, returns `false`.
    pub fn is_orthogonal(&self) -> bool {
        // we do not need to check v1y, v1z, and v2z as these must be zero
        self.v2x == 0.0 && self.v3x == 0.0 && self.v3y == 0.0
    }

    /// Check whether all dimensions of the simulation box are zero.
    pub fn is_zero(&self) -> bool {
        // we do not need to check v1y, v1z, and v2z as these must be zero
        self.x == 0.0 && self.y == 0.0 && self.z == 0.0 && self.is_orthogonal()
    }
}

impl Default for SimBox {
    /// Create a zero simulation box.
    fn default() -> Self {
        SimBox {
            v1x: 0.0f32,
            v2y: 0.0f32,
            v3z: 0.0f32,
            v1y: 0.0f32,
            v1z: 0.0f32,
            v2x: 0.0f32,
            v2z: 0.0f32,
            v3x: 0.0f32,
            v3y: 0.0f32,
        }
    }
}

/// Allows using .x, .y, and .z to reach the v1x, v2y, and v3z members of SimBox.
pub struct SimBoxDimensions {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Deref for SimBox {
    type Target = SimBoxDimensions;

    fn deref(&self) -> &Self::Target {
        unsafe { &*(self as *const SimBox as *const SimBoxDimensions) }
    }
}

/// Checks whether the simulation box exists and whether it is orthogonal.
pub(crate) fn simbox_check(simbox: Option<&SimBox>) -> Result<&SimBox, SimBoxError> {
    match simbox {
        Some(x) if x.is_orthogonal() => Ok(x),
        Some(_) => Err(SimBoxError::NotOrthogonal),
        None => Err(SimBoxError::DoesNotExist),
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
    fn from_lengths_angles_1() {
        let simbox =
            SimBox::from_lengths_angles([5.297, 4.863, 2.976].into(), [90.0, 90.0, 90.0].into());

        assert_approx_eq!(f32, simbox.v1x, 5.297, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2y, 4.863, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3z, 2.976, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v1y, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v1z, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2x, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2z, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3x, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3y, 0.0, epsilon = 0.00001);

        // convert back to lengths and angles
        let (lengths, angles) = simbox.to_lengths_angles();
        assert_approx_eq!(f32, lengths.x, 5.297, epsilon = 0.0001);
        assert_approx_eq!(f32, lengths.y, 4.863, epsilon = 0.0001);
        assert_approx_eq!(f32, lengths.z, 2.976, epsilon = 0.0001);

        assert_approx_eq!(f32, angles.x, 90.0, epsilon = 0.0001);
        assert_approx_eq!(f32, angles.y, 90.0, epsilon = 0.0001);
        assert_approx_eq!(f32, angles.z, 90.0, epsilon = 0.0001);
    }

    #[test]
    fn from_lengths_angles_2() {
        let simbox =
            SimBox::from_lengths_angles([5.297, 4.863, 2.976].into(), [120.0, 70.0, 80.0].into());

        assert_approx_eq!(f32, simbox.v1x, 5.297, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2y, 4.78912, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3z, 2.2277796, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v1y, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v1z, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2x, 0.8444507, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2z, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3x, 1.0178516, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3y, -1.6904297, epsilon = 0.00001);

        // convert back to lengths and angles
        let (lengths, angles) = simbox.to_lengths_angles();
        assert_approx_eq!(f32, lengths.x, 5.297, epsilon = 0.0001);
        assert_approx_eq!(f32, lengths.y, 4.863, epsilon = 0.0001);
        assert_approx_eq!(f32, lengths.z, 2.976, epsilon = 0.0001);

        assert_approx_eq!(f32, angles.x, 120.0, epsilon = 0.0001);
        assert_approx_eq!(f32, angles.y, 70.0, epsilon = 0.0001);
        assert_approx_eq!(f32, angles.z, 80.0, epsilon = 0.0001);
    }

    #[test]
    fn from_lengths_angles_3() {
        let simbox = SimBox::from_lengths_angles(
            [6.26832, 6.26832, 6.26832].into(),
            [60.0, 60.0, 90.0].into(),
        );

        assert_approx_eq!(f32, simbox.v1x, 6.26832, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2y, 6.26832, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3z, 4.43237, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v1y, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v1z, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2x, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2z, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3x, 3.13416, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3y, 3.13416, epsilon = 0.00001);

        // convert back to lengths and angles
        let (lengths, angles) = simbox.to_lengths_angles();
        assert_approx_eq!(f32, lengths.x, 6.26832, epsilon = 0.0001);
        assert_approx_eq!(f32, lengths.y, 6.26832, epsilon = 0.0001);
        assert_approx_eq!(f32, lengths.z, 6.26832, epsilon = 0.0001);

        assert_approx_eq!(f32, angles.x, 60.0, epsilon = 0.0001);
        assert_approx_eq!(f32, angles.y, 60.0, epsilon = 0.0001);
        assert_approx_eq!(f32, angles.z, 90.0, epsilon = 0.0001);
    }

    #[test]
    fn from_lengths_angles_4() {
        let simbox = SimBox::from_lengths_angles(
            [6.26832, 6.26832, 6.26832].into(),
            [70.53, 109.47, 70.53].into(),
        );

        assert_approx_eq!(f32, simbox.v1x, 6.26832, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2y, 5.90987, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3z, 5.11825, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v1y, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v1z, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2x, 2.08931, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2z, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3x, -2.08931, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3y, 2.95467, epsilon = 0.00001);
    }

    #[test]
    fn from_array9() {
        let arr = [
            6.26832, 5.90987, 5.11825, 0.0, 0.0, 2.08931, 0.0, -2.08931, 2.95467,
        ];

        let simbox = SimBox::from(arr);

        assert_approx_eq!(f32, simbox.v1x, 6.26832, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2y, 5.90987, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3z, 5.11825, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v1y, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v1z, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2x, 2.08931, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v2z, 0.0, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3x, -2.08931, epsilon = 0.00001);
        assert_approx_eq!(f32, simbox.v3y, 2.95467, epsilon = 0.00001);
    }

    #[test]
    #[should_panic(
        expected = "FATAL GROAN ERROR | SimBox::from | Unsupported Gromacs simulation box."
    )]
    fn from_array9_panics_1() {
        let arr = [
            6.26832, 5.90987, 5.11825, 1.3, 0.0, 2.08931, 0.0, -2.08931, 2.95467,
        ];
        let _simbox = SimBox::from(arr);
    }

    #[test]
    #[should_panic(
        expected = "FATAL GROAN ERROR | SimBox::from | Unsupported Gromacs simulation box."
    )]
    fn from_array9_panics_2() {
        let arr = [
            6.26832, 5.90987, 5.11825, 0.0, 2.5634, 2.08931, 0.0, -2.08931, 2.95467,
        ];
        let _simbox = SimBox::from(arr);
    }

    #[test]
    #[should_panic(
        expected = "FATAL GROAN ERROR | SimBox::from | Unsupported Gromacs simulation box."
    )]
    fn from_array9_panics_3() {
        let arr = [
            6.26832, 5.90987, 5.11825, 0.0, 0.0, 2.08931, -1.322, -2.08931, 2.95467,
        ];
        let _simbox = SimBox::from(arr);
    }

    #[test]
    #[should_panic(
        expected = "FATAL GROAN ERROR | SimBox::from | Unsupported Gromacs simulation box."
    )]
    fn from_array9_panics_4() {
        let arr = [
            6.26832, 5.90987, 5.11825, 1.3, 2.5634, 2.08931, -1.322, -2.08931, 2.95467,
        ];
        let _simbox = SimBox::from(arr);
    }
}
