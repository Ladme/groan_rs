// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of the SimBox structure and its methods.

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
    fn from(arr: [f32; 9]) -> Self {
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
    /// Check that the simulation box is orthogonal.
    /// ## Returns
    /// `true` if only the first three dimensions of the `SimBox` are non-zero. Otherwise, returns `false`.
    pub fn is_orthogonal(&self) -> bool {
        self.v1y == 0.0
            && self.v1z == 0.0
            && self.v2x == 0.0
            && self.v2z == 0.0
            && self.v3x == 0.0
            && self.v3y == 0.0
    }

    /// Check that the box is a valid simulation box, i.e. at least the first 3 dimensions of the box are positive.
    /// ## Returns
    /// `true` if the box is a valid simulation box. Otherwise, returns `false`.
    pub fn is_valid(&self) -> bool {
        self.x > 0.0 && self.y > 0.0 && self.z > 0.0
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
