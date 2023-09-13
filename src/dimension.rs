// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of the Dimension enum.

use std::fmt;

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Dimension {
    None,
    X,
    Y,
    Z,
    XY,
    XZ,
    YZ,
    XYZ,
}

impl From<[bool; 3]> for Dimension {
    /// Convert 3-member boolean array to Dimension enum.
    fn from(arr: [bool; 3]) -> Self {
        match arr {
            [false, false, false] => Dimension::None,
            [true, false, false] => Dimension::X,
            [false, true, false] => Dimension::Y,
            [false, false, true] => Dimension::Z,
            [true, true, false] => Dimension::XY,
            [true, false, true] => Dimension::XZ,
            [false, true, true] => Dimension::YZ,
            [true, true, true] => Dimension::XYZ,
        }
    }
}

impl From<Dimension> for [bool; 3] {
    /// Convert Dimension to a 3-member boolean array.
    fn from(dim: Dimension) -> Self {
        match dim {
            Dimension::None => [false, false, false],
            Dimension::X => [true, false, false],
            Dimension::Y => [false, true, false],
            Dimension::Z => [false, false, true],
            Dimension::XY => [true, true, false],
            Dimension::XZ => [true, false, true],
            Dimension::YZ => [false, true, true],
            Dimension::XYZ => [true, true, true],
        }
    }
}

impl fmt::Display for Dimension {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Dimension::None => write!(f, "None"),
            Dimension::X => write!(f, "X"),
            Dimension::Y => write!(f, "Y"),
            Dimension::Z => write!(f, "Z"),
            Dimension::XY => write!(f, "XY"),
            Dimension::XZ => write!(f, "XZ"),
            Dimension::YZ => write!(f, "YZ"),
            Dimension::XYZ => write!(f, "XYZ"),
        }
    }
}

impl Dimension {
    /// Return `true` if Dimension contains x-dimension.
    pub fn is_x(self) -> bool {
        matches!(self, Dimension::X | Dimension::XY | Dimension::XZ | Dimension::XYZ)
    }

    /// Return `true` if Dimension contains y-dimension.
    pub fn is_y(self) -> bool {
        matches!(self, Dimension::Y | Dimension::XY | Dimension::YZ | Dimension::XYZ)
    }

    /// Return `true` if Dimension contains z-dimension.
    pub fn is_z(self) -> bool {
        matches!(self, Dimension::Z | Dimension::XZ | Dimension::YZ | Dimension::XYZ)
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests_dimension {
    use super::*;

    #[test]
    fn from_array_x() {
        let dim: Dimension = [true, false, false].into();
        assert_eq!(dim, Dimension::X);
        assert!(dim.is_x());
        assert!(!dim.is_y());
        assert!(!dim.is_z());
    }

    #[test]
    fn from_array_y() {
        let dim: Dimension = [false, true, false].into();
        assert_eq!(dim, Dimension::Y);
        assert!(!dim.is_x());
        assert!(dim.is_y());
        assert!(!dim.is_z());
    }

    #[test]
    fn from_array_z() {
        let dim: Dimension = [false, false, true].into();
        assert_eq!(dim, Dimension::Z);
        assert!(!dim.is_x());
        assert!(!dim.is_y());
        assert!(dim.is_z());
    }

    #[test]
    fn from_array_xy() {
        let dim: Dimension = [true, true, false].into();
        assert_eq!(dim, Dimension::XY);
        assert!(dim.is_x());
        assert!(dim.is_y());
        assert!(!dim.is_z());
    }

    #[test]
    fn from_array_xz() {
        let dim: Dimension = [true, false, true].into();
        assert_eq!(dim, Dimension::XZ);
        assert!(dim.is_x());
        assert!(!dim.is_y());
        assert!(dim.is_z());
    }

    #[test]
    fn from_array_yz() {
        let dim: Dimension = [false, true, true].into();
        assert_eq!(dim, Dimension::YZ);
        assert!(!dim.is_x());
        assert!(dim.is_y());
        assert!(dim.is_z());
    }

    #[test]
    fn from_array_xyz() {
        let dim: Dimension = [true, true, true].into();
        assert_eq!(dim, Dimension::XYZ);
        assert!(dim.is_x());
        assert!(dim.is_y());
        assert!(dim.is_z());
    }

    #[test]
    fn from_array_none() {
        let dim: Dimension = [false, false, false].into();
        assert_eq!(dim, Dimension::None);
        assert!(!dim.is_x());
        assert!(!dim.is_y());
        assert!(!dim.is_z());
    }

    #[test]
    fn to_array_x() {
        let dim: [bool; 3] = Dimension::X.into();
        assert_eq!(dim, [true, false, false]);
    }

    #[test]
    fn to_array_y() {
        let dim: [bool; 3] = Dimension::Y.into();
        assert_eq!(dim, [false, true, false]);
    }

    #[test]
    fn to_array_z() {
        let dim: [bool; 3] = Dimension::Z.into();
        assert_eq!(dim, [false, false, true]);
    }

    #[test]
    fn to_array_xy() {
        let dim: [bool; 3] = Dimension::XY.into();
        assert_eq!(dim, [true, true, false]);
    }

    #[test]
    fn to_array_xz() {
        let dim: [bool; 3] = Dimension::XZ.into();
        assert_eq!(dim, [true, false, true]);
    }

    #[test]
    fn to_array_yz() {
        let dim: [bool; 3] = Dimension::YZ.into();
        assert_eq!(dim, [false, true, true]);
    }

    #[test]
    fn to_array_xyz() {
        let dim: [bool; 3] = Dimension::XYZ.into();
        assert_eq!(dim, [true, true, true]);
    }

    #[test]
    fn to_array_none() {
        let dim: [bool; 3] = Dimension::None.into();
        assert_eq!(dim, [false, false, false]);
    }
}
