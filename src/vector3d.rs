// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of methods for three-dimensional vector.

use crate::dimension::Dimension;
use crate::simbox::SimBox;

/// Function replicating the behavior of Python '%'.
fn floor_mod(x: f32, y: f32) -> f32 {
    (x % y + y) % y
}

/// Reciprocal of square root of 2, i.e. 1/sqrt(2).
const REC_SQRT2: f32 = 0.7071068f32;
/// Reciprocal of square root of 3, i.e. 1/sqrt(3).
const REC_SQRT3: f32 = 0.5773503f32;

/// Describes length and orientation of a vector in space or a position of a point in space.
#[derive(Debug, Clone, PartialEq)]
pub struct Vector3D {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl From<[f32; 3]> for Vector3D {
    fn from(arr: [f32; 3]) -> Self {
        Vector3D {
            x: arr[0],
            y: arr[1],
            z: arr[2],
        }
    }
}

impl From<Dimension> for Vector3D {
    /// Get unit vector oriented along the specified dimension(s).
    /// If `dim` is `Dimension::None`, returns a null vector.
    fn from(dim: Dimension) -> Self {
        match dim {
            Dimension::None => [0.0, 0.0, 0.0].into(),
            Dimension::X => [1.0, 0.0, 0.0].into(),
            Dimension::Y => [0.0, 1.0, 0.0].into(),
            Dimension::Z => [0.0, 0.0, 1.0].into(),
            Dimension::XY => [REC_SQRT2, REC_SQRT2, 0.0].into(),
            Dimension::XZ => [REC_SQRT2, 0.0, REC_SQRT2].into(),
            Dimension::YZ => [0.0, REC_SQRT2, REC_SQRT2].into(),
            Dimension::XYZ => [REC_SQRT3, REC_SQRT3, REC_SQRT3].into(),
        }
    }
}

impl Vector3D {
    /// Calculate length of the vector.
    ///
    /// ## Example
    /// ```
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let vector = Vector3D::from([1.0, 2.0, 3.0]);
    /// assert_approx_eq!(f32, vector.len(), 3.741657);
    /// ```
    pub fn len(&self) -> f32 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Convert vector to unit vector.
    ///
    /// ## Notes
    /// - Returns itself when applied to a null vector.
    ///
    /// ## Example
    /// ```
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let vector = Vector3D::from([1.0, 2.0, 3.0]).to_unit();
    ///
    /// assert_approx_eq!(f32, vector.x, 0.2672612);
    /// assert_approx_eq!(f32, vector.y, 0.5345225);
    /// assert_approx_eq!(f32, vector.z, 0.8017837);
    /// assert_approx_eq!(f32, vector.len(), 1.0);
    /// ```
    pub fn to_unit(mut self) -> Vector3D {
        let length = self.len();

        if length > 0.0 {
            self.x /= length;
            self.y /= length;
            self.z /= length;
        }

        self
    }

    /// Invert the vector. (Reverse direction of the vector.)
    ///
    /// ## Example
    /// ```
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let vector = Vector3D::from([1.0, -2.0, 3.0]).invert();
    ///
    /// assert_approx_eq!(f32, vector.x, -1.0);
    /// assert_approx_eq!(f32, vector.y,  2.0);
    /// assert_approx_eq!(f32, vector.z, -3.0);
    /// ```
    pub fn invert(mut self) -> Vector3D {
        self.x *= -1.0;
        self.y *= -1.0;
        self.z *= -1.0;

        self
    }

    /// Calculate the dot product of two vectors.
    ///
    /// ## Example
    /// ```
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let vector1 = Vector3D::from([4.0, 2.0, -1.0]);
    /// let vector2 = Vector3D::from([1.0, -3.0, 2.0]);
    ///
    /// assert_approx_eq!(f32, vector1.dot(&vector2), -4.0);
    /// ```
    pub fn dot(&self, vector: &Vector3D) -> f32 {
        self.x * vector.x + self.y * vector.y + self.z * vector.z
    }

    /// Calculate the cross product of two vectors.
    ///
    /// ## Example
    /// ```
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let vector1 = Vector3D::from([4.0, 2.0, -1.0]);
    /// let vector2 = Vector3D::from([1.0, -3.0, 2.0]);
    ///
    /// let cross = vector1.cross(&vector2);
    ///
    /// assert_approx_eq!(f32, cross.x, 1.0);
    /// assert_approx_eq!(f32, cross.y, -9.0);
    /// assert_approx_eq!(f32, cross.z, -14.0);
    /// ```
    pub fn cross(&self, vector: &Vector3D) -> Vector3D {
        Vector3D::from([
            self.y * vector.z - self.z * vector.y,
            self.z * vector.x - self.x * vector.z,
            self.x * vector.y - self.y * vector.x,
        ])
    }

    /// Calculate the angle between two vectors. Returns angle in radians.
    ///
    /// ## Example
    /// ```
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let vector1 = Vector3D::from([1.0, 2.0, 3.0]);
    /// let vector2 = Vector3D::from([3.0, 2.0, 1.0]);
    ///
    /// let angle = vector1.angle(&vector2);
    ///
    /// assert_approx_eq!(f32, angle, 0.77519345);
    /// ```
    ///
    /// ## Notes
    /// - Always returns a value between 0 and π.
    pub fn angle(&self, vector: &Vector3D) -> f32 {
        (self.dot(vector) / (self.len() * vector.len())).acos()
    }

    /// Shift vector in the direction of `orientation` by target distance.
    /// `orientation` does not have to be a unit vector.
    ///
    /// ## Warning
    /// Does not take periodic boundary conditions into consideration.
    /// If you want the vector to fit into simulation box, apply
    /// `Vector3D::wrap` after shifting.
    ///
    /// ## Example
    /// Shift a given point in space by 2 nm (in total) along x and y dimensions.
    /// ```
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let mut vector = Vector3D::from([1.0, 2.0, 3.0]);
    /// let orientation = Dimension::XY.into();
    ///
    /// // shifting along each of the x and y dimensions by √2 nm
    /// vector.shift(orientation, 2.0);
    ///
    /// assert_approx_eq!(f32, vector.x, 2.4142137);
    /// assert_approx_eq!(f32, vector.y, 3.4142137);
    /// assert_approx_eq!(f32, vector.z, 3.0);
    /// ```
    pub fn shift(&mut self, orientation: Vector3D, distance: f32) {
        let unit = orientation.to_unit();

        self.x += unit.x * distance;
        self.y += unit.y * distance;
        self.z += unit.z * distance;
    }

    /// Wrap coordinates of `Vector3D` so that each of them fits into the simulation box.
    /// Assumes orthogonal simulation box.
    ///
    /// ## Example
    /// ```
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let mut point = Vector3D::from([-0.5, 2.0, 4.2]);
    /// let simbox = SimBox::from([4.0, 4.0, 4.0]);
    ///
    /// point.wrap(&simbox);
    /// assert_approx_eq!(f32, point.x, 3.5, epsilon = 0.00001);
    /// assert_approx_eq!(f32, point.y, 2.0, epsilon = 0.00001);
    /// assert_approx_eq!(f32, point.z, 0.2, epsilon = 0.00001);
    /// ```
    pub fn wrap(&mut self, sbox: &SimBox) {
        self.x = Vector3D::wrap_coordinate(self.x, sbox.x);
        self.y = Vector3D::wrap_coordinate(self.y, sbox.y);
        self.z = Vector3D::wrap_coordinate(self.z, sbox.z);
    }

    ///Wrap a single coordinate into a simulation box.
    ///
    /// ## Note on performance
    /// You may think that the body of this function should look rather like this:
    /// ```ignore
    /// let wrapped = coor % box_len;
    /// if wrapped < 0.0 {
    ///    wrapped + box_len
    /// } else {
    ///    wrapped
    /// }
    /// ```
    /// But, taking a remainder of floats is actually a very costly operation. As in almost all cases,
    /// we will only perform AT MOST one iteration through one of the loops, the "loop" version of the wrapping
    /// function is usually much faster.
    /// This does not hold true, if the coordinate is VERY far away from the bounding box but that almost never happens.
    fn wrap_coordinate(coor: f32, box_len: f32) -> f32 {
        let mut wrapped = coor;

        while wrapped > box_len {
            wrapped -= box_len;
        }

        while wrapped < 0.0f32 {
            wrapped += box_len;
        }

        wrapped
    }

    /// Calculate distance between two points in the specified dimensions.
    /// In case a one-dimensional distance is requested, returns **oriented** distance
    /// between `self` and `point`.
    /// Takes periodic boundary conditions into consideration.
    ///
    /// ## Warning
    /// - Currently only works with orthogonal simulation boxes.
    ///
    /// ## Example 1: Two-dimensional distance
    /// Calculate distance between two points in the xy-plane.
    /// ```
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let point1 = Vector3D::from([1.0, 2.0, 3.0]);
    /// let point2 = Vector3D::from([3.5, 1.0, 2.0]);
    ///
    /// let simbox = SimBox::from([4.0, 4.0, 4.0]);
    ///
    /// let distance = point1.distance(&point2, Dimension::XY, &simbox);
    /// assert_approx_eq!(f32, distance, 1.802776);
    /// ```
    /// ## Example 2: One-dimensional oriented distance
    /// Calculate oriented distance between two points along the x-axis.
    /// ```
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let point1 = Vector3D::from([1.0, 2.0, 3.0]);
    /// let point2 = Vector3D::from([3.5, 1.0, 2.0]);
    ///
    /// let simbox = SimBox::from([4.0, 4.0, 4.0]);
    ///
    /// let distance1 = point1.distance(&point2, Dimension::X, &simbox);
    /// assert_approx_eq!(f32, distance1, 1.5);  // PBC are in play here, do not be confused by the sign
    /// let distance2 = point2.distance(&point1, Dimension::X, &simbox);
    /// assert_approx_eq!(f32, distance2, -1.5);
    /// ```
    pub fn distance(&self, point: &Vector3D, dim: Dimension, sbox: &SimBox) -> f32 {
        match dim {
            Dimension::None => 0.0,
            Dimension::X => Vector3D::min_image(self.x - point.x, sbox.x),
            Dimension::Y => Vector3D::min_image(self.y - point.y, sbox.y),
            Dimension::Z => Vector3D::min_image(self.z - point.z, sbox.z),
            Dimension::XY => {
                let dx = Vector3D::min_image(self.x - point.x, sbox.x);
                let dy = Vector3D::min_image(self.y - point.y, sbox.y);
                (dx * dx + dy * dy).sqrt()
            }
            Dimension::XZ => {
                let dx = Vector3D::min_image(self.x - point.x, sbox.x);
                let dz = Vector3D::min_image(self.z - point.z, sbox.z);
                (dx * dx + dz * dz).sqrt()
            }
            Dimension::YZ => {
                let dy = Vector3D::min_image(self.y - point.y, sbox.y);
                let dz = Vector3D::min_image(self.z - point.z, sbox.z);
                (dy * dy + dz * dz).sqrt()
            }
            Dimension::XYZ => {
                let dx = Vector3D::min_image(self.x - point.x, sbox.x);
                let dy = Vector3D::min_image(self.y - point.y, sbox.y);
                let dz = Vector3D::min_image(self.z - point.z, sbox.z);
                (dx * dx + dy * dy + dz * dz).sqrt()
            }
        }
    }

    /// Calculate distance between two points in the specified dimensions. **Ignore PBC.**
    /// In case a one-dimensional distance is requested, returns naive **oriented** distance
    /// between `self` and `point`.
    ///
    /// ## Example 1: Two-dimensional distance
    /// Calculate naive distance between two points in the xy-plane.
    /// ```
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let point1 = Vector3D::from([1.0, 2.0, 3.0]);
    /// let point2 = Vector3D::from([3.5, 1.0, 2.0]);
    ///
    /// let distance = point1.distance_naive(&point2, Dimension::XY);
    /// assert_approx_eq!(f32, distance, 2.692582);
    /// ```
    /// ## Example 2: One-dimensional oriented distance
    /// Calculate naive oriented distance between two points along the x-axis.
    /// ```
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let point1 = Vector3D::from([1.0, 2.0, 3.0]);
    /// let point2 = Vector3D::from([3.5, 1.0, 2.0]);
    ///
    /// let distance1 = point1.distance_naive(&point2, Dimension::X);
    /// assert_approx_eq!(f32, distance1, -2.5);
    /// let distance2 = point2.distance_naive(&point1, Dimension::X);
    /// assert_approx_eq!(f32, distance2, 2.5);
    /// ```
    pub fn distance_naive(&self, point: &Vector3D, dim: Dimension) -> f32 {
        match dim {
            Dimension::None => 0.0,
            Dimension::X => self.x - point.x,
            Dimension::Y => self.y - point.y,
            Dimension::Z => self.z - point.z,
            Dimension::XY => {
                let dx = self.x - point.x;
                let dy = self.y - point.y;
                (dx * dx + dy * dy).sqrt()
            }
            Dimension::XZ => {
                let dx = self.x - point.x;
                let dz = self.z - point.z;
                (dx * dx + dz * dz).sqrt()
            }
            Dimension::YZ => {
                let dy = self.y - point.y;
                let dz = self.z - point.z;
                (dy * dy + dz * dz).sqrt()
            }
            Dimension::XYZ => {
                let dx = self.x - point.x;
                let dy = self.y - point.y;
                let dz = self.z - point.z;
                (dx * dx + dy * dy + dz * dz).sqrt()
            }
        }
    }

    /// Calculate shortest vector connecting `self` with `point` taking PBC into consideration.
    ///
    /// ## Warning
    /// Currently only works with rectangular simulation boxes.
    ///
    /// ## Notes
    /// - In case the real point and its periodic image are equidistant from `self`,
    /// this function will return either of the two shortest vectors. It is not defined
    /// which vector will be returned.
    ///
    /// ## Example
    /// ```
    /// use groan_rs::prelude::*;
    /// use float_cmp::assert_approx_eq;
    ///
    /// let point1 = Vector3D::from([1.0, 2.0, 3.0]);
    /// let point2 = Vector3D::from([3.0, 2.0, 1.0]);
    /// let simbox = SimBox::from([3.5, 5.0, 5.0]);
    ///
    /// let vec = point1.vector_to(&point2, &simbox);
    ///
    /// assert_approx_eq!(f32, vec.x, -1.5);
    /// assert_approx_eq!(f32, vec.y, 0.0);
    /// assert_approx_eq!(f32, vec.z, -2.0);
    /// ```
    pub fn vector_to(&self, point: &Vector3D, sbox: &SimBox) -> Vector3D {
        let box_half = Vector3D {
            x: sbox.x / 2.0,
            y: sbox.y / 2.0,
            z: sbox.z / 2.0,
        };

        Vector3D::from([
            floor_mod(point.x - self.x + box_half.x, sbox.x) - box_half.x,
            floor_mod(point.y - self.y + box_half.y, sbox.y) - box_half.y,
            floor_mod(point.z - self.z + box_half.z, sbox.z) - box_half.z,
        ])
    }

    /// Takes the distance and returns a new distance that is modified according to the minimum image convention.
    fn min_image(dx: f32, box_len: f32) -> f32 {
        let half_box = box_len / 2.0;
        let mut new_dx = dx;

        while new_dx > half_box {
            new_dx -= box_len;
        }

        while new_dx < -half_box {
            new_dx += box_len;
        }

        new_dx
    }

    /// Apply `Dimension` as a filter for `Vector3D`.
    ///
    /// All dimensions of the `Vector3D` which do not match the `Dimension` are set to 0.
    ///
    /// ## Example
    /// ```
    /// use groan_rs::prelude::*;
    ///
    /// let mut vec: Vector3D = [1.0, 2.0, 3.0].into();
    ///
    /// vec.filter(Dimension::XZ);
    ///
    /// assert_eq!(vec.x, 1.0);
    /// assert_eq!(vec.y, 0.0);
    /// assert_eq!(vec.z, 3.0);
    /// ```
    pub fn filter(&mut self, dim: Dimension) {
        if !dim.is_x() {
            self.x = 0.0;
        }

        if !dim.is_y() {
            self.y = 0.0;
        }

        if !dim.is_z() {
            self.z = 0.0;
        }
    }
}

impl Default for Vector3D {
    /// Create a zero vector.
    fn default() -> Self {
        Vector3D {
            x: 0.0f32,
            y: 0.0f32,
            z: 0.0f32,
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

    #[test]
    fn len() {
        let vec = Vector3D::from([4.3, 5.6, 1.2]);
        assert_approx_eq!(f32, vec.len(), 7.161704);
    }

    #[test]
    fn len_null() {
        let vec = Vector3D::from([0.0, 0.0, 0.0]);
        assert_approx_eq!(f32, vec.len(), 0.0);
    }

    #[test]
    fn to_unit() {
        let vec = Vector3D::from([4.3, 5.6, 1.2]).to_unit();

        assert_approx_eq!(f32, vec.x, 0.6004158);
        assert_approx_eq!(f32, vec.y, 0.7819368);
        assert_approx_eq!(f32, vec.z, 0.16755791);
        assert_approx_eq!(f32, vec.len(), 1.0);
    }

    #[test]
    fn to_unit_null() {
        let vec = Vector3D::from([0.0, 0.0, 0.0]).to_unit();

        assert_approx_eq!(f32, vec.x, 0.0);
        assert_approx_eq!(f32, vec.y, 0.0);
        assert_approx_eq!(f32, vec.z, 0.0);
        assert_approx_eq!(f32, vec.len(), 0.0);
    }

    #[test]
    fn to_unit_small() {
        let vec = Vector3D::from([0.13, 0.0, 0.0]).to_unit();

        assert_approx_eq!(f32, vec.x, 1.0);
        assert_approx_eq!(f32, vec.y, 0.0);
        assert_approx_eq!(f32, vec.z, 0.0);
        assert_approx_eq!(f32, vec.len(), 1.0);
    }

    #[test]
    fn dot_1() {
        let vector1 = Vector3D::from([1.0, 0.0, 0.0]);
        let vector2 = Vector3D::from([0.0, 1.0, 1.0]);

        assert_approx_eq!(f32, vector1.dot(&vector2), 0.0);
    }

    #[test]
    fn dot_2() {
        let vector1 = Vector3D::from([2.0, 3.0, 4.0]);
        let vector2 = Vector3D::from([1.0, 2.0, 3.0]);

        assert_approx_eq!(f32, vector1.dot(&vector2), 20.0);
    }

    #[test]
    fn dot_3() {
        let vector1 = Vector3D::from([-2.0, 0.0, 5.0]);
        let vector2 = Vector3D::from([3.0, 1.0, -4.0]);

        assert_approx_eq!(f32, vector1.dot(&vector2), -26.0);
    }

    #[test]
    fn dot_4() {
        let vector1 = Vector3D::from([-2.0, 0.0, 5.0]);
        let vector2 = Vector3D::from([-3.0, 1.0, -4.0]);

        assert_approx_eq!(f32, vector1.dot(&vector2), -14.0);
    }

    #[test]
    fn dot_5() {
        let vector1 = Vector3D::from([-2.5, 0.3, 5.1]);
        let vector2 = Vector3D::from([-3.9, 1.1, -4.2]);

        assert_approx_eq!(f32, vector1.dot(&vector2), -11.34);
    }

    #[test]
    fn cross_1() {
        let vector1 = Vector3D::from([1.0, 0.0, 0.0]);
        let vector2 = Vector3D::from([0.0, 1.0, 1.0]);

        let cross = vector1.cross(&vector2);
        assert_approx_eq!(f32, cross.x, 0.0);
        assert_approx_eq!(f32, cross.y, -1.0);
        assert_approx_eq!(f32, cross.z, 1.0);
    }

    #[test]
    fn cross_2() {
        let vector1 = Vector3D::from([2.0, 3.0, 4.0]);
        let vector2 = Vector3D::from([1.0, 2.0, 3.0]);

        let cross = vector1.cross(&vector2);
        assert_approx_eq!(f32, cross.x, 1.0);
        assert_approx_eq!(f32, cross.y, -2.0);
        assert_approx_eq!(f32, cross.z, 1.0);
    }

    #[test]
    fn cross_3() {
        let vector1 = Vector3D::from([-2.0, 0.0, 5.0]);
        let vector2 = Vector3D::from([3.0, 1.0, -4.0]);

        let cross = vector1.cross(&vector2);
        assert_approx_eq!(f32, cross.x, -5.0);
        assert_approx_eq!(f32, cross.y, 7.0);
        assert_approx_eq!(f32, cross.z, -2.0);
    }

    #[test]
    fn cross_4() {
        let vector1 = Vector3D::from([-2.0, 0.0, 5.0]);
        let vector2 = Vector3D::from([-3.0, 1.0, -4.0]);

        let cross = vector1.cross(&vector2);
        assert_approx_eq!(f32, cross.x, -5.0);
        assert_approx_eq!(f32, cross.y, -23.0);
        assert_approx_eq!(f32, cross.z, -2.0);
    }

    #[test]
    fn cross_5() {
        let vector1 = Vector3D::from([-2.5, 0.3, 5.1]);
        let vector2 = Vector3D::from([-3.9, 1.1, -4.2]);

        let cross = vector1.cross(&vector2);
        assert_approx_eq!(f32, cross.x, -6.87);
        assert_approx_eq!(f32, cross.y, -30.39);
        assert_approx_eq!(f32, cross.z, -1.58);
    }

    #[test]
    fn angle_1() {
        let vector1 = Vector3D::from([2.0, 0.0, 0.0]);
        let vector2 = Vector3D::from([0.0, 2.0, 0.0]);

        assert_approx_eq!(f32, vector1.angle(&vector2), 1.57079632);
    }

    #[test]
    fn angle_2() {
        let vector1 = Vector3D::from([2.0, 0.0, 0.0]);
        let vector2 = Vector3D::from([0.0, -2.0, 0.0]);

        assert_approx_eq!(f32, vector1.angle(&vector2), 1.57079632);
    }

    #[test]
    fn angle_3() {
        let vector1 = Vector3D::from([1.0, 0.0, 0.0]);
        let vector2 = Vector3D::from([0.0, 0.0, 7.0]);

        assert_approx_eq!(f32, vector1.angle(&vector2), 1.57079632);
    }

    #[test]
    fn angle_4() {
        let vector1 = Vector3D::from([1.0, 0.0, 0.0]);
        let vector2 = Vector3D::from([3.0, 0.0, 3.0]);

        assert_approx_eq!(f32, vector1.angle(&vector2), 0.785398163);
    }

    #[test]
    fn angle_5() {
        let vector1 = Vector3D::from([1.0, 0.0, 0.0]);
        let vector2 = Vector3D::from([4.0, 0.0, 0.0]);

        assert_approx_eq!(f32, vector1.angle(&vector2), 0.0);
    }

    #[test]
    fn angle_6() {
        let vector1 = Vector3D::from([1.0, 0.0, 0.0]);
        let vector2 = Vector3D::from([-4.0, 0.0, 0.0]);

        assert_approx_eq!(f32, vector1.angle(&vector2), 3.14159265);
    }

    #[test]
    fn angle_7() {
        let vector1 = Vector3D::from([1.0, -1.0, 3.5]);
        let vector2 = Vector3D::from([1.2, 2.4, -0.7]);

        assert_approx_eq!(f32, vector1.angle(&vector2), 1.9269546);
    }

    #[test]
    fn shift_x() {
        let mut vector = Vector3D::from([-2.5, 0.3, 5.1]);
        let orientation = Vector3D::from(Dimension::X);

        vector.shift(orientation, 1.5);

        assert_approx_eq!(f32, vector.x, -1.0);
        assert_approx_eq!(f32, vector.y, 0.3);
        assert_approx_eq!(f32, vector.z, 5.1);
    }

    #[test]
    fn shift_xyz() {
        let mut vector = Vector3D::from([-2.5, 0.3, 5.1]);
        let orientation = Vector3D::from(Dimension::XYZ);

        vector.shift(orientation, 3.5);

        assert_approx_eq!(f32, vector.x, -0.479274);
        assert_approx_eq!(f32, vector.y, 2.320726);
        assert_approx_eq!(f32, vector.z, 7.120726);

        assert_approx_eq!(
            f32,
            vector.distance_naive(&[-2.5, 0.3, 5.1].into(), Dimension::XYZ),
            3.5
        );
    }

    #[test]
    fn shift_arbitrary() {
        let mut vector = Vector3D::from([-2.5, 0.3, 5.1]);
        let orientation = Vector3D::from([1.0, 0.5, 2.0]);

        vector.shift(orientation, 4.2);

        assert_approx_eq!(f32, vector.x, -0.666970); // shifted by 1.833
        assert_approx_eq!(f32, vector.y, 1.216515); // shifted by 0.917
        assert_approx_eq!(f32, vector.z, 8.766060); // shifted by 3.666

        assert_approx_eq!(
            f32,
            vector.distance_naive(&[-2.5, 0.3, 5.1].into(), Dimension::XYZ),
            4.2
        );
    }

    #[test]
    fn shift_arbitrary_negative() {
        let mut vector = Vector3D::from([-2.5, 0.3, 5.1]);
        let orientation = Vector3D::from([1.0, 0.5, 2.0]);

        vector.shift(orientation, -4.2);

        assert_approx_eq!(f32, vector.x, -4.333030); // shifted by -1.833
        assert_approx_eq!(f32, vector.y, -0.616515); // shifted by -0.917
        assert_approx_eq!(f32, vector.z, 1.433940); // shifted by -3.666

        assert_approx_eq!(
            f32,
            vector.distance_naive(&[-2.5, 0.3, 5.1].into(), Dimension::XYZ),
            4.2
        );
    }

    #[test]
    fn wrap() {
        let mut vector1 = Vector3D::from([-1.0, 1.5, 3.0]);
        let mut vector2 = Vector3D::from([2.0, 2.2, -0.3]);
        let mut vector3 = Vector3D::from([-54.2, 77.8, 124.5]);
        let simbox = SimBox::from([2.0, 2.0, 2.0]);

        vector1.wrap(&simbox);
        assert_approx_eq!(f32, vector1.x, 1.0);
        assert_approx_eq!(f32, vector1.y, 1.5);
        assert_approx_eq!(f32, vector1.z, 1.0);

        vector2.wrap(&simbox);
        assert_approx_eq!(f32, vector2.x, 2.0);
        assert_approx_eq!(f32, vector2.y, 0.2);
        assert_approx_eq!(f32, vector2.z, 1.7);

        vector3.wrap(&simbox);
        assert_approx_eq!(f32, vector3.x, 1.8, epsilon = 0.00001);
        assert_approx_eq!(f32, vector3.y, 1.8, epsilon = 0.00001);
        assert_approx_eq!(f32, vector3.z, 0.5, epsilon = 0.00001);
    }

    #[test]
    fn distance_x() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            point1.distance(&point2, Dimension::X, &simbox),
            1.5,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance(&point1, Dimension::X, &simbox),
            -1.5,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_y() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            point1.distance(&point2, Dimension::Y, &simbox),
            -0.2,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance(&point1, Dimension::Y, &simbox),
            0.2,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_z() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            point1.distance(&point2, Dimension::Z, &simbox),
            -1.8,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance(&point1, Dimension::Z, &simbox),
            1.8,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_xy() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            point1.distance(&point2, Dimension::XY, &simbox),
            1.51327,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance(&point1, Dimension::XY, &simbox),
            1.51327,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_xz() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            point1.distance(&point2, Dimension::XZ, &simbox),
            2.34307,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance(&point1, Dimension::XZ, &simbox),
            2.34307,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_yz() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            point1.distance(&point2, Dimension::YZ, &simbox),
            1.81108,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance(&point1, Dimension::YZ, &simbox),
            1.81108,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_xyz() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(
            f32,
            point1.distance(&point2, Dimension::XYZ, &simbox),
            2.351595,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance(&point1, Dimension::XYZ, &simbox),
            2.351595,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_xyz_outofbox() {
        let point1 = Vector3D::from([-1.0, 4.5, 2.3]);
        let point2 = Vector3D::from([3.5, -0.5, 4.2]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(f32, point1.distance(&point2, Dimension::X, &simbox), -0.5);
        assert_approx_eq!(f32, point1.distance(&point2, Dimension::Y, &simbox), 1.0);
        assert_approx_eq!(f32, point1.distance(&point2, Dimension::Z, &simbox), -1.9);
    }

    #[test]
    fn distance_none() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert_approx_eq!(f32, point1.distance(&point2, Dimension::None, &simbox), 0.0);
        assert_approx_eq!(f32, point2.distance(&point1, Dimension::None, &simbox), 0.0);
    }

    #[test]
    fn distance_naive_x() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        assert_approx_eq!(
            f32,
            point1.distance_naive(&point2, Dimension::X),
            -2.5,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance_naive(&point1, Dimension::X),
            2.5,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_y() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        assert_approx_eq!(
            f32,
            point1.distance_naive(&point2, Dimension::Y),
            3.8,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance_naive(&point1, Dimension::Y),
            -3.8,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_z() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        assert_approx_eq!(
            f32,
            point1.distance_naive(&point2, Dimension::Z),
            2.2,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance_naive(&point1, Dimension::Z),
            -2.2,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_xy() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        assert_approx_eq!(
            f32,
            point1.distance_naive(&point2, Dimension::XY),
            4.548626,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance_naive(&point1, Dimension::XY),
            4.548626,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_xz() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        assert_approx_eq!(
            f32,
            point1.distance_naive(&point2, Dimension::XZ),
            3.330165,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance_naive(&point1, Dimension::XZ),
            3.330165,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_yz() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        assert_approx_eq!(
            f32,
            point1.distance_naive(&point2, Dimension::YZ),
            4.390900,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance_naive(&point1, Dimension::YZ),
            4.390900,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_xyz() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        assert_approx_eq!(
            f32,
            point1.distance_naive(&point2, Dimension::XYZ),
            5.052722,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance_naive(&point1, Dimension::XYZ),
            5.052722,
            epsilon = 0.00001
        );
    }

    #[test]
    fn distance_naive_none() {
        let point1 = Vector3D::from([1.0, 3.9, 2.6]);
        let point2 = Vector3D::from([3.5, 0.1, 0.4]);

        assert_approx_eq!(
            f32,
            point1.distance_naive(&point2, Dimension::None),
            0.0,
            epsilon = 0.00001
        );
        assert_approx_eq!(
            f32,
            point2.distance_naive(&point1, Dimension::None),
            0.0,
            epsilon = 0.00001
        );
    }

    #[test]
    fn vector_to_nopbc() {
        let p1 = Vector3D::from([4.0, 4.0, 5.0]);
        let p2 = Vector3D::from([5.0, 5.0, 3.0]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        let vec = p1.vector_to(&p2, &simbox);

        assert_approx_eq!(f32, vec.x, 1.0);
        assert_approx_eq!(f32, vec.y, 1.0);
        assert_approx_eq!(f32, vec.z, -2.0);
    }

    #[test]
    fn vector_to_one_dim() {
        let p1 = Vector3D::from([3.0, 0.0, 7.0]);
        let p2 = Vector3D::from([1.0, 2.0, 1.0]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        let vec = p1.vector_to(&p2, &simbox);

        assert_approx_eq!(f32, vec.x, -2.0);
        assert_approx_eq!(f32, vec.y, 2.0);
        assert_approx_eq!(f32, vec.z, 4.0);
    }

    #[test]
    fn vector_to_two_dim() {
        let p1 = Vector3D::from([1.0, 2.0, 5.0]);
        let p2 = Vector3D::from([9.0, 8.0, 6.0]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        let vec = p1.vector_to(&p2, &simbox);

        assert_approx_eq!(f32, vec.x, -2.0);
        assert_approx_eq!(f32, vec.y, -4.0);
        assert_approx_eq!(f32, vec.z, 1.0);
    }

    #[test]
    fn vector_to_all_dim() {
        let p1 = Vector3D::from([8.0, 9.0, 2.0]);
        let p2 = Vector3D::from([1.0, 3.0, 9.0]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        let vec = p1.vector_to(&p2, &simbox);

        assert_approx_eq!(f32, vec.x, 3.0);
        assert_approx_eq!(f32, vec.y, 4.0);
        assert_approx_eq!(f32, vec.z, -3.0);
    }

    #[test]
    fn vector_to_same() {
        let p1 = Vector3D::from([0.0, 3.0, 10.0]);
        let p2 = Vector3D::from([10.0, 3.0, 0.0]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        let vec = p1.vector_to(&p2, &simbox);

        assert_approx_eq!(f32, vec.x, 0.0);
        assert_approx_eq!(f32, vec.y, 0.0);
        assert_approx_eq!(f32, vec.z, 0.0);
    }

    #[test]
    fn vector_to_equidistant() {
        let p1 = Vector3D::from([7.0, 4.0, 3.0]);
        let p2 = Vector3D::from([2.0, 5.0, 2.0]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        let vec = p1.vector_to(&p2, &simbox);

        assert_approx_eq!(f32, vec.x.abs(), 5.0);
        assert_approx_eq!(f32, vec.y, 1.0);
        assert_approx_eq!(f32, vec.z, -1.0);
    }

    #[test]
    fn filter_vec_x() {
        let mut vec: Vector3D = [4.3, 1.8, 2.7].into();
        vec.filter(Dimension::X);
        assert_eq!(vec, [4.3, 0.0, 0.0].into());
    }

    #[test]
    fn filter_vec_y() {
        let mut vec: Vector3D = [4.3, 1.8, 2.7].into();
        vec.filter(Dimension::Y);
        assert_eq!(vec, [0.0, 1.8, 0.0].into());
    }

    #[test]
    fn filter_vec_z() {
        let mut vec: Vector3D = [4.3, 1.8, 2.7].into();
        vec.filter(Dimension::Z);
        assert_eq!(vec, [0.0, 0.0, 2.7].into());
    }

    #[test]
    fn filter_vec_xy() {
        let mut vec: Vector3D = [4.3, 1.8, 2.7].into();
        vec.filter(Dimension::XY);
        assert_eq!(vec, [4.3, 1.8, 0.0].into());
    }

    #[test]
    fn filter_vec_xz() {
        let mut vec: Vector3D = [4.3, 1.8, 2.7].into();
        vec.filter(Dimension::XZ);
        assert_eq!(vec, [4.3, 0.0, 2.7].into());
    }

    #[test]
    fn filter_vec_yz() {
        let mut vec: Vector3D = [4.3, 1.8, 2.7].into();
        vec.filter(Dimension::YZ);
        assert_eq!(vec, [0.0, 1.8, 2.7].into());
    }

    #[test]
    fn filter_vec_xyz() {
        let mut vec: Vector3D = [4.3, 1.8, 2.7].into();
        vec.filter(Dimension::XYZ);
        assert_eq!(vec, [4.3, 1.8, 2.7].into());
    }

    #[test]
    fn filter_vec_none() {
        let mut vec: Vector3D = [4.3, 1.8, 2.7].into();
        vec.filter(Dimension::None);
        assert_eq!(vec, [0.0, 0.0, 0.0].into());
    }

    #[test]
    fn dim_to_vector_none() {
        let vec = Vector3D::from(Dimension::None);
        assert_eq!(vec, [0.0, 0.0, 0.0].into());
    }

    #[test]
    fn dim_to_vector_x() {
        let vec = Vector3D::from(Dimension::X);
        assert_eq!(vec, [1.0, 0.0, 0.0].into());
    }

    #[test]
    fn dim_to_vector_y() {
        let vec = Vector3D::from(Dimension::Y);
        assert_eq!(vec, [0.0, 1.0, 0.0].into());
    }

    #[test]
    fn dim_to_vector_z() {
        let vec = Vector3D::from(Dimension::Z);
        assert_eq!(vec, [0.0, 0.0, 1.0].into());
    }

    #[test]
    fn dim_to_vector_xy() {
        let vec = Vector3D::from(Dimension::XY);
        assert_eq!(vec, [REC_SQRT2, REC_SQRT2, 0.0].into());
        assert_approx_eq!(f32, vec.len(), 1.0);
    }

    #[test]
    fn dim_to_vector_xz() {
        let vec = Vector3D::from(Dimension::XZ);
        assert_eq!(vec, [REC_SQRT2, 0.0, REC_SQRT2].into());
        assert_approx_eq!(f32, vec.len(), 1.0);
    }

    #[test]
    fn dim_to_vector_yz() {
        let vec = Vector3D::from(Dimension::YZ);
        assert_eq!(vec, [0.0, REC_SQRT2, REC_SQRT2].into());
        assert_approx_eq!(f32, vec.len(), 1.0);
    }

    #[test]
    fn dim_to_vector_xyz() {
        let vec = Vector3D::from(Dimension::XYZ);
        assert_eq!(vec, [REC_SQRT3, REC_SQRT3, REC_SQRT3].into());
        assert_approx_eq!(f32, vec.len(), 1.0);
    }
}
