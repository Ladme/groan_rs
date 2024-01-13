// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of shapes for geometry selection.

use crate::structures::{dimension::Dimension, simbox::SimBox, vector3d::Vector3D};

/// Structure describing a sphere for geometry selections.
#[derive(Debug, Clone)]
pub struct Sphere {
    /// Coordinates of the center of the sphere.
    position: Vector3D,
    /// Radius of the sphere.
    radius: f32,
}

/// Structure describing a rectangular box for geometry selections.
#[derive(Debug, Clone)]
pub struct Rectangular {
    /// Coordinates of the box origin.
    position: Vector3D,
    /// Length of the box side along the x-axis.
    x: f32,
    /// Length of the box side along the y-axis.
    y: f32,
    /// Length of the box side along the z-axis.
    z: f32,
}

/// Structure describing a cylinder for geometry selections.
#[derive(Debug, Clone)]
pub struct Cylinder {
    /// Coordinates of the cylinder base center.
    position: Vector3D,
    /// Radius of the cylinder.
    radius: f32,
    /// Height of the cylinder.
    height: f32,
    /// Orientation of the cylinder in space.
    orientation: Dimension,
    /// Plane in which the cylinder base is placed.
    plane: Dimension,
}

/// Any structure implementing this trait can be used for geometry selection.
pub trait Shape {
    /// Returns `true` if target point is inside the `Shape`. Else returns `false`.
    fn inside(&self, point: &Vector3D, simbox: &SimBox) -> bool;
}

impl Sphere {
    /// Construct a new Sphere.
    ///
    /// ## Arguments
    /// - `position` - Coordinates of the center of the sphere.
    /// - `radius` - Radius of the sphere.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// // constructs a sphere with radius 2.0 nm
    /// let sphere = Sphere::new(
    ///     [1.0, 2.0, 3.0].into(), // position of the sphere center
    ///     2.0,                    // radius of the sphere
    /// );
    /// ```
    pub fn new(position: Vector3D, radius: f32) -> Self {
        Sphere { position, radius }
    }

    /// Get position of the center of the sphere.
    pub fn get_position(&self) -> &Vector3D {
        &self.position
    }

    /// Get radius of the sphere.
    pub fn get_radius(&self) -> f32 {
        self.radius
    }
}

impl Shape for Sphere {
    /// Check if point is inside the sphere.
    fn inside(&self, point: &Vector3D, simbox: &SimBox) -> bool {
        point.distance(&self.position, Dimension::XYZ, simbox) < self.radius
    }
}

impl Rectangular {
    /// Construct a new Rectangular box.
    ///
    /// ## Arguments
    /// - `position` - Coordinates of the box origin.
    /// - `x` - Length of the box side along the x-axis.
    /// - `y` - Length of the box side along the y-axis.
    /// - `z` - Length of the box side along the z-axis.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// // constructs a rectangular box of specified size
    /// let rect = Rectangular::new(
    ///     [1.0, 2.0, 3.0].into(), // position of box origin
    ///     1.0,    // size of the box along x-dimension
    ///     2.0,    // size of the box along y-dimension
    ///     3.0,    // size of the box along z-dimension
    /// );
    /// ```
    pub fn new(position: Vector3D, x: f32, y: f32, z: f32) -> Self {
        Rectangular { position, x, y, z }
    }

    /// Get coordinate of the box origin.
    pub fn get_position(&self) -> &Vector3D {
        &self.position
    }

    /// Get length of the box side along the x-axis.
    pub fn get_x(&self) -> f32 {
        self.x
    }

    /// Get length of the box side along the y-axis.
    pub fn get_y(&self) -> f32 {
        self.y
    }

    /// Get length of the box side along the z-axis.
    pub fn get_z(&self) -> f32 {
        self.z
    }
}

impl Shape for Rectangular {
    /// Check if point is inside the rectangular box.
    fn inside(&self, point: &Vector3D, simbox: &SimBox) -> bool {
        let mut dx = point.distance(&self.position, Dimension::X, simbox);
        if dx < 0.0 {
            dx += simbox.x;
        }
        let mut dy = point.distance(&self.position, Dimension::Y, simbox);
        if dy < 0.0 {
            dy += simbox.y;
        }
        let mut dz = point.distance(&self.position, Dimension::Z, simbox);
        if dz < 0.0 {
            dz += simbox.z;
        }

        dx <= self.x && dy <= self.y && dz <= self.z
    }
}

impl Cylinder {
    /// Construct a new Cylinder.
    ///
    /// ## Arguments
    /// - `position` - Coordinates of the center of the cylinder's base.
    /// - `orientation` - Orientation of the cylinder in space. Only X, Y, and Z are supported.
    /// - `radius` - Radius of the cylinder.
    /// - `height` - Height of the cylinder along the main axis of the cylinder.
    ///
    /// ## Panics
    /// Panics if `orientation` is not `Dimension::X`, `Dimension::Y`, nor `Dimension::Z`.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// // constructs a cylinder oriented along the z-dimension of the simulation box
    /// let cylinder = Cylinder::new(
    ///     [1.0, 2.0, 3.0].into(), // position
    ///     1.5,                    // radius
    ///     2.0,                    // height
    ///     Dimension::Z,           // orientation
    /// );
    /// ```
    pub fn new(position: Vector3D, radius: f32, height: f32, orientation: Dimension) -> Self {
        let plane = match orientation {
            Dimension::X => Dimension::YZ,
            Dimension::Y => Dimension::XZ,
            Dimension::Z => Dimension::XY,
            d => panic!(
                "FATAL GROAN ERROR | Cylinder::new | Unsupported orientation dimension '{}'.",
                d
            ),
        };

        Cylinder {
            position,
            radius,
            height,
            orientation,
            plane,
        }
    }

    /// Get position of the center of the cylinder's base.
    pub fn get_position(&self) -> &Vector3D {
        &self.position
    }

    /// Get radius of the cylinder.
    pub fn get_radius(&self) -> f32 {
        self.radius
    }

    /// Get height of the cylinder
    pub fn get_height(&self) -> f32 {
        self.height
    }

    /// Get orientation of the cylinder
    pub fn get_orientation(&self) -> Dimension {
        self.orientation
    }
}

impl Shape for Cylinder {
    /// Check if point is inside the cylinder.
    fn inside(&self, point: &Vector3D, simbox: &SimBox) -> bool {
        let mut distance_axis = point.distance(&self.position, self.orientation, simbox);

        if distance_axis < 0.0 {
            match self.orientation {
                Dimension::X => distance_axis += simbox.x,
                Dimension::Y => distance_axis += simbox.y,
                Dimension::Z => distance_axis += simbox.z,
                d => panic!("FATAL GROAN ERROR | Cylinder::inside | Orientation dimension '{}' should never occur in a cylinder shape.", d),
            }
        }

        if distance_axis > self.height
            || point.distance(&self.position, self.plane, simbox) > self.radius
        {
            return false;
        }

        true
    }
}

#[cfg(test)]
mod tests_sphere {
    use super::*;
    use float_cmp::assert_approx_eq;
    use rand::Rng;

    #[test]
    fn new() {
        let sphere = Sphere::new([1.0, 2.0, 3.0].into(), 1.5);

        assert_approx_eq!(f32, sphere.get_position().x, 1.0);
        assert_approx_eq!(f32, sphere.get_position().y, 2.0);
        assert_approx_eq!(f32, sphere.get_position().z, 3.0);
        assert_approx_eq!(f32, sphere.get_radius(), 1.5);
    }

    #[test]
    fn inside_nopbc() {
        let sphere = Sphere::new([1.0, 2.0, 3.0].into(), 1.5);

        let point = Vector3D::from([2.0, 2.5, 2.4]);
        let simbox = SimBox::from([5.0, 5.0, 5.0]);

        assert!(sphere.inside(&point, &simbox));
    }

    #[test]
    fn inside_pbc() {
        let sphere = Sphere::new([1.0, 2.0, 4.5].into(), 1.5);

        let point = Vector3D::from([4.8, 2.1, 0.3]);
        let simbox = SimBox::from([5.0, 5.0, 5.0]);

        assert!(sphere.inside(&point, &simbox));
    }

    #[test]
    fn not_inside() {
        let sphere = Sphere::new([1.0, 2.0, 4.5].into(), 1.5);

        let point = Vector3D::from([4.0, 2.1, 0.3]);
        let simbox = SimBox::from([5.0, 5.0, 5.0]);

        assert!(!sphere.inside(&point, &simbox));
    }

    #[test]
    fn inside_random() {
        let sphere_center = Vector3D::from([1.0, 2.0, 3.0]);
        let sphere_radius = 2.5;

        let sphere = Sphere::new(sphere_center.clone(), sphere_radius);
        let simbox = SimBox::from([5.0, 5.0, 5.0]);
        let mut rng = rand::thread_rng();

        for _ in 0..100 {
            let x = rng.gen_range(0.0..5.0);
            let y = rng.gen_range(0.0..5.0);
            let z = rng.gen_range(0.0..5.0);

            let point = Vector3D::from([x, y, z]);

            assert_eq!(
                sphere.inside(&point, &simbox),
                point.distance(&sphere_center, Dimension::XYZ, &simbox) < sphere_radius
            );
        }
    }
}

#[cfg(test)]
mod tests_rectangular {
    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn new() {
        let rect = Rectangular::new([1.0, 2.0, 3.0].into(), 2.0, 1.8, 0.7);

        assert_approx_eq!(f32, rect.get_position().x, 1.0);
        assert_approx_eq!(f32, rect.get_position().y, 2.0);
        assert_approx_eq!(f32, rect.get_position().z, 3.0);

        assert_approx_eq!(f32, rect.get_x(), 2.0);
        assert_approx_eq!(f32, rect.get_y(), 1.8);
        assert_approx_eq!(f32, rect.get_z(), 0.7);
    }

    #[test]
    fn inside_nopbc_1() {
        let rect = Rectangular::new([1.0, 2.0, 3.0].into(), 3.0, 2.0, 1.0);

        let point = Vector3D::from([3.1, 3.8, 3.9]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(rect.inside(&point, &simbox));
    }

    #[test]
    fn inside_nopbc_2() {
        let rect = Rectangular::new([2.5, 3.1, 0.3].into(), 1.2, 1.3, 5.0);

        let point = Vector3D::from([2.6, 4.3, 4.9]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(rect.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_nopbc_1() {
        let rect = Rectangular::new([1.0, 2.0, 3.0].into(), 3.0, 2.0, 1.0);

        let point = Vector3D::from([4.1, 3.8, 3.9]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!rect.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_nopbc_2() {
        let rect = Rectangular::new([1.0, 2.0, 3.0].into(), 3.0, 2.0, 1.0);

        let point = Vector3D::from([2.1, 1.9, 3.9]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!rect.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_nopbc_3() {
        let rect = Rectangular::new([1.0, 2.0, 3.0].into(), 3.0, 2.0, 1.0);

        let point = Vector3D::from([2.1, 2.5, 4.1]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!rect.inside(&point, &simbox));
    }

    #[test]
    fn inside_pbc_1() {
        let rect = Rectangular::new([1.0, 2.0, 3.0].into(), 4.0, 2.0, 1.5);

        let point = Vector3D::from([0.5, 3.8, 3.3]);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(rect.inside(&point, &simbox));
    }

    #[test]
    fn inside_pbc_2() {
        let rect = Rectangular::new([1.0, 2.0, 3.0].into(), 1.0, 4.0, 1.5);

        let point = Vector3D::from([1.3, 1.2, 3.5]);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(rect.inside(&point, &simbox));
    }

    #[test]
    fn inside_pbc_3() {
        let rect = Rectangular::new([1.0, 2.0, 3.0].into(), 1.0, 2.0, 1.5);

        let point = Vector3D::from([1.9, 2.2, 0.0]);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(rect.inside(&point, &simbox));
    }
}

#[cfg(test)]
mod tests_cylinder {
    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn new() {
        let cyl = Cylinder::new([1.0, 2.0, 3.0].into(), 0.8, 4.3, Dimension::Y);

        assert_approx_eq!(f32, cyl.get_position().x, 1.0);
        assert_approx_eq!(f32, cyl.get_position().y, 2.0);
        assert_approx_eq!(f32, cyl.get_position().z, 3.0);
        assert_eq!(cyl.get_orientation(), Dimension::Y);
        assert_approx_eq!(f32, cyl.get_radius(), 0.8);
        assert_approx_eq!(f32, cyl.get_height(), 4.3);
    }

    #[test]
    fn inside_x_nopbc() {
        let cylinder = Cylinder::new([2.0, 1.0, 3.0].into(), 2.0, 4.0, Dimension::X);

        let point = Vector3D::from([4.2, 1.8, 2.2]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_x_nopbc() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 2.0, 4.0, Dimension::X);

        let point = Vector3D::from([2.9, 3.8, 2.2]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!cylinder.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_x_nopbc2() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 2.0, 4.0, Dimension::X);

        let point = Vector3D::from([3.1, 4.6, 1.2]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_x_pbc_1() {
        let cylinder = Cylinder::new([2.0, 1.0, 3.0].into(), 2.0, 3.0, Dimension::X);

        let point = Vector3D::from([0.3, 1.4, 2.2]);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_x_pbc_2() {
        let cylinder = Cylinder::new([2.0, 1.0, 3.0].into(), 2.0, 3.0, Dimension::X);

        let point = Vector3D::from([2.4, 3.8, 2.8]);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_y_nopbc() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 4.0, 4.0, Dimension::Y);

        let point = Vector3D::from([5.2, 3.8, 3.2]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_y_nopbc() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 2.0, 4.0, Dimension::Y);

        let point = Vector3D::from([4.2, 7.3, 2.2]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!cylinder.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_y_nopbc2() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 2.0, 4.0, Dimension::Y);

        let point = Vector3D::from([1.1, 5.4, 3.7]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_y_pbc_1() {
        let cylinder = Cylinder::new([1.0, 2.0, 3.0].into(), 2.0, 3.0, Dimension::Y);

        let point = Vector3D::from([0.7, 0.8, 3.4]);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_y_pbc_2() {
        let cylinder = Cylinder::new([1.0, 2.0, 3.5].into(), 2.0, 3.0, Dimension::Y);

        let point = Vector3D::from([3.5, 3.4, 0.0]);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_z_nopbc() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 2.0, 2.0, Dimension::Z);

        let point = Vector3D::from([4.0, 3.8, 4.8]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_z_nopbc() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 2.0, 4.0, Dimension::Z);

        let point = Vector3D::from([4.0, 3.8, 7.2]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!cylinder.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_z_nopbc2() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 2.0, 4.0, Dimension::Z);

        let point = Vector3D::from([4.9, 3.7, 6.8]);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_z_pbc_1() {
        let cylinder = Cylinder::new([1.0, 2.0, 3.0].into(), 2.0, 3.0, Dimension::Z);

        let point = Vector3D::from([1.4, 1.5, 1.5]);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_z_pbc_2() {
        let cylinder = Cylinder::new([3.0, 2.0, 3.0].into(), 2.0, 3.0, Dimension::Z);

        let point = Vector3D::from([0.3, 1.5, 3.1]);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(cylinder.inside(&point, &simbox));
    }
}
