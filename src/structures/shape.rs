// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Implementation of shapes for geometry selection.

use crate::structures::{dimension::Dimension, simbox::SimBox, vector3d::Vector3D};

/// Structure describing a sphere for geometry selections.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(deny_unknown_fields))]
pub struct Sphere {
    /// Coordinates of the center of the sphere.
    position: Vector3D,
    /// Radius of the sphere.
    radius: f32,
}

/// Structure describing a rectangular box for geometry selections.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(deny_unknown_fields))]
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
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(deny_unknown_fields))]
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

/// Structure describing a triangular prism for geometry selections.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(deny_unknown_fields))]
pub struct TriangularPrism {
    /// Coordinates of the vertices of the base.
    base1: Vector3D,
    base2: Vector3D,
    base3: Vector3D,
    /// Height of the prism.
    height: f32,
    /// Orientation of the triangular prism in space.
    orientation: Dimension,
    /// Plane in which the triangular prism base is placed.
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
    ///     Vector3D::new(1.0, 2.0, 3.0), // position of the sphere center
    ///     2.0,                           // radius of the sphere
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
    /// Check if point is inside the sphere. Takes PBC into consideration.
    ///
    /// For PBC-unsafe variant, see [Sphere::inside_naive].
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
    ///     Vector3D::new(1.0, 2.0, 3.0), // position of box origin
    ///     1.0,                          // size of the box along x-dimension
    ///     2.0,                          // size of the box along y-dimension
    ///     3.0,                          // size of the box along z-dimension
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
    /// Check if point is inside the rectangular box. Takes PBC into consideration.
    ///
    /// For PBC-unsafe variant, see [Rectangular::inside_naive].
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
    ///     Vector3D::new(1.0, 2.0, 3.0), // position
    ///     1.5,                          // radius
    ///     2.0,                          // height
    ///     Dimension::Z,                 // orientation
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
    /// Check if point is inside the cylinder. Takes PBC into consideration.
    ///
    /// For PBC-unsafe variant, see [Cylinder::inside_naive].
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

impl TriangularPrism {
    /// Construct a new Triangular Prism.
    ///
    /// ## Arguments
    /// - `base1`, `base2`, `base3` - Coordinates of the vertices of the prism's base.
    /// - `height` - Perpendicular distance between the two bases of the prism.
    ///
    /// ## Panics
    /// Panics if the base of the prism is *not* entirely localized in the XY, XZ, or YZ plane.
    /// Panics if the base can't be constructed.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// // constructs a triangular prism oriented along the x-dimension of the simulation box
    /// let prism = TriangularPrism::new(
    ///     Vector3D::new(1.0, 2.0, 3.0),  // first vertex of the base
    ///     Vector3D::new(1.0, 4.5, 5.3),  // second vertex of the base
    ///     Vector3D::new(1.0, 2.3, 1.5),  // third vertex of the base
    ///     4.2                            // height of the prism
    /// );
    /// ```
    ///
    /// ## Notes
    /// - Periodic boundary conditions are not applied to the base of the triangular prism.
    ///   In other words, if you define the base of a triangular prism with a point
    ///   that is outside of the box boundaries, you will only select atoms located in the area
    ///   enclosed by the base itself, not by the periodic image of the base.
    ///
    /// ```text
    ///     ╔════════════════╗
    ///     ║                ║
    ///     ║       A        ║
    ///     ║    ###         ║
    ///     ║#######         ║
    ///  C  ║#####B       C'!║
    ///     ║                ║
    ///     ║                ║
    ///     ╚════════════════╝
    ///
    /// ```
    ///   `A`, `B`, `C` specify vertices of the base, `#` roughly specifies areas that are considered
    ///   to be inside the base. `C'` specifies periodic image of `C`.
    ///   Areas marked with `!` are *NOT* included in the prism, even though they could be if we considered PBC.
    ///
    /// - The above-described behavior may change in future versions of the `groan_rs` library!
    ///
    /// - Since the upper base of the prism is defined relative to the bottom base (using the `height` attribute),
    ///   periodic boundary conditions *ARE* used in this case. In other words,
    ///   if you define a prism with a height that would position the upper base outside of the box,
    ///   the atoms will be selected based on the periodic image of the upper base.
    ///
    /// ```text
    ///     ╔════════════════╗
    ///     ║      #######   ║
    ///     ║      #######   ║
    ///     ║      A#####B   ║
    ///     ║                ║
    ///     ║                ║
    ///     ║      #######   ║
    ///     ║      #######   ║
    ///     ╚════════════════╝
    ///
    /// ```
    ///   (The vertex `C` is not seen.)
    pub fn new(base1: Vector3D, base2: Vector3D, base3: Vector3D, height: f32) -> Self {
        let bases = [
            (base1.x, base2.x, base3.x, Dimension::X, Dimension::YZ),
            (base1.y, base2.y, base3.y, Dimension::Y, Dimension::XZ),
            (base1.z, base2.z, base3.z, Dimension::Z, Dimension::XY),
        ];

        let mut orientation = None;
        let mut plane = Dimension::None;

        for (val1, val2, val3, current_orientation, current_plane) in bases.iter() {
            if val1 == val2 && val2 == val3 {
                if orientation.is_some() {
                    panic!("FATAL GROAN ERROR | TriangularPrism::new | Base of the requested TriangularPrism can not be constructed.");
                }
                orientation = Some(*current_orientation);
                plane = *current_plane;
            }
        }

        if orientation.is_none() {
            panic!("FATAL GROAN ERROR | TriangularPrism::new | Base of the requested TriangularPrism does not lie in xy, xz, nor yz plane.");
        }

        // determine position of the bases

        TriangularPrism {
            base1,
            base2,
            base3,
            height,
            orientation: orientation.unwrap(),
            plane,
        }
    }

    /// Get coordinates of the first vertex of the base.
    pub fn get_base1(&self) -> &Vector3D {
        &self.base1
    }

    /// Get coordinates of the second vertex of the base.
    pub fn get_base2(&self) -> &Vector3D {
        &self.base2
    }

    /// Get coordinates of the third vertex of the base.
    pub fn get_base3(&self) -> &Vector3D {
        &self.base3
    }

    /// Get height of the triangular prism.
    pub fn get_height(&self) -> f32 {
        self.height
    }

    /// Get orientation of the triangular prism.
    pub fn get_orientation(&self) -> Dimension {
        self.orientation
    }

    /// Get plane of the base of the triangular prism.
    pub fn get_plane(&self) -> Dimension {
        self.plane
    }

    fn sign(point1: &Vector3D, point2: &Vector3D, point3: &Vector3D, plane: Dimension) -> f32 {
        match plane {
            Dimension::XY => {
                (point1.x - point3.x) * (point2.y - point3.y)
                    - (point2.x - point3.x) * (point1.y - point3.y)
            }
            Dimension::XZ => {
                (point1.x - point3.x) * (point2.z - point3.z)
                    - (point2.x - point3.x) * (point1.z - point3.z)
            }
            Dimension::YZ => {
                (point1.y - point3.y) * (point2.z - point3.z)
                    - (point2.y - point3.y) * (point1.z - point3.z)
            }
            _ => panic!(
                "FATAL GROAN ERROR | TriangularPrism::sign | This dimension should never occur."
            ),
        }
    }
}

impl Shape for TriangularPrism {
    /// Check if point is inside the triangular prism.
    ///
    /// Adapted from `https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle`
    fn inside(&self, point: &Vector3D, simbox: &SimBox) -> bool {
        let mut distance_from_base =
            point.distance(self.get_base1(), self.get_orientation(), simbox);

        if distance_from_base < 0.0 {
            match self.get_orientation() {
                Dimension::X => distance_from_base += simbox.x,
                Dimension::Y => distance_from_base += simbox.y,
                Dimension::Z => distance_from_base += simbox.z,
                d => panic!("FATAL GROAN ERROR | TriangularPrism::inside | Orientation dimension '{}' should never occur in a triangular prism.", d),
            }
        }

        if distance_from_base >= self.get_height() {
            return false;
        }

        let d1 = TriangularPrism::sign(point, self.get_base1(), self.get_base2(), self.get_plane());
        let d2 = TriangularPrism::sign(point, self.get_base2(), self.get_base3(), self.get_plane());
        let d3 = TriangularPrism::sign(point, self.get_base3(), self.get_base1(), self.get_plane());

        let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
        let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);

        !(has_neg && has_pos)
    }
}

/// Trait implemented by structures that can be used for geometry selection without PBC.
pub trait NaiveShape {
    /// Returns `true` if target point is inside the `NaiveShape`. **Does not take periodic boundary conditions into consideration.**
    fn inside_naive(&self, point: &Vector3D) -> bool;
}

impl NaiveShape for Sphere {
    /// Check if point is inside the sphere. **Does not take PBC into consideration.**
    ///
    /// For PBC-safe variant, see [Sphere::inside].
    fn inside_naive(&self, point: &Vector3D) -> bool {
        point.distance_naive(&self.position, Dimension::XYZ) < self.radius
    }
}

impl NaiveShape for Cylinder {
    /// Check if point is inside the cylinder. **Does not take PBC into consideration.**
    ///
    /// For PBC-safe variant, see [Cylinder::inside].
    fn inside_naive(&self, point: &Vector3D) -> bool {
        point.distance_naive(&self.position, self.orientation) < self.height
            && point.distance_naive(&self.position, self.plane) < self.radius
    }
}

impl NaiveShape for Rectangular {
    /// Check if point is inside the rectangular box.
    ///
    /// For PBC-safe variant, see [Rectangular::inside].
    fn inside_naive(&self, point: &Vector3D) -> bool {
        let dx = point.distance_naive(&self.position, Dimension::X);
        let dy = point.distance_naive(&self.position, Dimension::Y);
        let dz = point.distance_naive(&self.position, Dimension::Z);

        dx <= self.x && dy <= self.y && dz <= self.z
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

        let point = Vector3D::new(2.0, 2.5, 2.4);
        let simbox = SimBox::from([5.0, 5.0, 5.0]);

        assert!(sphere.inside(&point, &simbox));
    }

    #[test]
    fn inside_pbc() {
        let sphere = Sphere::new([1.0, 2.0, 4.5].into(), 1.5);

        let point = Vector3D::new(4.8, 2.1, 0.3);
        let simbox = SimBox::from([5.0, 5.0, 5.0]);

        assert!(sphere.inside(&point, &simbox));
    }

    #[test]
    fn not_inside() {
        let sphere = Sphere::new([1.0, 2.0, 4.5].into(), 1.5);

        let point = Vector3D::new(4.0, 2.1, 0.3);
        let simbox = SimBox::from([5.0, 5.0, 5.0]);

        assert!(!sphere.inside(&point, &simbox));
    }

    #[test]
    fn inside_random() {
        let sphere_center = Vector3D::new(1.0, 2.0, 3.0);
        let sphere_radius = 2.5;

        let sphere = Sphere::new(sphere_center.clone(), sphere_radius);
        let simbox = SimBox::from([5.0, 5.0, 5.0]);
        let mut rng = rand::thread_rng();

        for _ in 0..100 {
            let x = rng.gen_range(0.0..5.0);
            let y = rng.gen_range(0.0..5.0);
            let z = rng.gen_range(0.0..5.0);

            let point = Vector3D::new(x, y, z);

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

        let point = Vector3D::new(3.1, 3.8, 3.9);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(rect.inside(&point, &simbox));
    }

    #[test]
    fn inside_nopbc_2() {
        let rect = Rectangular::new([2.5, 3.1, 0.3].into(), 1.2, 1.3, 5.0);

        let point = Vector3D::new(2.6, 4.3, 4.9);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(rect.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_nopbc_1() {
        let rect = Rectangular::new([1.0, 2.0, 3.0].into(), 3.0, 2.0, 1.0);

        let point = Vector3D::new(4.1, 3.8, 3.9);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!rect.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_nopbc_2() {
        let rect = Rectangular::new([1.0, 2.0, 3.0].into(), 3.0, 2.0, 1.0);

        let point = Vector3D::new(2.1, 1.9, 3.9);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!rect.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_nopbc_3() {
        let rect = Rectangular::new([1.0, 2.0, 3.0].into(), 3.0, 2.0, 1.0);

        let point = Vector3D::new(2.1, 2.5, 4.1);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!rect.inside(&point, &simbox));
    }

    #[test]
    fn inside_pbc_1() {
        let rect = Rectangular::new([1.0, 2.0, 3.0].into(), 4.0, 2.0, 1.5);

        let point = Vector3D::new(0.5, 3.8, 3.3);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(rect.inside(&point, &simbox));
    }

    #[test]
    fn inside_pbc_2() {
        let rect = Rectangular::new([1.0, 2.0, 3.0].into(), 1.0, 4.0, 1.5);

        let point = Vector3D::new(1.3, 1.2, 3.5);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(rect.inside(&point, &simbox));
    }

    #[test]
    fn inside_pbc_3() {
        let rect = Rectangular::new([1.0, 2.0, 3.0].into(), 1.0, 2.0, 1.5);

        let point = Vector3D::new(1.9, 2.2, 0.0);
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

        let point = Vector3D::new(4.2, 1.8, 2.2);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_x_nopbc() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 2.0, 4.0, Dimension::X);

        let point = Vector3D::new(2.9, 3.8, 2.2);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!cylinder.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_x_nopbc2() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 2.0, 4.0, Dimension::X);

        let point = Vector3D::new(3.1, 4.6, 1.2);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_x_pbc_1() {
        let cylinder = Cylinder::new([2.0, 1.0, 3.0].into(), 2.0, 3.0, Dimension::X);

        let point = Vector3D::new(0.3, 1.4, 2.2);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_x_pbc_2() {
        let cylinder = Cylinder::new([2.0, 1.0, 3.0].into(), 2.0, 3.0, Dimension::X);

        let point = Vector3D::new(2.4, 3.8, 2.8);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_y_nopbc() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 4.0, 4.0, Dimension::Y);

        let point = Vector3D::new(5.2, 3.8, 3.2);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_y_nopbc() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 2.0, 4.0, Dimension::Y);

        let point = Vector3D::new(4.2, 7.3, 2.2);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!cylinder.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_y_nopbc2() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 2.0, 4.0, Dimension::Y);

        let point = Vector3D::new(1.1, 5.4, 3.7);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_y_pbc_1() {
        let cylinder = Cylinder::new([1.0, 2.0, 3.0].into(), 2.0, 3.0, Dimension::Y);

        let point = Vector3D::new(0.7, 0.8, 3.4);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_y_pbc_2() {
        let cylinder = Cylinder::new([1.0, 2.0, 3.5].into(), 2.0, 3.0, Dimension::Y);

        let point = Vector3D::new(3.5, 3.4, 0.0);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_z_nopbc() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 2.0, 2.0, Dimension::Z);

        let point = Vector3D::new(4.0, 3.8, 4.8);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_z_nopbc() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 2.0, 4.0, Dimension::Z);

        let point = Vector3D::new(4.0, 3.8, 7.2);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!cylinder.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_z_nopbc2() {
        let cylinder = Cylinder::new([3.0, 3.0, 3.0].into(), 2.0, 4.0, Dimension::Z);

        let point = Vector3D::new(4.9, 3.7, 6.8);
        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_z_pbc_1() {
        let cylinder = Cylinder::new([1.0, 2.0, 3.0].into(), 2.0, 3.0, Dimension::Z);

        let point = Vector3D::new(1.4, 1.5, 1.5);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(cylinder.inside(&point, &simbox));
    }

    #[test]
    fn inside_z_pbc_2() {
        let cylinder = Cylinder::new([3.0, 2.0, 3.0].into(), 2.0, 3.0, Dimension::Z);

        let point = Vector3D::new(0.3, 1.5, 3.1);
        let simbox = SimBox::from([4.0, 4.0, 4.0]);

        assert!(cylinder.inside(&point, &simbox));
    }
}

#[cfg(test)]
mod tests_triprism {
    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn new_z() {
        let base1 = Vector3D::new(3.0, 4.0, 2.0);
        let base2 = Vector3D::new(7.0, 5.0, 2.0);
        let base3 = Vector3D::new(4.0, 3.0, 2.0);

        let prism = TriangularPrism::new(base1, base2, base3, 4.3);

        assert_approx_eq!(f32, prism.get_base1().x, 3.0);
        assert_approx_eq!(f32, prism.get_base1().y, 4.0);
        assert_approx_eq!(f32, prism.get_base1().z, 2.0);

        assert_approx_eq!(f32, prism.get_base2().x, 7.0);
        assert_approx_eq!(f32, prism.get_base2().y, 5.0);
        assert_approx_eq!(f32, prism.get_base2().z, 2.0);

        assert_approx_eq!(f32, prism.get_base3().x, 4.0);
        assert_approx_eq!(f32, prism.get_base3().y, 3.0);
        assert_approx_eq!(f32, prism.get_base3().z, 2.0);

        assert_approx_eq!(f32, prism.get_height(), 4.3);
        assert_eq!(prism.get_orientation(), Dimension::Z);
        assert_eq!(prism.get_plane(), Dimension::XY);
    }

    #[test]
    fn new_y() {
        let base1 = Vector3D::new(3.0, 3.0, 3.0);
        let base2 = Vector3D::new(7.0, 3.0, 2.0);
        let base3 = Vector3D::new(4.0, 3.0, 5.0);

        let prism = TriangularPrism::new(base1, base2, base3, 4.3);

        assert_eq!(prism.get_orientation(), Dimension::Y);
        assert_eq!(prism.get_plane(), Dimension::XZ);
    }

    #[test]
    fn new_x() {
        let base1 = Vector3D::new(5.0, 7.0, 3.0);
        let base2 = Vector3D::new(5.0, 0.0, 2.0);
        let base3 = Vector3D::new(5.0, 4.0, 5.0);

        let prism = TriangularPrism::new(base1, base2, base3, 4.3);

        assert_eq!(prism.get_orientation(), Dimension::X);
        assert_eq!(prism.get_plane(), Dimension::YZ);
    }

    #[test]
    #[should_panic(
        expected = "FATAL GROAN ERROR | TriangularPrism::new | Base of the requested TriangularPrism does not lie in xy, xz, nor yz plane."
    )]
    fn invalid_base_orientation() {
        let base1 = Vector3D::new(3.0, 4.0, 2.0);
        let base2 = Vector3D::new(7.0, 5.0, 1.8);
        let base3 = Vector3D::new(4.0, 3.0, 2.0);

        let _ = TriangularPrism::new(base1, base2, base3, 4.3);
    }

    #[test]
    #[should_panic(
        expected = "FATAL GROAN ERROR | TriangularPrism::new | Base of the requested TriangularPrism can not be constructed."
    )]
    fn invalid_base() {
        let base1 = Vector3D::new(3.0, 4.0, 2.0);
        let base2 = Vector3D::new(7.0, 4.0, 2.0);
        let base3 = Vector3D::new(4.0, 4.0, 2.0);

        let _ = TriangularPrism::new(base1, base2, base3, 4.3);
    }

    #[test]
    fn inside_x_nopbc() {
        let base1 = Vector3D::new(5.0, 7.0, 3.0);
        let base2 = Vector3D::new(5.0, 0.0, 2.0);
        let base3 = Vector3D::new(5.0, 4.0, 5.0);

        let prism = TriangularPrism::new(base1, base2, base3, 4.3);

        let point = Vector3D::new(9.1, 4.8, 3.6);

        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(prism.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_x_nopbc() {
        let base1 = Vector3D::new(5.0, 7.0, 3.0);
        let base2 = Vector3D::new(5.0, 0.0, 2.0);
        let base3 = Vector3D::new(5.0, 4.0, 5.0);

        let prism = TriangularPrism::new(base1, base2, base3, 4.3);

        let point = Vector3D::new(9.7, 4.8, 3.6);

        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!prism.inside(&point, &simbox));
    }

    #[test]
    fn inside_x_pbc() {
        let base1 = Vector3D::new(5.0, 7.0, 3.0);
        let base2 = Vector3D::new(5.0, 0.0, 2.0);
        let base3 = Vector3D::new(5.0, 4.0, 5.0);

        let prism = TriangularPrism::new(base1, base2, base3, 4.3);

        let point = Vector3D::new(0.3, 4.8, 3.6);

        let simbox = SimBox::from([8.0, 8.0, 8.0]);

        assert!(prism.inside(&point, &simbox));
    }

    #[test]
    fn inside_y_nopbc() {
        let base1 = Vector3D::new(3.0, 3.0, 3.0);
        let base2 = Vector3D::new(7.0, 3.0, 2.0);
        let base3 = Vector3D::new(4.0, 3.0, 5.0);

        let prism = TriangularPrism::new(base1, base2, base3, 4.3);

        let point = Vector3D::new(4.8, 5.6, 3.6);

        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(prism.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_y_nopbc() {
        let base1 = Vector3D::new(3.0, 3.0, 3.0);
        let base2 = Vector3D::new(7.0, 3.0, 2.0);
        let base3 = Vector3D::new(4.0, 3.0, 5.0);

        let prism = TriangularPrism::new(base1, base2, base3, 4.3);

        let point = Vector3D::new(5.5, 5.6, 3.6);

        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!prism.inside(&point, &simbox));
    }

    #[test]
    fn inside_y_pbc() {
        let base1 = Vector3D::new(3.0, 3.0, 3.0);
        let base2 = Vector3D::new(7.0, 3.0, 2.0);
        let base3 = Vector3D::new(4.0, 3.0, 5.0);

        let prism = TriangularPrism::new(base1, base2, base3, 4.3);

        let point = Vector3D::new(4.8, 2.1, 3.6);

        let simbox = SimBox::from([10.0, 5.0, 10.0]);

        assert!(prism.inside(&point, &simbox));
    }

    #[test]
    fn inside_z_nopbc() {
        let base1 = Vector3D::new(3.0, 4.0, 2.0);
        let base2 = Vector3D::new(7.0, 5.0, 2.0);
        let base3 = Vector3D::new(4.0, 3.0, 2.0);

        let prism = TriangularPrism::new(base1, base2, base3, 4.3);

        let point = Vector3D::new(4.8, 3.6, 2.1);

        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(prism.inside(&point, &simbox));
    }

    #[test]
    fn not_inside_z_nopbc() {
        let base1 = Vector3D::new(3.0, 4.0, 2.0);
        let base2 = Vector3D::new(7.0, 5.0, 2.0);
        let base3 = Vector3D::new(4.0, 3.0, 2.0);

        let prism = TriangularPrism::new(base1, base2, base3, 4.3);

        let point = Vector3D::new(4.8, 3.4, 2.1);

        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(!prism.inside(&point, &simbox));
    }

    #[test]
    fn inside_z_pbc() {
        let base1 = Vector3D::new(3.0, 4.0, 8.0);
        let base2 = Vector3D::new(7.0, 5.0, 8.0);
        let base3 = Vector3D::new(4.0, 3.0, 8.0);

        let prism = TriangularPrism::new(base1, base2, base3, 4.3);

        let point = Vector3D::new(4.8, 3.6, 2.1);

        let simbox = SimBox::from([10.0, 10.0, 10.0]);

        assert!(prism.inside(&point, &simbox));
    }
}
