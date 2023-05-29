
use crate::simbox::SimBox;

#[derive(Debug)]
pub struct Vector3D {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl From<[f32; 3]> for Vector3D {
    fn from(arr: [f32; 3]) -> Self {
        Vector3D { x: arr[0], y: arr[1], z: arr[2] }
    }
}

impl Vector3D {

    /// Wraps coordinates of Vector3D so that each of them fits into the simulation box.
    /// Assumes orthogonal simulation box.
    pub fn wrap(&mut self, sbox: &SimBox) {
        self.x = Vector3D::wrap_coordinate(self.x, sbox.x);
        self.y = Vector3D::wrap_coordinate(self.y, sbox.y);
        self.z = Vector3D::wrap_coordinate(self.z, sbox.z);
    }

    /// Wraps a single coordinate into the simulation box.
    fn wrap_coordinate(coor: f32, box_len: f32) -> f32 {
        let wrapped = coor % box_len;
        if wrapped < 0.0 { wrapped + box_len }
        else { wrapped }
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn test_wrap() {

        let mut vector1 = Vector3D::from([-1.0, 1.5, 3.0]);
        let mut vector2 = Vector3D::from([2.0, 2.2, -0.3]);
        let mut vector3 = Vector3D::from([-54.2, 77.8, 124.5]);
        let simbox = SimBox::from([2.0, 2.0, 2.0]);

        vector1.wrap(&simbox);
        assert_approx_eq!(f32, vector1.x, 1.0);
        assert_approx_eq!(f32, vector1.y, 1.5);
        assert_approx_eq!(f32, vector1.z, 1.0);

        vector2.wrap(&simbox);
        assert_approx_eq!(f32, vector2.x, 0.0);
        assert_approx_eq!(f32, vector2.y, 0.2);
        assert_approx_eq!(f32, vector2.z, 1.7);

        vector3.wrap(&simbox);
        assert_approx_eq!(f32, vector3.x, 1.8, epsilon = 0.00001);
        assert_approx_eq!(f32, vector3.y, 1.8, epsilon = 0.00001);
        assert_approx_eq!(f32, vector3.z, 0.5, epsilon = 0.00001);
    }
}