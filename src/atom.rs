// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of the Atom structure and associated methods.

use crate::vector3d::Vector3D;
use crate::simbox::SimBox;

#[derive(Debug)]
pub struct Atom {
    residue_number: u32,
    residue_name: String,
    atom_number: u32,
    gmx_atom_number: u64,
    atom_name: String,
    position: Vector3D,
    velocity: Vector3D,
    force: Vector3D,
}

impl Atom {

    /// Create new Atom structure with the specified properties.
    pub fn new(
            residue_number: u32, residue_name: &str, 
            atom_number: u32, gmx_atom_number: u64, atom_name: &str, 
            position: Vector3D, velocity: Vector3D, force: Vector3D) -> Self {
                Atom { 
                    residue_number, 
                    residue_name: residue_name.to_string(), 
                    atom_number, 
                    gmx_atom_number, 
                    atom_name: atom_name.to_string(), 
                    position, 
                    velocity, 
                    force }
            }
    
    /// Get the number of the residue to which the atom belongs.
    pub fn get_residue_number(&self) -> u32 { self.residue_number }
    
    /// Get the name of the residue to which the atom belongs.
    pub fn get_residue_name(&self) -> &str { &self.residue_name }

    /// Get the number of the atom as presented in gro file.
    pub fn get_atom_number(&self) -> u32 { self.atom_number }

    /// Get the number of the atom as used by gromacs.
    pub fn get_gmx_atom_number(&self) -> u64 { self.gmx_atom_number }

    /// Get the name of the atom.
    pub fn get_atom_name(&self) -> &str { &self.atom_name }

    /// Get the coordinates of the atom.
    pub fn get_position(&self) -> &Vector3D { &self.position }

    /// Get the velocity vector of the atom.
    pub fn get_velocity(&self) -> &Vector3D { &self.velocity }

    /// Get the vector of the total force acting on the atom.
    pub fn get_force(&self) -> &Vector3D { &self.force }

    /// Translates the position of the atom by the provided Vector3D.
    /// Wraps the atom to the simulation box.
    pub fn translate(&mut self, translate: &Vector3D, sbox: &SimBox) {
        self.position.x += translate.x;
        self.position.y += translate.y;
        self.position.z += translate.z;

        self.position.wrap(sbox);
    }

    /// Translates the position of the atom by the provided Vector3D.
    /// DOES NOT wrap the atom to the simulation box.
    pub fn translate_nopbc(&mut self, translate: &Vector3D) {
        self.position.x += translate.x;
        self.position.y += translate.y;
        self.position.z += translate.z;
    }

}



#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn test_new_atom() {

        let atom = Atom::new(45, "GLY", 123, 123, "BB", [15.123, 14.321, 9.834].into(), [-3.432, 0.184, 1.234].into(), [5.1235, 2.3451, -0.32145].into());

        assert_eq!(atom.get_residue_number(), 45);
        assert_eq!(atom.get_residue_name(), "GLY");
        assert_eq!(atom.get_atom_number(), 123);
        assert_eq!(atom.get_gmx_atom_number(), 123);
        assert_eq!(atom.get_atom_name(), "BB");

        assert_approx_eq!(f32, atom.get_position().x, 15.123);
        assert_approx_eq!(f32, atom.get_position().y, 14.321);
        assert_approx_eq!(f32, atom.get_position().z, 9.834);

        assert_approx_eq!(f32, atom.get_velocity().x, -3.432);
        assert_approx_eq!(f32, atom.get_velocity().y, 0.184);
        assert_approx_eq!(f32, atom.get_velocity().z, 1.234);

        assert_approx_eq!(f32, atom.get_force().x, 5.1235);
        assert_approx_eq!(f32, atom.get_force().y, 2.3451);
        assert_approx_eq!(f32, atom.get_force().z, -0.32145);

    }
}