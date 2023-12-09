// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Small functions for testing purposes.

#[cfg(test)]
pub(crate) mod utilities {
    use crate::structures::atom::Atom;
    use float_cmp::assert_approx_eq;

    pub(crate) fn compare_atoms(atom1: &Atom, atom2: &Atom) {
        assert_eq!(atom1.get_residue_number(), atom2.get_residue_number());
        assert_eq!(atom1.get_residue_name(), atom2.get_residue_name());
        assert_eq!(atom1.get_atom_number(), atom2.get_atom_number());
        assert_eq!(atom1.get_atom_name(), atom2.get_atom_name());
        assert_eq!(atom1.get_chain(), atom2.get_chain());
        assert_eq!(atom1.get_mass(), atom2.get_mass());
        assert_eq!(atom1.get_element_name(), atom2.get_element_name());
        assert_eq!(atom1.get_element_symbol(), atom2.get_element_symbol());
        assert_eq!(atom1.get_vdw(), atom2.get_vdw());
        assert_eq!(
            atom1.get_expected_max_bonds(),
            atom2.get_expected_max_bonds()
        );
        assert_eq!(atom1.get_bonded(), atom2.get_bonded());

        if let (Some(pos1), Some(pos2)) = (atom1.get_position(), atom2.get_position()) {
            assert_approx_eq!(f32, pos1.x, pos2.x);
            assert_approx_eq!(f32, pos1.y, pos2.y);
            assert_approx_eq!(f32, pos1.z, pos2.z);
        } else {
            assert!(
                atom1.get_position().is_none() && atom2.get_position().is_none(),
                "Positions are not both None"
            );
        }

        if let (Some(vel1), Some(vel2)) = (atom1.get_velocity(), atom2.get_velocity()) {
            assert_approx_eq!(f32, vel1.x, vel2.x);
            assert_approx_eq!(f32, vel1.y, vel2.y);
            assert_approx_eq!(f32, vel1.z, vel2.z);
        } else {
            assert!(
                atom1.get_velocity().is_none() && atom2.get_velocity().is_none(),
                "Velocities are not both None"
            );
        }

        if let (Some(force1), Some(force2)) = (atom1.get_force(), atom2.get_force()) {
            assert_approx_eq!(f32, force1.x, force2.x);
            assert_approx_eq!(f32, force1.y, force2.y);
            assert_approx_eq!(f32, force1.z, force2.z);
        } else {
            assert!(
                atom1.get_force().is_none() && atom2.get_force().is_none(),
                "Forces are not both None"
            );
        }
    }
}
