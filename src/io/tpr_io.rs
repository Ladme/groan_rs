// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of functions for reading tpr files.

use std::path::Path;

use crate::{
    errors::{AtomError, ParseTprError},
    structures::{atom::Atom, simbox::SimBox},
    system::System,
};

/// Read topology from a tpr file. Construct a system.
///
/// ## Notes
/// - Does not create groups apart from the default `all` and `All`.
/// - The atoms and residues in the tpr file are numbered sequentially,
/// no matter their original numbering in the gro/pdb file.
pub fn read_tpr(filename: impl AsRef<Path>) -> Result<System, ParseTprError> {
    let tpr = match minitpr::TprFile::parse(filename) {
        Ok(x) => x,
        Err(e) => return Err(ParseTprError::CouldNotRead(format!("{}", e))),
    };

    let atoms: Vec<Atom> = tpr.topology.atoms.into_iter().map(Into::into).collect();

    let simbox = tpr.simbox.map(minitpr::SimBox::into);

    let mut system = System::new(&tpr.system_name, atoms, simbox);
    system.set_lambda(tpr.header.lambda as f32);

    for bond in tpr.topology.bonds {
        match system.add_bond(bond.atom1, bond.atom2) {
            Ok(_) => (),
            Err(error) => match error {
                AtomError::InvalidBond(i1, i2) => return Err(ParseTprError::InvalidBond(i1, i2)),
                AtomError::OutOfRange(x) => panic!("FATAL GROAN ERROR | tpr_io::read_tpr | Bond could not be created for atom index {}. This atom does not exist. `minitpr` [version {}] probably failed to parse the tpr file correctly.", x, minitpr::MINITPR_VERSION),
                _ => panic!("FATAL GROAN ERROR | tpr_io::read_tpr | Bond could not be created for an unknown reason (`system.add_bond` returned an unexpected error type)"),
            }
        }
    }

    Ok(system)
}

impl From<minitpr::Atom> for Atom {
    /// Convert `Atom` from `minitpr` crate to native `groan_rs` Atom.
    fn from(value: minitpr::Atom) -> Self {
        let mut atom = Self::new(
            value.residue_number as usize,
            &value.residue_name,
            value.atom_number as usize,
            &value.atom_name,
        )
        .with_mass(value.mass as f32)
        .with_charge(value.charge as f32);

        if let Some(pos) = value.position {
            atom.set_position(pos.into())
        }

        if let Some(vel) = value.velocity {
            atom.set_velocity(vel.into())
        }

        if let Some(force) = value.force {
            atom.set_force(force.into())
        }

        if let Some(element) = value.element {
            atom.set_element_name(&element.name().to_string().to_lowercase());
            atom.set_element_symbol(element.symbol());
        }

        atom
    }
}

impl From<minitpr::SimBox> for SimBox {
    /// Convert SimBox from `minitpr` crate to native `groan_rs` Simbox.
    fn from(value: minitpr::SimBox) -> Self {
        Self::from([
            value.simbox[0][0] as f32,
            value.simbox[1][1] as f32,
            value.simbox[2][2] as f32,
            value.simbox[0][1] as f32,
            value.simbox[0][2] as f32,
            value.simbox[1][0] as f32,
            value.simbox[1][2] as f32,
            value.simbox[2][0] as f32,
            value.simbox[2][1] as f32,
        ])
    }
}

#[cfg(test)]
mod tests {
    use float_cmp::assert_approx_eq;

    use crate::{io::gro_io::read_gro, structures::vector3d::Vector3D};

    use super::*;

    #[test]
    fn convert_atoms() {
        let atom1 = minitpr::Atom {
            atom_name: String::from("BB"),
            atom_number: 1,
            residue_name: String::from("SER"),
            residue_number: 1,
            mass: 72.0,
            charge: 1.0,
            element: Some(minitpr::Element::C),
            position: Some([1.497, 2.432, 0.439]),
            velocity: None,
            force: Some([17.3424, -7.2974, 0.3412]),
        };

        let atom2 = minitpr::Atom {
            atom_name: String::from("SC1"),
            atom_number: 2,
            residue_name: String::from("SER"),
            residue_number: 1,
            mass: 36.0,
            charge: 0.0,
            element: None,
            position: None,
            velocity: Some([-1.424, 0.439, -2.432]),
            force: None,
        };

        let converted1: Atom = atom1.into();
        let converted2: Atom = atom2.into();

        assert_eq!(converted1.get_atom_name(), "BB");
        assert_eq!(converted1.get_atom_number(), 1);
        assert_eq!(converted1.get_residue_name(), "SER");
        assert_eq!(converted1.get_residue_number(), 1);
        assert_approx_eq!(f32, converted1.get_mass().unwrap(), 72.0);
        assert_approx_eq!(f32, converted1.get_charge().unwrap(), 1.0);
        assert_eq!(converted1.get_element_name(), Some("carbon"));
        assert_eq!(converted1.get_element_symbol(), Some("C"));
        assert_eq!(
            converted1.get_position(),
            Some(&Vector3D::new(1.497, 2.432, 0.439))
        );
        assert!(!converted1.has_velocity());
        assert_eq!(
            converted1.get_force(),
            Some(&Vector3D::new(17.3424, -7.2974, 0.3412))
        );

        assert_eq!(converted2.get_atom_name(), "SC1");
        assert_eq!(converted2.get_atom_number(), 2);
        assert_eq!(converted2.get_residue_name(), "SER");
        assert_eq!(converted2.get_residue_number(), 1);
        assert_approx_eq!(f32, converted2.get_mass().unwrap(), 36.0);
        assert_approx_eq!(f32, converted2.get_charge().unwrap(), 0.0);
        assert_eq!(converted2.get_element_name(), None);
        assert_eq!(converted2.get_element_symbol(), None);
        assert!(!converted2.has_position());
        assert_eq!(
            converted2.get_velocity(),
            Some(&Vector3D::new(-1.424, 0.439, -2.432))
        );
        assert!(!converted2.has_force());
    }

    #[test]
    fn convert_box() {
        let tpr_simbox = minitpr::SimBox {
            simbox: [[13.452, 0.0, 0.0], [0.0, 13.452, 0.0], [0.0, 0.0, 15.321]],
            simbox_rel: [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            simbox_v: [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        };

        let converted: SimBox = tpr_simbox.into();

        assert_approx_eq!(f32, converted.x, 13.452);
        assert_approx_eq!(f32, converted.y, 13.452);
        assert_approx_eq!(f32, converted.z, 15.321);

        assert_approx_eq!(f32, converted.v1y, 0.0);
        assert_approx_eq!(f32, converted.v1z, 0.0);

        assert_approx_eq!(f32, converted.v2x, 0.0);
        assert_approx_eq!(f32, converted.v2z, 0.0);

        assert_approx_eq!(f32, converted.v3x, 0.0);
        assert_approx_eq!(f32, converted.v3y, 0.0);
    }

    #[test]
    fn parse_tpr() {
        let system_tpr = read_tpr("test_files/example.tpr").unwrap();
        let system_gro = read_gro("test_files/example.gro").unwrap();

        assert_eq!(system_tpr.get_n_atoms(), system_gro.get_n_atoms());
        assert_eq!(system_tpr.get_lambda(), 0.0);

        for (a1, a2) in system_tpr.atoms_iter().zip(system_gro.atoms_iter()) {
            crate::test_utilities::utilities::compare_atoms_tpr_with_gro(a1, a2);
            assert!(!a1.has_force());
        }

        // check bonds
        system_tpr.get_atom_as_ref(0).unwrap().get_bonded().isin(1);
        system_tpr.get_atom_as_ref(1).unwrap().get_bonded().isin(0);
        system_tpr.get_atom_as_ref(1).unwrap().get_bonded().isin(2);
        system_tpr.get_atom_as_ref(1).unwrap().get_bonded().isin(4);
        system_tpr
            .get_atom_as_ref(854)
            .unwrap()
            .get_bonded()
            .isin(853);
        system_tpr
            .get_atom_as_ref(854)
            .unwrap()
            .get_bonded()
            .isin(855);
    }
}
