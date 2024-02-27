// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of functions for reading tpr files.

use std::path::Path;

use crate::{
    errors::{AtomError, ParseTprError},
    structures::{atom::Atom, simbox::SimBox},
    system::general::System,
};

/// Read topology from a tpr file. Construct a system.
/// Does not read intermolecular bonds, does not read positions/velocities/forces,
/// does not create groups (other than the default `all` and `All`).
/// Note that the atoms and residues in the tpr file are numbered sequentially,
/// no matter their original numbering in the gro/pdb file.
pub(crate) fn read_tpr(filename: impl AsRef<Path>) -> Result<System, ParseTprError> {
    let tpr = match minitpr::TprFile::parse(filename) {
        Ok(x) => x,
        Err(e) => return Err(ParseTprError::CouldNotRead(format!("{}", e))),
    };

    let atoms: Vec<Atom> = tpr
        .topology
        .atoms
        .into_iter()
        .map(|x| -> Atom { Atom::from_minitpr_atom(x) })
        .collect();

    let simbox = tpr.simbox.map(SimBox::from_minitpr_simbox);

    let mut system = System::new(&tpr.system_name, atoms, simbox);

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

impl Atom {
    /// Convert Atom from `minitpr` to native groan_rs Atom.
    fn from_minitpr_atom(tpr_atom: minitpr::Atom) -> Atom {
        let atom = Atom::new(
            tpr_atom.residue_number as usize,
            &tpr_atom.residue_name,
            tpr_atom.atom_number as usize,
            &tpr_atom.atom_name,
        )
        .with_mass(tpr_atom.mass as f32)
        .with_charge(tpr_atom.charge as f32);

        match tpr_atom.element {
            None => atom,
            Some(x) => atom
                .with_element_name(&x.name().to_string().to_lowercase())
                .with_element_symbol(x.symbol()),
        }
    }
}

impl SimBox {
    /// Convert SimBox from `minitpr` to native groan_rs Simbox.
    fn from_minitpr_simbox(tpr_simbox: minitpr::SimBox) -> SimBox {
        SimBox::from([
            tpr_simbox.simbox[0][0] as f32,
            tpr_simbox.simbox[1][1] as f32,
            tpr_simbox.simbox[2][2] as f32,
            tpr_simbox.simbox[0][1] as f32,
            tpr_simbox.simbox[0][2] as f32,
            tpr_simbox.simbox[1][0] as f32,
            tpr_simbox.simbox[1][2] as f32,
            tpr_simbox.simbox[2][0] as f32,
            tpr_simbox.simbox[2][1] as f32,
        ])
    }
}

#[cfg(test)]
mod tests {
    use float_cmp::assert_approx_eq;

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
        };

        let atom2 = minitpr::Atom {
            atom_name: String::from("SC1"),
            atom_number: 2,
            residue_name: String::from("SER"),
            residue_number: 1,
            mass: 36.0,
            charge: 0.0,
            element: None,
        };

        let converted1 = Atom::from_minitpr_atom(atom1);
        let converted2 = Atom::from_minitpr_atom(atom2);

        assert_eq!(converted1.get_atom_name(), "BB");
        assert_eq!(converted1.get_atom_number(), 1);
        assert_eq!(converted1.get_residue_name(), "SER");
        assert_eq!(converted1.get_residue_number(), 1);
        assert_approx_eq!(f32, converted1.get_mass().unwrap(), 72.0);
        assert_approx_eq!(f32, converted1.get_charge().unwrap(), 1.0);
        assert_eq!(converted1.get_element_name(), Some("carbon"));
        assert_eq!(converted1.get_element_symbol(), Some("C"));
        assert!(!converted1.has_position());
        assert!(!converted1.has_velocity());
        assert!(!converted1.has_force());

        assert_eq!(converted2.get_atom_name(), "SC1");
        assert_eq!(converted2.get_atom_number(), 2);
        assert_eq!(converted2.get_residue_name(), "SER");
        assert_eq!(converted2.get_residue_number(), 1);
        assert_approx_eq!(f32, converted2.get_mass().unwrap(), 36.0);
        assert_approx_eq!(f32, converted2.get_charge().unwrap(), 0.0);
        assert_eq!(converted2.get_element_name(), None);
        assert_eq!(converted2.get_element_symbol(), None);
        assert!(!converted2.has_position());
        assert!(!converted2.has_velocity());
        assert!(!converted2.has_force());
    }

    #[test]
    fn convert_box() {
        let tpr_simbox = minitpr::SimBox {
            simbox: [[13.452, 0.0, 0.0], [0.0, 13.452, 0.0], [0.0, 0.0, 15.321]],
            simbox_rel: [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
            simbox_v: [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        };

        let converted = SimBox::from_minitpr_simbox(tpr_simbox);

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

        assert_eq!(system_tpr.get_n_atoms(), 16844);

        let atom1 = system_tpr.get_atom_as_ref(0).unwrap();
        let atom2 = system_tpr.get_atom_as_ref(16843).unwrap();

        assert_eq!(atom1.get_atom_name(), "BB");
        assert_eq!(atom1.get_atom_number(), 1);
        assert_eq!(atom1.get_residue_name(), "GLY");
        assert_eq!(atom1.get_residue_number(), 1);
        assert_approx_eq!(f32, atom1.get_charge().unwrap(), 1.0);
        assert_approx_eq!(f32, atom1.get_mass().unwrap(), 72.0);
        assert!(atom1.get_element_name().is_none());
        assert!(atom1.get_element_symbol().is_none());
        assert!(!atom1.has_position());
        assert!(!atom1.has_velocity());
        assert!(!atom1.has_force());

        assert_eq!(atom2.get_atom_name(), "CL");
        assert_eq!(atom2.get_atom_number(), 16844);
        assert_eq!(atom2.get_residue_name(), "ION");
        assert_eq!(atom2.get_residue_number(), 11180);
        assert_approx_eq!(f32, atom2.get_charge().unwrap(), -1.0);
        assert_approx_eq!(f32, atom2.get_mass().unwrap(), 35.453);
        assert!(atom2.get_element_name().is_none());
        assert!(atom2.get_element_symbol().is_none());
        assert!(!atom2.has_position());
        assert!(!atom2.has_velocity());
        assert!(!atom2.has_force());

        let bonded1 = atom1.get_bonded();
        assert_eq!(bonded1.get_n_atoms(), 1);
        assert!(bonded1.isin(1));

        let bonded2 = atom2.get_bonded();
        assert!(bonded2.is_empty());

        let bonded_special = system_tpr.get_atom_as_ref(2).unwrap().get_bonded();
        assert_eq!(bonded_special.get_n_atoms(), 2);
        assert!(bonded_special.isin(1));
        assert!(bonded_special.isin(3));
    }
}
