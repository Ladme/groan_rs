// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of functions for reading and writing pqr files.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::structures::atom::Atom;
use crate::structures::vector3d::Vector3D;
use crate::system::general::System;

use crate::errors::ParsePqrError;

/// Read a pqr file and construct a System structure.
///
/// ## Supported keywords
/// This function can handle lines starting with ATOM, HETATM, ENDMDL, END.
/// All other lines are ignored. Note that PQR file format does not support
/// providing information about the simulation box and the name of the system.
///
/// ## Notes
/// - The fields on each line must be separated by whitespace. The expected format is this:
///
/// `record_name atom_number atom_name residue_name (chain_id) residue_number x y z charge vdw_radius`
///
/// `chain_id` is an optional property and can be skipped.
///
/// - Reading ends once `ENDMDL`, `END`, or the end of file is reached.
/// - Name of the system is set to "Unknown". Simulation box of the system is undefined.
/// - The implementation of this function is based on the MDAnalysis' implementation,
/// see `https://docs.mdanalysis.org/2.0.0/documentation_pages/coordinates/PQR.html`.
pub fn read_pqr(filename: impl AsRef<Path>) -> Result<System, ParsePqrError> {
    let file = match File::open(filename.as_ref()) {
        Ok(x) => x,
        Err(_) => return Err(ParsePqrError::FileNotFound(Box::from(filename.as_ref()))),
    };

    let reader = BufReader::new(file);

    let mut atoms: Vec<Atom> = Vec::new();

    for raw_line in reader.lines() {
        let line = match raw_line {
            Ok(x) => x,
            Err(_) => return Err(ParsePqrError::LineNotFound(Box::from(filename.as_ref()))),
        };

        // parse ATOM/HETATM line
        if (line.len() >= 4 && line[0..4] == *"ATOM")
            || (line.len() >= 6 && line[0..6] == *"HETATM")
        {
            atoms.push(line_as_atom(&line)?);
        }
        // END or ENDMDL is reached => stop reading
        else if line.len() >= 3 && line[0..3] == *"END" {
            break;
        }
    }

    Ok(System::new("Unknown", atoms, None))
}

fn line_as_atom(line: &str) -> Result<Atom, ParsePqrError> {
    // split the line
    let split: Vec<&str> = line.split_whitespace().collect();

    // decide whether chain_id is available
    let convert = match split.len() {
        11 => 0,
        10 => 1,
        _ => return Err(ParsePqrError::ParseAtomLineErr(line.to_string())),
    };

    // parse properties of the atoms
    let atom_number = split[1]
        .parse::<usize>()
        .map_err(|_| ParsePqrError::ParseAtomLineErr(line.to_string()))?;
    let atom_name = split[2];
    let residue_name = split[3];

    let chain = if convert == 0 {
        if split[4].len() == 1 {
            Some(split[4].chars().next().expect(
                "FATAL GROAN ERROR | pqr_io::line_as_atom | Chain ID should be one character.",
            ))
        } else {
            return Err(ParsePqrError::ParseAtomLineErr(line.to_string()));
        }
    } else {
        None
    };

    let residue_number = split[5 - convert]
        .parse::<usize>()
        .map_err(|_| ParsePqrError::ParseAtomLineErr(line.to_string()))?;

    let x = parse_float(line, split[6 - convert])? / 10.0;
    let y = parse_float(line, split[7 - convert])? / 10.0;
    let z = parse_float(line, split[8 - convert])? / 10.0;

    let charge = parse_float(line, split[9 - convert])?;
    let vdw = parse_float(line, split[10 - convert])? / 10.0;

    let atom = Atom::new(residue_number, residue_name, atom_number, atom_name)
        .with_position(Vector3D::new(x, y, z))
        .with_charge(charge)
        .with_vdw(vdw);

    // add chain information, if available
    match chain {
        Some(x) => Ok(atom.with_chain(x)),
        None => Ok(atom),
    }
}

fn parse_float(line: &str, string: &str) -> Result<f32, ParsePqrError> {
    string
        .parse::<f32>()
        .map_err(|_| ParsePqrError::ParseAtomLineErr(line.to_string()))
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests_read {
    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn read_simple() {
        let system = read_pqr("test_files/example.pqr").unwrap();

        assert_eq!(system.get_name(), "Unknown");
        assert_eq!(system.get_n_atoms(), 50);
        assert!(system.get_box_as_ref().is_none());

        let atoms = system.get_atoms_as_ref();

        // check the first atom
        let first = &atoms[0];
        assert_eq!(first.get_residue_number(), 1);
        assert_eq!(first.get_residue_name(), "THR");
        assert_eq!(first.get_atom_name(), "BB");
        assert_eq!(first.get_atom_number(), 1);
        assert_eq!(first.get_chain().unwrap(), 'A');
        assert_eq!(first.get_charge().unwrap(), 1.0);
        assert_eq!(first.get_vdw().unwrap(), 0.28);

        assert_approx_eq!(f32, first.get_position().unwrap().x, 1.660);
        assert_approx_eq!(f32, first.get_position().unwrap().y, 2.061);
        assert_approx_eq!(f32, first.get_position().unwrap().z, 3.153);

        // check atom somewhere in the middle
        let middle = &atoms[24];
        assert_eq!(middle.get_residue_number(), 11);
        assert_eq!(middle.get_residue_name(), "LEU");
        assert_eq!(middle.get_atom_name(), "SC1");
        assert_eq!(middle.get_atom_number(), 25);
        assert_eq!(middle.get_chain().unwrap(), 'B');
        assert_eq!(middle.get_charge().unwrap(), 0.3);
        assert_eq!(middle.get_vdw().unwrap(), 0.35);

        assert_approx_eq!(f32, middle.get_position().unwrap().x, 3.161);
        assert_approx_eq!(f32, middle.get_position().unwrap().y, 2.868);
        assert_approx_eq!(f32, middle.get_position().unwrap().z, 2.797);

        // check the last atom
        let last = &atoms[49];
        assert_eq!(last.get_residue_number(), 21);
        assert_eq!(last.get_residue_name(), "LYS");
        assert_eq!(last.get_atom_name(), "SC2");
        assert_eq!(last.get_atom_number(), 50);
        assert_eq!(last.get_chain().unwrap(), 'C');
        assert_eq!(last.get_charge().unwrap(), -1.3);
        assert_eq!(last.get_vdw().unwrap(), 0.07);

        assert_approx_eq!(f32, last.get_position().unwrap().x, 4.706);
        assert_approx_eq!(f32, last.get_position().unwrap().y, 4.447);
        assert_approx_eq!(f32, last.get_position().unwrap().z, 2.813);

        // check that the velocity and force of all atoms is zero
        for atom in atoms.iter() {
            assert_eq!(atom.get_velocity(), None);
            assert_eq!(atom.get_force(), None);
        }
    }

    #[test]
    fn read_endmdl() {
        let system = read_pqr("test_files/example_endmdl.pqr").unwrap();

        assert_eq!(system.get_name(), "Unknown");
        assert!(system.get_box_as_ref().is_none());
        assert_eq!(system.get_n_atoms(), 17);

        assert_eq!(system.atoms_iter().nth(0).unwrap().get_atom_number(), 1);
        assert_eq!(system.atoms_iter().nth(16).unwrap().get_atom_number(), 17);
    }

    #[test]
    fn read_end() {
        let system = read_pqr("test_files/example_end.pqr").unwrap();

        assert_eq!(system.get_name(), "Unknown");
        assert!(system.get_box_as_ref().is_none());
        assert_eq!(system.get_n_atoms(), 17);

        assert_eq!(system.atoms_iter().nth(0).unwrap().get_atom_number(), 1);
        assert_eq!(system.atoms_iter().nth(16).unwrap().get_atom_number(), 17);
    }

    #[test]
    fn read_nochain() {
        let system_chain = read_pqr("test_files/example.pqr").unwrap();
        let system_nochain = read_pqr("test_files/example_nochain.pqr").unwrap();

        assert_eq!(system_chain.get_name(), system_nochain.get_name());
        assert!(!system_chain.has_box());
        assert!(!system_nochain.has_box());

        for (ac, anc) in system_chain.atoms_iter().zip(system_nochain.atoms_iter()) {
            assert_eq!(ac.get_residue_number(), anc.get_residue_number());
            assert_eq!(ac.get_residue_name(), anc.get_residue_name());
            assert_eq!(ac.get_atom_number(), anc.get_atom_number());
            assert_eq!(ac.get_atom_name(), anc.get_atom_name());
            assert_eq!(ac.get_position().unwrap().x, anc.get_position().unwrap().x);
            assert_eq!(ac.get_position().unwrap().y, anc.get_position().unwrap().y);
            assert_eq!(ac.get_position().unwrap().z, anc.get_position().unwrap().z);
            assert_eq!(ac.get_charge(), anc.get_charge());
            assert_eq!(ac.get_vdw(), anc.get_vdw());

            assert_eq!(ac.get_velocity(), anc.get_velocity());
            assert_eq!(ac.get_force(), anc.get_force());

            assert_eq!(anc.get_chain(), None);
        }
    }

    #[test]
    fn read_weird_format() {
        let system1 = read_pqr("test_files/example.pqr").unwrap();
        let system2 = read_pqr("test_files/example_weird_format.pqr").unwrap();

        assert_eq!(system1.get_name(), system2.get_name());
        assert!(!system1.has_box());
        assert!(!system2.has_box());

        for (a1, a2) in system1.atoms_iter().zip(system2.atoms_iter()) {
            assert_eq!(a1.get_residue_number(), a2.get_residue_number());
            assert_eq!(a1.get_residue_name(), a2.get_residue_name());
            assert_eq!(a1.get_atom_number(), a2.get_atom_number());
            assert_eq!(a1.get_atom_name(), a2.get_atom_name());
            assert_eq!(a1.get_position().unwrap().x, a2.get_position().unwrap().x);
            assert_eq!(a1.get_position().unwrap().y, a2.get_position().unwrap().y);
            assert_eq!(a1.get_position().unwrap().z, a2.get_position().unwrap().z);
            assert_eq!(a1.get_charge(), a2.get_charge());
            assert_eq!(a1.get_vdw(), a2.get_vdw());

            assert_eq!(a1.get_velocity(), a2.get_velocity());
            assert_eq!(a1.get_force(), a2.get_force());

            assert_eq!(a1.get_chain(), a2.get_chain());
        }
    }

    #[test]
    fn read_mixed_chain() {
        let system1 = read_pqr("test_files/example.pqr").unwrap();
        let system2 = read_pqr("test_files/example_mixchain.pqr").unwrap();

        assert_eq!(system1.get_name(), system2.get_name());
        assert!(!system1.has_box());
        assert!(!system2.has_box());

        for (a1, a2) in system1.atoms_iter().zip(system2.atoms_iter()) {
            assert_eq!(a1.get_residue_number(), a2.get_residue_number());
            assert_eq!(a1.get_residue_name(), a2.get_residue_name());
            assert_eq!(a1.get_atom_number(), a2.get_atom_number());
            assert_eq!(a1.get_atom_name(), a2.get_atom_name());
            assert_eq!(a1.get_position().unwrap().x, a2.get_position().unwrap().x);
            assert_eq!(a1.get_position().unwrap().y, a2.get_position().unwrap().y);
            assert_eq!(a1.get_position().unwrap().z, a2.get_position().unwrap().z);
            assert_eq!(a1.get_charge(), a2.get_charge());
            assert_eq!(a1.get_vdw(), a2.get_vdw());

            assert_eq!(a1.get_velocity(), a2.get_velocity());
            assert_eq!(a1.get_force(), a2.get_force());

            if a2.get_atom_number() == 12
                || a2.get_atom_number() == 27
                || a2.get_atom_number() == 43
            {
                assert!(a2.get_chain().is_none());
            } else {
                assert_eq!(a1.get_chain(), a2.get_chain());
            }
        }
    }

    macro_rules! read_pqr_fails {
        ($name:ident, $file:expr, $variant:path, $expected:expr) => {
            #[test]
            fn $name() {
                let file = $file;
                match read_pqr(file) {
                    Err($variant(e)) => assert_eq!(e, $expected),
                    Ok(_) => panic!("Parsing should have failed, but it succeeded."),
                    Err(e) => panic!("Parsing successfully failed but incorrect error type `{:?}` was returned.", e),
                }
            }
        };
    }

    read_pqr_fails!(
        read_invalid_chain,
        "test_files/example_invalid_chain.pqr",
        ParsePqrError::ParseAtomLineErr,
        "ATOM     21  SC1 PHE BX 10      28.140  28.580  33.590  0.0000 0.0000"
    );

    read_pqr_fails!(
        read_invalid_coord,
        "test_files/example_invalid_coord.pqr",
        ParsePqrError::ParseAtomLineErr,
        "ATOM     32  BB  VAL C  15      35.650  33,730  31.780  0.0000 0.0000"
    );

    read_pqr_fails!(
        read_invalid_vdw,
        "test_files/example_invalid_vdw.pqr",
        ParsePqrError::ParseAtomLineErr,
        "ATOM      8  BB  SER A   4      20.040  24.210  30.630  0.0000 0.0O00"
    );
}
