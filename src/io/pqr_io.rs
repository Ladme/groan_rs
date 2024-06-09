// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of functions for reading and writing pqr files.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::structures::atom::Atom;
use crate::structures::simbox::SimBox;
use crate::structures::vector3d::Vector3D;
use crate::system::System;

use crate::errors::{ParsePdbError, ParsePqrError, WritePqrError};
use crate::io::pdb_io;

/// Read a pqr file and construct a System structure.
///
/// ## Supported keywords
/// This function can handle lines starting with ATOM, HETATM, TITLE, ENDMDL, END, and CRYST1.
/// All other lines are ignored.
///
/// ## Notes
/// - TITLE and CRYST1 are expected to have the same format as in PDB files and are optional.
///
/// - The fields of each ATOM and HETATM line must be separated by whitespace. The expected format is this:
///
/// `record_name atom_number atom_name residue_name (chain_id) residue_number x y z charge vdw_radius`
///
/// `chain_id` is an optional property and can be skipped.
///
/// - In case multiple TITLE lines are provided, the **last one** is used as the
/// name of the system. If no TITLE line is provided, "Unknown" is used as the name.
///
/// - In case multiple CRYST1 lines are provided, information from the **last one** is used.
/// If no CRYST1 line is provided, the simulation box is undefined.
///
/// - Reading ends once `ENDMDL`, `END`, or the end of file is reached.
/// - The implementation of this function is based on the MDAnalysis' implementation,
/// see `https://docs.mdanalysis.org/2.0.0/documentation_pages/coordinates/PQR.html`
/// but unlike MDAnalysis, this function also reads CRYST1 and TITLE keywords.
pub fn read_pqr(filename: impl AsRef<Path>) -> Result<System, ParsePqrError> {
    let file = match File::open(filename.as_ref()) {
        Ok(x) => x,
        Err(_) => return Err(ParsePqrError::FileNotFound(Box::from(filename.as_ref()))),
    };

    let reader = BufReader::new(file);

    let mut atoms: Vec<Atom> = Vec::new();
    let mut title = "Unknown".to_string();
    let mut simbox = None;

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
        // parse TITLE line
        else if line.len() >= 5 && line[0..5] == *"TITLE" {
            title = pdb_io::line_as_title(&line);
        }
        // parse CRYST1 line
        else if line.len() >= 6 && line[0..6] == *"CRYST1" {
            simbox = match pdb_io::line_as_box(&line) {
                Ok(x) => Some(x),
                Err(ParsePdbError::ParseBoxLineErr(x)) => return Err(ParsePqrError::ParseBoxLineErr(x)),
                Err(e) => panic!("FATAL GROAN ERROR | pqr_io::read_pqr | Unexpected error type {} returned from pdb_io::line_as_box", e),
            }
        }
        // END or ENDMDL is reached => stop reading
        else if line.len() >= 3 && line[0..3] == *"END" {
            break;
        }
    }

    Ok(System::new(&title, atoms, simbox))
}

/// Contains specification of the precision of the PQR file.
#[derive(Debug, Clone, PartialEq)]
pub struct PqrPrecision {
    pub position: usize,
    pub charge: usize,
    pub vdw: usize,
}

impl PqrPrecision {
    /// Create new `PqrPrecision` structure with specified values of precision.
    pub fn new(position: usize, charge: usize, vdw: usize) -> Self {
        PqrPrecision {
            position,
            charge,
            vdw,
        }
    }
}

/// ## Methods for writing pqr files.
impl System {
    /// Write all atoms of the system into a pqr file with the given name.
    ///
    /// ## Parameters
    /// - `precision` parameter specifies the number of decimal places to be printed for
    /// position, charge and radius.
    /// - If not provided, the default values are used. In such case, position coordinates
    /// are written with 3 decimal places, and charge and radius are both written with 4 decimal places.
    ///
    /// ## Returns
    /// `Ok` if writing has been successful. Otherwise `WritePqrError`.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// use groan_rs::io::pqr_io;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// // we use the default PQR file precision
    /// if let Err(e) = system.write_pqr("system.pqr", None) {
    ///     eprintln!("{}", e);
    ///     return;
    /// }
    /// ```
    ///
    /// ## Notes
    /// - Unlike many other programs, `groan_rs` library also uses `CRYST1` and `TITLE` in pqr files.
    /// The pqr file will thus contain information about box dimensions and the name of the system.
    pub fn write_pqr(
        &self,
        filename: impl AsRef<Path>,
        precision: Option<PqrPrecision>,
    ) -> Result<(), WritePqrError> {
        match self.group_write_pqr("all", filename, precision) {
            Ok(_) => Ok(()),
            Err(WritePqrError::GroupNotFound(_)) => {
                panic!(
                    "FATAL GROAN ERROR | System::write_pqr | Default group 'all' does not exist."
                )
            }
            Err(e) => Err(e),
        }
    }

    /// Write atoms of the specified group into a pqr file with the given name.
    ///
    /// ## Parameters
    /// - `precision` parameter specifies the number of decimal places to be printed for
    /// position, charge and radius.
    /// - If not provided, the default values are used. In such case, position coordinates
    /// are written with 3 decimal places, and charge and radius are both written with 4 decimal places.
    ///
    /// ## Returns
    /// `Ok` if writing has been successful. Otherwise `WritePqrError`.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// use groan_rs::io::pqr_io;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// // set precision of the pqr file (optional step)
    /// let precision = pqr_io::PqrPrecision::new(5, 3, 3);
    ///
    /// if let Err(e) = system.group_write_pqr("Protein", "protein.pqr", Some(precision)) {
    ///     eprintln!("{}", e);
    ///     return;
    /// }
    /// ```
    /// ## Notes
    /// - Unlike many other programs, `groan_rs` library also uses `CRYST1` and `TITLE` in pqr files.
    /// The pqr file will thus contain information about box dimensions and the name of the system.
    pub fn group_write_pqr(
        &self,
        group_name: &str,
        filename: impl AsRef<Path>,
        precision: Option<PqrPrecision>,
    ) -> Result<(), WritePqrError> {
        if !self.group_exists(group_name) {
            return Err(WritePqrError::GroupNotFound(group_name.to_string()));
        }

        let output = File::create(&filename)
            .map_err(|_| WritePqrError::CouldNotCreate(Box::from(filename.as_ref())))?;

        let mut writer = BufWriter::new(output);

        // create PqrPrecision structure if not provided
        let precision = precision.unwrap_or(PqrPrecision::new(3, 4, 4));

        let title = match group_name {
            "all" => self.get_name().to_owned(),
            _ => format!("Group `{}` from {}", group_name, self.get_name()),
        };

        write_header(&mut writer, &title, self.get_box_as_ref())?;

        for atom in self.group_iter(group_name).expect(
            "FATAL GROAN ERROR | System::group_write_pqr | Group should exist but it does not.",
        ) {
            atom.write_pqr(&mut writer, &precision)?;
        }

        writeln!(writer, "TER\nEND").map_err(|_| WritePqrError::CouldNotWrite)?;

        Ok(())
    }
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

/// Write a header for a PQR file.
/// If `simbox` is `None`, CRYST line is not written.
fn write_header(
    writer: &mut BufWriter<File>,
    title: &str,
    simbox: Option<&SimBox>,
) -> Result<(), WritePqrError> {
    writeln!(writer, "TITLE     {}", title).map_err(|_| WritePqrError::CouldNotWrite)?;

    if let Some(simbox) = simbox {
        let (lengths, angles) = simbox.to_lengths_angles();
        writeln!(
            writer,
            "CRYST1{:>9.3}{:>9.3}{:>9.3}{:>7.2}{:>7.2}{:>7.2} P 1           1",
            lengths.x * 10.0,
            lengths.y * 10.0,
            lengths.z * 10.0,
            angles.x,
            angles.y,
            angles.z,
        )
        .map_err(|_| WritePqrError::CouldNotWrite)?;
    }

    Ok(())
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

        assert_eq!(system.get_name(), "Buforin II peptide P11L");
        assert_eq!(system.get_n_atoms(), 50);

        // check box size
        let simbox = system.get_box_as_ref().unwrap();
        assert_approx_eq!(f32, simbox.x, 6.0861);
        assert_approx_eq!(f32, simbox.y, 6.0861);
        assert_approx_eq!(f32, simbox.z, 6.0861);

        assert_eq!(simbox.v1y, 0.0f32);
        assert_eq!(simbox.v1z, 0.0f32);
        assert_eq!(simbox.v2x, 0.0f32);

        assert_eq!(simbox.v2z, 0.0f32);
        assert_eq!(simbox.v3x, 0.0f32);
        assert_eq!(simbox.v3y, 0.0f32);

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

        assert_eq!(system.atoms_iter().next().unwrap().get_atom_number(), 1);
        assert_eq!(system.atoms_iter().nth(16).unwrap().get_atom_number(), 17);
    }

    #[test]
    fn read_end() {
        let system = read_pqr("test_files/example_end.pqr").unwrap();

        assert_eq!(system.get_name(), "Unknown");
        assert!(system.get_box_as_ref().is_none());
        assert_eq!(system.get_n_atoms(), 17);

        assert_eq!(system.atoms_iter().next().unwrap().get_atom_number(), 1);
        assert_eq!(system.atoms_iter().nth(16).unwrap().get_atom_number(), 17);
    }

    #[test]
    fn read_nochain() {
        let system_chain = read_pqr("test_files/example.pqr").unwrap();
        let system_nochain = read_pqr("test_files/example_nochain.pqr").unwrap();

        assert_eq!(system_chain.get_name(), system_nochain.get_name());
        assert_eq!(
            system_chain.get_box_as_ref().unwrap().x,
            system_nochain.get_box_as_ref().unwrap().x
        );
        assert_eq!(
            system_chain.get_box_as_ref().unwrap().y,
            system_nochain.get_box_as_ref().unwrap().y
        );
        assert_eq!(
            system_chain.get_box_as_ref().unwrap().z,
            system_nochain.get_box_as_ref().unwrap().z
        );

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
        assert_eq!(
            system1.get_box_as_ref().unwrap().x,
            system2.get_box_as_ref().unwrap().x
        );
        assert_eq!(
            system1.get_box_as_ref().unwrap().y,
            system2.get_box_as_ref().unwrap().y
        );
        assert_eq!(
            system1.get_box_as_ref().unwrap().z,
            system2.get_box_as_ref().unwrap().z
        );

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
        assert_eq!(
            system1.get_box_as_ref().unwrap().x,
            system2.get_box_as_ref().unwrap().x
        );
        assert_eq!(
            system1.get_box_as_ref().unwrap().y,
            system2.get_box_as_ref().unwrap().y
        );
        assert_eq!(
            system1.get_box_as_ref().unwrap().z,
            system2.get_box_as_ref().unwrap().z
        );

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

    read_pqr_fails!(
        read_invalid_box,
        "test_files/example_invalid_box.pqr",
        ParsePqrError::ParseBoxLineErr,
        "CRYST1   60.861   60.861.3 60.861  90.00  90.00  90.00 P 1           1"
    );
}

#[cfg(test)]
mod tests_write {
    use super::*;
    use tempfile::NamedTempFile;

    #[test]
    fn write() {
        let system = System::from_file("test_files/example.pqr").unwrap();

        let pqr_output = NamedTempFile::new().unwrap();
        let path_to_output = pqr_output.path();

        if system.write_pqr(path_to_output, None).is_err() {
            panic!("Writing pqr file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example.pqr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_nochain() {
        let mut system = System::from_file("test_files/example.pqr").unwrap();

        for atom in system.atoms_iter_mut() {
            atom.reset_chain();
        }

        let pqr_output = NamedTempFile::new().unwrap();
        let path_to_output = pqr_output.path();

        if system.write_pqr(path_to_output, None).is_err() {
            panic!("Writing pqr file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example_nochain.pqr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_nobox() {
        let system = System::from_file("test_files/example_nobox.pqr").unwrap();

        let pqr_output = NamedTempFile::new().unwrap();
        let path_to_output = pqr_output.path();

        if system.write_pqr(path_to_output, None).is_err() {
            panic!("Writing pqr file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example_nobox.pqr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group() {
        let mut system = System::from_file("test_files/example.pqr").unwrap();

        system
            .group_create("Selected", "resname ARG PHE LEU")
            .unwrap();

        let pqr_output = NamedTempFile::new().unwrap();
        let path_to_output = pqr_output.path();

        if system
            .group_write_pqr("Selected", path_to_output, None)
            .is_err()
        {
            panic!("Writing pqr file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example_group.pqr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_precision() {
        let system = System::from_file("test_files/example.pqr").unwrap();

        let pqr_output = NamedTempFile::new().unwrap();
        let path_to_output = pqr_output.path();

        let precision = PqrPrecision::new(6, 0, 2);

        if system.write_pqr(path_to_output, Some(precision)).is_err() {
            panic!("Writing pqr file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example_precision.pqr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_large() {
        let mut system = System::from_file("test_files/example.pqr").unwrap();

        system
            .get_atom_as_ref_mut(3)
            .unwrap()
            .set_atom_number(12753);
        system
            .get_atom_as_ref_mut(28)
            .unwrap()
            .set_atom_number(127533497463);
        system
            .get_atom_as_ref_mut(29)
            .unwrap()
            .set_atom_number(999999);
        system
            .get_atom_as_ref_mut(31)
            .unwrap()
            .set_atom_name("SC1234");
        system
            .get_atom_as_ref_mut(2)
            .unwrap()
            .set_residue_name("ARGG");
        system.get_atom_as_ref_mut(17).unwrap().reset_chain();
        system
            .get_atom_as_ref_mut(17)
            .unwrap()
            .set_residue_number(29345);
        system
            .get_atom_as_ref_mut(13)
            .unwrap()
            .set_position_x(14.32);
        system
            .get_atom_as_ref_mut(12)
            .unwrap()
            .set_position_x(214.32134);
        system
            .get_atom_as_ref_mut(12)
            .unwrap()
            .set_position_y(16.21);
        system
            .get_atom_as_ref_mut(11)
            .unwrap()
            .set_position_z(9423.32);
        system.get_atom_as_ref_mut(42).unwrap().set_charge(11.32);
        system.get_atom_as_ref_mut(43).unwrap().set_charge(-11.32);
        system.get_atom_as_ref_mut(45).unwrap().set_vdw(1.477);

        let pqr_output = NamedTempFile::new().unwrap();
        let path_to_output = pqr_output.path();

        if system.write_pqr(path_to_output, None).is_err() {
            panic!("Writing pqr file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example_large.pqr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn read_write_large() {
        let system = System::from_file("test_files/example_large.pqr").unwrap();

        let pqr_output = NamedTempFile::new().unwrap();
        let path_to_output = pqr_output.path();

        if system.write_pqr(path_to_output, None).is_err() {
            panic!("Writing pqr file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example_large.pqr").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
}
