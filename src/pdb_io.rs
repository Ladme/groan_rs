// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of functions for reading and writing pdb files.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::atom::Atom;
use crate::errors::{ParsePdbError, WritePdbError};
use crate::simbox::SimBox;
use crate::system::System;
use crate::vector3d::Vector3D;

/// Read a pdb file and construct a System structure.
///
/// ## Warning
/// Currently only supports orthogonal simulation boxes!
///
/// ## Supported keywords
/// This function can handle lines starting with ATOM, HETATM, TITLE and CRYST1.
/// All other lines are ignored.
///
/// ## Notes
/// - In case multiple TITLE lines are provided, the **last one** is used as the
/// name of the system. If no TITLE line is provided, "Unknown" is used as the name.
///
/// - In case multiple CRYST1 lines are provided, information from the **last one** is used.
/// If no CRYST1 line is provided, the simulation box size is set to 0 in all dimensions.
/// Note again that only orthogonal simulation boxes are supported.
///
pub fn read_pdb(filename: impl AsRef<Path>) -> Result<System, ParsePdbError> {
    let file = match File::open(filename.as_ref()) {
        Ok(x) => x,
        Err(_) => return Err(ParsePdbError::FileNotFound(Box::from(filename.as_ref()))),
    };

    let reader = BufReader::new(file);

    let mut atoms: Vec<Atom> = Vec::new();
    let mut title = "Unknown".to_string();
    let mut simbox = SimBox::from([0.0, 0.0, 0.0]);

    for raw_line in reader.lines() {
        let line = match raw_line {
            Ok(x) => x,
            Err(_) => return Err(ParsePdbError::LineNotFound(Box::from(filename.as_ref()))),
        };

        // parse ATOM/HETATM line
        if (line.len() >= 4 && line[0..4] == "ATOM".to_string())
            || (line.len() >= 6 && line[0..6] == "HETATM".to_string())
        {
            atoms.push(line_as_atom(&line)?);
        }
        // parse TITLE line
        else if line.len() >= 5 && line[0..5] == "TITLE".to_string() {
            title = line_as_title(&line)?;
        }
        // parse CRYST1 line
        else if line.len() >= 6 && line[0..6] == "CRYST1".to_string() {
            simbox = line_as_box(&line)?;
        }
    }

    Ok(System::new(&title, atoms, simbox))
}

impl System {

    /// Write all atoms of the `System` into a pdb file with the given name.
    ///
    /// ## Returns
    /// `Ok` if writing has been successful. Otherwise `WritePdbError`.
    ///
    /// ## Warning
    /// This function only supports orthogonal rectangular simulation boxes.
    /// 
    /// ## Example
    /// ```no_run
    /// use groan_rs::System;
    ///
    /// let system = System::from_file("system.gro").unwrap();
    /// if let Err(e) = system.write_pdb("system.pdb") {
    ///     eprintln!("{}", e);
    ///     return;
    /// }
    /// ```
    /// 
    /// ## Notes
    /// - The chain identifier will not be written for any atom.
    pub fn write_pdb(&self, filename: impl AsRef<Path>) -> Result<(), WritePdbError> {
        match self.group_write_pdb("all", filename) {
            Ok(_) => Ok(()),
            Err(WritePdbError::GroupNotFound(_)) => panic!("Groan error. Default group 'all' does not exist."),
            Err(e) => Err(e),
        }
    }


    /// Write atoms of the specified group into a pdb file with the given name.
    ///
    /// ## Returns
    /// `Ok` if writing has been successful. Otherwise `WritePdbError`.
    ///
    /// ## Warning
    /// This function only supports orthogonal rectangular simulation boxes.
    /// 
    /// ## Example
    /// ```no_run
    /// use groan_rs::System;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// if let Err(e) = system.group_write_pdb("Protein", "protein.pdb") {
    ///     eprintln!("{}", e);
    ///     return;
    /// }
    /// ```
    /// 
    /// ## Notes
    /// - The chain identifier will not be written for any atom.
    pub fn group_write_pdb(&self, group_name: &str, filename: impl AsRef<Path>) -> Result<(), WritePdbError> {
        if !self.group_exists(group_name) {
            return Err(WritePdbError::GroupNotFound(group_name.to_string()));
        }

        let output = File::create(&filename)
            .map_err(|_| WritePdbError::CouldNotCreate(Box::from(filename.as_ref())))?;

        let mut writer = BufWriter::new(output);

        let title = match group_name {
            "all" => self.get_name().to_owned(),
            _ => format!("Group `{}` from {}", group_name, self.get_name()),
        };

        write_header(&mut writer, &title, self.get_box_as_ref())?;

        for atom in self.group_iter(group_name).unwrap() {
            atom.write_pdb(&mut writer)?;
        }

        write_line(&mut writer, "TER\nENDMDL")?;

        writer.flush().map_err(|_| WritePdbError::CouldNotWrite())?;

        Ok(())
    }
}

/// Parse a single line from a pdb file as an atom.
///
/// ## Notes
/// - Parses lines starting with ATOM or HETATM.
fn line_as_atom(line: &str) -> Result<Atom, ParsePdbError> {
    // check line length
    if line.len() < 54 {
        return Err(ParsePdbError::ParseAtomLineErr(line.to_string()));
    }

    // parsing atom number
    let atom_number = line[6..11]
        .trim()
        .parse::<usize>()
        .map_err(|_| ParsePdbError::ParseAtomLineErr(line.to_string()))?;

    // parsing atom name
    let atom_name = line[12..16].trim().to_string();
    if atom_name.is_empty() {
        return Err(ParsePdbError::ParseAtomLineErr(line.to_string()));
    }

    // parsing residue name
    let residue_name = line[17..21].trim().to_string();
    if residue_name.is_empty() {
        return Err(ParsePdbError::ParseAtomLineErr(line.to_string()));
    }

    // parsing residue number
    let residue_number = line[22..26]
        .trim()
        .parse::<usize>()
        .map_err(|_| ParsePdbError::ParseAtomLineErr(line.to_string()))?;

    // parsing position
    let mut curr = 30usize;
    let mut position = [0.0, 0.0, 0.0];
    for i in 0..3 {
        position[i] = line[curr..curr + 8]
            .trim()
            .parse::<f32>()
            .map(|x| x / 10.0)
            .map_err(|_| ParsePdbError::ParseAtomLineErr(line.to_string()))?;

        curr += 8;
    }

    Ok(Atom::new(
        residue_number,
        &residue_name,
        atom_number,
        &atom_name,
        position.into(),
        Vector3D::default(),
        Vector3D::default(),
    ))
}

/// Parse a single line as a simulation box.
///
/// ## Notes
/// - Parses a line starting with CRYST1.
fn line_as_box(line: &str) -> Result<SimBox, ParsePdbError> {
    // check line length
    if line.len() < 54 {
        return Err(ParsePdbError::ParseBoxLineErr(line.to_string()));
    }

    let mut boxsize = [0.0, 0.0, 0.0];
    let mut curr = 6usize;
    for i in 0..3 {
        boxsize[i] = line[curr..curr + 9]
            .trim()
            .parse::<f32>()
            .map(|x| x / 10.0)
            .map_err(|_| ParsePdbError::ParseBoxLineErr(line.to_string()))?;

        curr += 9;
    }

    // check that the box is orthogonal
    for _ in 0..3 {
        let value = line[curr..curr + 7]
            .trim()
            .parse::<f32>()
            .map_err(|_| ParsePdbError::ParseBoxLineErr(line.to_string()))?;

        if value != 90.0 {
            return Err(ParsePdbError::NonOrthogonalBox(line.to_string()));
        }

        curr += 7;
    }

    Ok(SimBox::from(boxsize))
}

/// Parse a single line as a title.
///
/// ## Notes
/// - Parses a line starting with TITLE.
/// - In case the TITLE line is empty, 'Unknown' is used as the name for the system.
fn line_as_title(line: &str) -> Result<String, ParsePdbError> {
    let title = line[5..].trim().to_string();
    if title.is_empty() {
        return Ok("Unknown".to_string());
    }

    Ok(title)
}

fn write_line<W: Write>(writer: &mut W, line: &str) -> Result<(), WritePdbError> {
    writeln!(writer, "{}", line).map_err(|_| WritePdbError::CouldNotWrite())
}

/// Writer a header for a PDB file.
///
/// ## Warning
/// Currently only supports orthogonal simulation boxes.
fn write_header(
    writer: &mut BufWriter<File>,
    title: &str,
    simbox: &SimBox,
) -> Result<(), WritePdbError> {
    write_line(writer, &format!("TITLE     {}", title))?;

    write_line(writer, "REMARK    THIS IS A SIMULATION BOX")?;

    write_line(
        writer,
        &format!(
            "CRYST1{:>9.3}{:>9.3}{:>9.3}{:>7.2}{:>7.2}{:>7.2} P 1           1",
            simbox.x * 10.0,
            simbox.y * 10.0,
            simbox.z * 10.0,
            90.0,
            90.0,
            90.0
        ),
    )?;

    write_line(writer, "MODEL        1")?;

    Ok(())
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests_read {
    use super::*;
    use float_cmp::approx_eq;

    #[test]
    fn read_simple() {
        let system = read_pdb("test_files/example.pdb").unwrap();

        assert_eq!(system.get_name(), "Buforin II peptide P11L");
        assert_eq!(system.get_n_atoms(), 50);

        // check box size
        let simbox = system.get_box_as_ref();
        assert!(approx_eq!(f32, simbox.x, 6.0861));
        assert!(approx_eq!(f32, simbox.y, 6.0861));
        assert!(approx_eq!(f32, simbox.z, 6.0861));

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

        assert!(approx_eq!(f32, first.get_position().x, 1.660));
        assert!(approx_eq!(f32, first.get_position().y, 2.061));
        assert!(approx_eq!(f32, first.get_position().z, 3.153));

        // check atom somewhere in the middle
        let middle = &atoms[24];
        assert_eq!(middle.get_residue_number(), 11);
        assert_eq!(middle.get_residue_name(), "LEU");
        assert_eq!(middle.get_atom_name(), "SC1");
        assert_eq!(middle.get_atom_number(), 25);

        assert!(approx_eq!(f32, middle.get_position().x, 3.161));
        assert!(approx_eq!(f32, middle.get_position().y, 2.868));
        assert!(approx_eq!(f32, middle.get_position().z, 2.797));

        // check the last atom
        let last = &atoms[49];
        assert_eq!(last.get_residue_number(), 21);
        assert_eq!(last.get_residue_name(), "LYS");
        assert_eq!(last.get_atom_name(), "SC2");
        assert_eq!(last.get_atom_number(), 50);

        assert!(approx_eq!(f32, last.get_position().x, 4.706));
        assert!(approx_eq!(f32, last.get_position().y, 4.447));
        assert!(approx_eq!(f32, last.get_position().z, 2.813));

        // check that the velocity and force of all atoms is zero
        for atom in atoms.iter() {
            assert_eq!(atom.get_velocity().x, 0.0f32);
            assert_eq!(atom.get_velocity().y, 0.0f32);
            assert_eq!(atom.get_velocity().z, 0.0f32);

            assert_eq!(atom.get_force().x, 0.0f32);
            assert_eq!(atom.get_force().y, 0.0f32);
            assert_eq!(atom.get_force().z, 0.0f32);
        }
    }

    #[test]
    fn read_nochain() {
        let system_chain = read_pdb("test_files/example.pdb").unwrap();
        let system_nochain = read_pdb("test_files/example_nochain.pdb").unwrap();

        assert_eq!(system_chain.get_name(), system_nochain.get_name());
        assert_eq!(system_chain.get_box_as_ref().x, system_nochain.get_box_as_ref().x);
        assert_eq!(system_chain.get_box_as_ref().y, system_nochain.get_box_as_ref().y);
        assert_eq!(system_chain.get_box_as_ref().z, system_nochain.get_box_as_ref().z);

        for (ac, anc) in system_chain.atoms_iter().zip(system_nochain.atoms_iter()) {
            assert_eq!(ac.get_residue_number(), anc.get_residue_number());
            assert_eq!(ac.get_residue_name(), anc.get_residue_name());
            assert_eq!(ac.get_atom_number(), anc.get_atom_number());
            assert_eq!(ac.get_atom_name(), anc.get_atom_name());
            assert_eq!(ac.get_position().x, anc.get_position().x);
            assert_eq!(ac.get_position().y, anc.get_position().y);
            assert_eq!(ac.get_position().z, anc.get_position().z);

            assert_eq!(ac.get_velocity(), anc.get_velocity());
            assert_eq!(ac.get_force(), anc.get_force());
        }
    }

    #[test]
    fn read_hetatm() {
        let system = read_pdb("test_files/example_hetatm.pdb").unwrap();

        assert_eq!(system.get_name(), "Buforin II peptide P11L");
        assert_eq!(system.get_n_atoms(), 50);

        // check box size
        let simbox = system.get_box_as_ref();
        assert!(approx_eq!(f32, simbox.x, 6.0861));
        assert!(approx_eq!(f32, simbox.y, 6.0861));
        assert!(approx_eq!(f32, simbox.z, 6.0861));

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

        assert!(approx_eq!(f32, first.get_position().x, 1.660));
        assert!(approx_eq!(f32, first.get_position().y, 2.061));
        assert!(approx_eq!(f32, first.get_position().z, 3.153));

        // check atom somewhere in the middle
        let middle = &atoms[24];
        assert_eq!(middle.get_residue_number(), 11);
        assert_eq!(middle.get_residue_name(), "LEU");
        assert_eq!(middle.get_atom_name(), "SC1");
        assert_eq!(middle.get_atom_number(), 25);

        assert!(approx_eq!(f32, middle.get_position().x, 3.161));
        assert!(approx_eq!(f32, middle.get_position().y, 2.868));
        assert!(approx_eq!(f32, middle.get_position().z, 2.797));

        // check the last atom
        let last = &atoms[49];
        assert_eq!(last.get_residue_number(), 21);
        assert_eq!(last.get_residue_name(), "LYS");
        assert_eq!(last.get_atom_name(), "SC2");
        assert_eq!(last.get_atom_number(), 50);

        assert!(approx_eq!(f32, last.get_position().x, 4.706));
        assert!(approx_eq!(f32, last.get_position().y, 4.447));
        assert!(approx_eq!(f32, last.get_position().z, 2.813));

        // check that the velocity and force of all atoms is zero
        for atom in atoms.iter() {
            assert_eq!(atom.get_velocity().x, 0.0f32);
            assert_eq!(atom.get_velocity().y, 0.0f32);
            assert_eq!(atom.get_velocity().z, 0.0f32);

            assert_eq!(atom.get_force().x, 0.0f32);
            assert_eq!(atom.get_force().y, 0.0f32);
            assert_eq!(atom.get_force().z, 0.0f32);
        }
    }

    #[test]
    fn read_no_title() {
        let system = read_pdb("test_files/example_notitle.pdb").unwrap();

        assert_eq!(system.get_name(), "Unknown");
        assert_eq!(system.get_n_atoms(), 50);

        // check box size
        let simbox = system.get_box_as_ref();
        assert!(approx_eq!(f32, simbox.x, 6.0861));
        assert!(approx_eq!(f32, simbox.y, 6.0861));
        assert!(approx_eq!(f32, simbox.z, 6.0861));
    }

    #[test]
    fn read_empty_title() {
        let system = read_pdb("test_files/example_empty_title.pdb").unwrap();

        assert_eq!(system.get_name(), "Unknown");
        assert_eq!(system.get_n_atoms(), 50);

        // check box size
        let simbox = system.get_box_as_ref();
        assert!(approx_eq!(f32, simbox.x, 6.0861));
        assert!(approx_eq!(f32, simbox.y, 6.0861));
        assert!(approx_eq!(f32, simbox.z, 6.0861));
    }

    #[test]
    fn read_no_box() {
        let system = read_pdb("test_files/example_nobox.pdb").unwrap();

        assert_eq!(system.get_name(), "Buforin II peptide P11L");
        assert_eq!(system.get_n_atoms(), 50);

        // check box size
        let simbox = system.get_box_as_ref();
        assert!(approx_eq!(f32, simbox.x, 0.0));
        assert!(approx_eq!(f32, simbox.y, 0.0));
        assert!(approx_eq!(f32, simbox.z, 0.0));
    }

    #[test]
    fn read_multiple_titles() {
        let system = read_pdb("test_files/example_multiple_titles.pdb").unwrap();

        assert_eq!(system.get_name(), "Third title");
        assert_eq!(system.get_n_atoms(), 50);

        // check box size
        let simbox = system.get_box_as_ref();
        assert!(approx_eq!(f32, simbox.x, 6.0861));
        assert!(approx_eq!(f32, simbox.y, 6.0861));
        assert!(approx_eq!(f32, simbox.z, 6.0861));
    }

    #[test]
    fn read_multiple_boxes() {
        let system = read_pdb("test_files/example_multiple_boxes.pdb").unwrap();

        assert_eq!(system.get_name(), "Buforin II peptide P11L");
        assert_eq!(system.get_n_atoms(), 50);

        // check box size
        let simbox = system.get_box_as_ref();
        assert!(approx_eq!(f32, simbox.x, 5.0861));
        assert!(approx_eq!(f32, simbox.y, 5.0861));
        assert!(approx_eq!(f32, simbox.z, 5.0861));
    }

    macro_rules! read_pdb_fails {
        ($name:ident, $file:expr, $variant:path, $expected:expr) => {
            #[test]
            fn $name() {
                let file = $file;
                match read_pdb(file) {
                    Err($variant(e)) => assert_eq!(e, $expected),
                    Ok(_) => panic!("Parsing should have failed, but it succeeded."),
                    Err(e) => panic!("Parsing successfully failed but incorrect error type `{:?}` was returned.", e),
                }
            }
        };
    }

    read_pdb_fails!(
        read_nonorthogonal_box,
        "test_files/example_nonorthogonal.pdb",
        ParsePdbError::NonOrthogonalBox,
        "CRYST1   60.861   60.861   60.861  90.00  89.00  90.00 P 1           1"
    );

    read_pdb_fails!(
        read_invalid_box,
        "test_files/example_invalid_box.pdb",
        ParsePdbError::ParseBoxLineErr,
        "CRYST1   60.861   60.f61   60.861  90.00  90.00  90.00 P 1           1"
    );

    read_pdb_fails!(
        read_invalid_box2,
        "test_files/example_invalid_box2.pdb",
        ParsePdbError::ParseBoxLineErr,
        "CRYST1   60.861   60.861   60.861  90.00  90.00  90.0O P 1           1"
    );

    read_pdb_fails!(
        read_short_box,
        "test_files/example_short_box.pdb",
        ParsePdbError::ParseBoxLineErr,
        "CRYST1   60.861   60.861   60.861  90.00  90.00  90.0"
    );

    read_pdb_fails!(
        read_short_atom,
        "test_files/example_short_atom.pdb",
        ParsePdbError::ParseAtomLineErr,
        "ATOM     35  SC1 HIS A  16      34.580  36.530  27.8"
    );

    read_pdb_fails!(
        read_invalid_atom,
        "test_files/example_invalid_atom.pdb",
        ParsePdbError::ParseAtomLineErr,
        "ATOM      30  SC1 ARG A  14      32.540  35.200  34.040  1.00  0.00            "
    );
}

#[cfg(test)]
mod tests_write {
    use super::*;
    use file_diff;
    use tempfile::NamedTempFile;

    #[test]
    fn write() {
        let system = System::from_file("test_files/example_novelocities.gro").unwrap();

        let pdb_output = NamedTempFile::new().unwrap();
        let path_to_output = pdb_output.path();

        if let Err(_) = system.write_pdb(path_to_output) {
            panic!("Writing pdb file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example_nochain.pdb").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_fails() {
        let system = System::from_file("test_files/example.gro").unwrap();

        match system.write_pdb("Xhfguiaghqueiowhd/nonexistent.ndx") {
            Err(WritePdbError::CouldNotCreate(e)) => {
                assert_eq!(e, Box::from(Path::new("Xhfguiaghqueiowhd/nonexistent.ndx")))
            }
            Ok(_) => panic!("Writing should have failed, but it did not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }
    }

    #[test]
    fn write_wrap() {
        let atom1 = Atom::new(
            158,
            "THR",
            1,
            "BBBBT",
            [0.0, 0.0, 0.0].into(),
            [0.0, 0.0, 0.0].into(),
            [0.0, 0.0, 0.0].into(),
        );

        let atom2 = Atom::new(
            158,
            "THR",
            99999,
            "SC1",
            [0.0, 0.0, 0.0].into(),
            [0.0, 0.0, 0.0].into(),
            [0.0, 0.0, 0.0].into(),
        );

        let atom3 = Atom::new(
            10003,
            "ARG",
            100000,
            "BB",
            [0.0, 0.0, 0.0].into(),
            [0.0, 0.0, 0.0].into(),
            [0.0, 0.0, 0.0].into(),
        );

        let atom4 = Atom::new(
            10003,
            "ARGGT",
            200001,
            "SC1",
            [0.0, 0.0, 0.0].into(),
            [0.0, 0.0, 0.0].into(),
            [0.0, 0.0, 0.0].into(),
        );

        let atom5 = Atom::new(
            10003,
            "ARG",
            200005,
            "SC2",
            [0.0, 0.0, 0.0].into(),
            [0.0, 0.0, 0.0].into(),
            [0.0, 0.0, 0.0].into(),
        );

        let atoms = vec![atom1, atom2, atom3, atom4, atom5];
        let simbox = SimBox::from([1.0, 1.0, 1.0]);

        let system = System::new("Expected atom and residue wrapping", atoms, simbox);

        let pdb_output = NamedTempFile::new().unwrap();
        let path_to_output = pdb_output.path();

        if let Err(_) = system.write_pdb(path_to_output) {
            panic!("Writing pdb file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/wrapping_expected.pdb").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }
    
    #[test]
    fn write_group() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.read_ndx("test_files/index.ndx").unwrap();

        let pdb_output = NamedTempFile::new().unwrap();
        let path_to_output = pdb_output.path();

        if let Err(_) = system.group_write_pdb("Protein", path_to_output) {
            panic!("Writing pdb file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/protein.pdb").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group_fails() {
        let system = System::from_file("test_files/example.gro").unwrap();

        let pdb_output = NamedTempFile::new().unwrap();
        let path_to_output = pdb_output.path();

        match system.group_write_pdb("Protein", path_to_output) {
            Err(WritePdbError::GroupNotFound(e)) => assert_eq!(e, "Protein"),
            Ok(_) => panic!("Writing should have failed, but it did not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }
    }
}