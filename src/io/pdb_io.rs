// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of functions for reading and writing pdb files.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::errors::{ParsePdbConnectivityError, ParsePdbError, WritePdbError};
use crate::structures::{atom::Atom, simbox::SimBox, vector3d::Vector3D};
use crate::system::general::System;

/// Read a pdb file and construct a System structure. Does not read connectivity.
///
/// ## Warning
/// Currently only supports orthogonal simulation boxes!
///
/// ## Supported keywords
/// This function can handle lines starting with ATOM, HETATM, TITLE, ENDMDL, END, and CRYST1.
/// All other lines are ignored.
///
/// ## Notes
/// - Reading ends once `ENDMDL`, `END`, or the end of file is reached.
///
/// - In case multiple TITLE lines are provided, the **last one** is used as the
/// name of the system. If no TITLE line is provided, "Unknown" is used as the name.
///
/// - In case multiple CRYST1 lines are provided, information from the **last one** is used.
/// If no CRYST1 line is provided, the simulation box size is set to 0 in all dimensions.
/// Note again that only orthogonal simulation boxes are supported.
///
/// - If you want to load connectivity from the PDB file,
/// use `System::add_bonds_from_pdb` after constructing the `System` structure.
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
        if (line.len() >= 4 && line[0..4] == *"ATOM")
            || (line.len() >= 6 && line[0..6] == *"HETATM")
        {
            atoms.push(line_as_atom(&line)?);
        }
        // parse TITLE line
        else if line.len() >= 5 && line[0..5] == *"TITLE" {
            title = line_as_title(&line)?;
        }
        // parse CRYST1 line
        else if line.len() >= 6 && line[0..6] == *"CRYST1" {
            simbox = line_as_box(&line)?;
        }
        // END or ENDMDL is reached => stop reading
        else if line.len() >= 3 && line[0..3] == *"END" {
            break;
        }
    }

    Ok(System::new(&title, atoms, simbox))
}

/// ## Advanced methods for working with PDB files.
impl System {
    /// Read connectivity information from a PDB file.
    ///
    /// ## Returns
    /// - `Ok` if the connectivity was read correctly.
    /// - `ParsePdbConnectivityError::NoBondsWarning` is the file was read successfully
    /// but no bonds have been found. The `System` structure is still modified.
    /// - Other `ParsePdbConnectivityError` if an error occured.
    /// If an error occurs, the `System` structure is not changed.
    ///
    /// ## Example
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// // read structure from a pdb file
    /// let mut system = System::from_file("system.pdb").unwrap();
    ///
    /// // add connectivity information from a different pdb file
    /// if let Err(e) = system.add_bonds_from_pdb("connectivity.pdb") {
    ///     eprintln!("{}", e);
    ///     return;
    /// }
    ///
    /// // perform analysis...
    /// ```
    ///
    /// ## Warning
    /// The connectivity block in the PDB file uses atom numbers specified in the PDB file.
    /// If you renumber the atoms of your system after loading the structure and BEFORE
    /// adding the bonds from the PDB file, you may get completely nonsensical result!
    ///
    /// In other words, **don't** do this:
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.pdb").unwrap();
    /// system.atoms_renumber(); // <--- DON'T!
    ///
    /// system.add_bonds_from_pdb("system.pdb").unwrap();
    /// ```
    ///
    /// ## Notes
    /// - Reading ends once `END`` keyword or the end of file is reached.
    /// - Unlike with `read_pdb`, the `ENDMDL` keyword is ignored by this function.
    /// - Note that the PDB file with the connectivity information can but does not have to contain information about the atoms.
    /// All lines other than those starting with the `CONECT` keyword are ignored.
    /// - This function can read `CONECT` lines of any length, i.e. it is not limited by the traditional requirement
    /// of using at most 4 bonds in a single `CONECT` line.
    pub fn add_bonds_from_pdb(
        &mut self,
        filename: impl AsRef<Path>,
    ) -> Result<(), ParsePdbConnectivityError> {
        // sanity checking the system
        // if the system contains multiple atoms with the same atom number, the connectivity can't be used
        // this is because atom numbers are used in the CONECT lines:
        // if there are duplicates, it is impossible to decide what atoms are bonded to what atoms
        if self.has_duplicate_atom_numbers() {
            return Err(ParsePdbConnectivityError::DuplicateAtomNumbers);
        }

        // open the pdb file
        let file = match File::open(filename.as_ref()) {
            Ok(x) => x,
            Err(_) => {
                return Err(ParsePdbConnectivityError::FileNotFound(Box::from(
                    filename.as_ref(),
                )))
            }
        };

        let reader = BufReader::new(file);

        // create a mapping from atom numbers to their indices for quick lookup
        let atom_number_to_index: HashMap<usize, usize> = self
            .atoms_iter()
            .enumerate()
            .map(|(index, atom)| (atom.get_atom_number(), index))
            .collect();

        // create vector for temporary storing of connectivity information
        let mut temp_bonded = vec![Vec::new(); self.get_n_atoms()];

        for raw_line in reader.lines() {
            let line = match raw_line {
                Ok(x) => x,
                Err(_) => {
                    return Err(ParsePdbConnectivityError::LineNotFound(Box::from(
                        filename.as_ref(),
                    )))
                }
            };

            if line.len() >= 6 && line[0..6] == *"CONECT" {
                line_as_conect(&mut temp_bonded, &line, &atom_number_to_index)?;
            } else if line.trim().len() == 3 && line[0..3] == *"END" {
                break;
            }
        }

        // validate that all bonded atoms agree that they are bonded
        // i.e. all bonds must be bidirectional edges of the graph describing the system topology
        for atom_index in 1..self.get_n_atoms() {
            let bonded1 = temp_bonded
                .get(atom_index)
                .expect("FATAL GROAN ERROR | System::add_bonds_from_pdb (1) | Invalid atom index.");

            for &bonded_index in bonded1.iter() {
                let bonded2 = temp_bonded.get(bonded_index).expect(
                    "FATAL GROAN ERROR | System::add_bonds_from_pdb (2) | Invalid atom index.",
                );

                if !bonded2.contains(&atom_index) {
                    let atom1 = self.get_atom_as_ref(atom_index).expect(
                        "FATAL GROAN ERROR | System::add_bonds_from_pdb (3) | Invalid atom index.",
                    );
                    let atom2 = self.get_atom_as_ref(bonded_index).expect(
                        "FATAL GROAN ERROR | System::add_bonds_from_pdb (4) | Invalid atom index.",
                    );
                    return Err(ParsePdbConnectivityError::BondingInconsistency(
                        atom1.get_atom_number(),
                        atom2.get_atom_number(),
                    ));
                }
            }
        }

        let mut empty = true;
        // transform `temp_bonded` to information in the `System` structure
        for (bonded, atom) in temp_bonded.into_iter().zip(self.atoms_iter_mut()) {
            // safety:
            // a) we know that all indices are valid as we obtain them from `atom_number_to_index` hashmap
            // b) we apply the `set_bonded` method to all atoms
            unsafe {
                if !bonded.is_empty() {
                    empty = false;
                }

                atom.set_bonded(bonded);
            }
        }

        if empty {
            // no bonds have been detected in the pdb file
            Err(ParsePdbConnectivityError::NoBondsWarning(Box::from(
                filename.as_ref(),
            )))
        } else {
            Ok(())
        }
    }
}

/// ## Methods for writing pdb files.
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
    /// use groan_rs::prelude::*;
    ///
    /// let system = System::from_file("system.gro").unwrap();
    /// if let Err(e) = system.write_pdb("system.pdb", false) {
    ///     eprintln!("{}", e);
    ///     return;
    /// }
    /// ```
    ///
    /// ## Notes on connectivity
    /// - The function will attempt to write connectivity into the PDB file if `write_connectivity` is true.
    /// - Connectivity block can only be written for systems with fewer than 100,000 atoms
    /// and for systems in which each atom has a unique number. In case these requirements
    /// are not fulfilled and the `write_connectivity` is `true`, an error is returned
    /// and no output PDB file is written.
    /// - Even though `groan_rs` library can read `CONECT` lines of any length,
    /// this function prints at most 4 bonds on a single `CONECT` line as is traditionally requested.
    /// If the atom is bonded to more atoms, multiple `CONECT` lines will be written for it.
    pub fn write_pdb(
        &self,
        filename: impl AsRef<Path>,
        write_connectivity: bool,
    ) -> Result<(), WritePdbError> {
        match self.group_write_pdb("all", filename, write_connectivity) {
            Ok(_) => Ok(()),
            Err(WritePdbError::GroupNotFound(_)) => {
                panic!(
                    "FATAL GROAN ERROR | System::write_pdb | Default group 'all' does not exist."
                )
            }
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
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// if let Err(e) = system.group_write_pdb("Protein", "protein.pdb", false) {
    ///     eprintln!("{}", e);
    ///     return;
    /// }
    /// ```
    ///
    /// ## Notes on connectivity
    /// - The function will attempt to write connectivity into the PDB file if `write_connectivity` is true.
    /// - Connectivity block can only be written for systems with fewer than 100,000 atoms
    /// and for systems in which each atom has a unique number. In case these requirements
    /// are not fulfilled and the `write_connectivity` is `true`, an error is returned
    /// and no output PDB file is written.
    /// - Even though `groan_rs` library can read `CONECT` lines of any length,
    /// this function prints at most 4 bonds on a single `CONECT` line as is traditionally requested.
    /// If the atom is bonded to more atoms, multiple `CONECT` lines will be written for it.
    pub fn group_write_pdb(
        &self,
        group_name: &str,
        filename: impl AsRef<Path>,
        write_connectivity: bool,
    ) -> Result<(), WritePdbError> {
        if !self.group_exists(group_name) {
            return Err(WritePdbError::GroupNotFound(group_name.to_string()));
        }

        if write_connectivity {
            // connectivity block can't be printed if the number of atoms in the system is higher than 99,999
            if self.get_n_atoms() > 99_999 {
                return Err(WritePdbError::ConectTooLarge(self.get_n_atoms()));
            }
            // or if any of the atoms in the structure have duplicate numbers
            if self.has_duplicate_atom_numbers() {
                return Err(WritePdbError::ConectDuplicateAtomNumbers);
            }
        }

        let output = File::create(&filename)
            .map_err(|_| WritePdbError::CouldNotCreate(Box::from(filename.as_ref())))?;

        let mut writer = BufWriter::new(output);

        let title = match group_name {
            "all" => self.get_name().to_owned(),
            _ => format!("Group `{}` from {}", group_name, self.get_name()),
        };

        write_header(&mut writer, &title, self.get_box_as_ref())?;

        for atom in self.group_iter(group_name).expect(
            "FATAL GROAN ERROR | System::group_write_pdb | Group should exist but it does not.",
        ) {
            atom.write_pdb(&mut writer)?;
        }

        write_line(&mut writer, "TER\nENDMDL")?;

        if write_connectivity {
            write_connectivity_section(self, &mut writer, group_name)?;
        }

        write_line(&mut writer, "END")?;

        writer.flush().map_err(|_| WritePdbError::CouldNotWrite)?;

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

    // parsing chain
    let chain = line.chars().nth(21).filter(|&x| !x.is_whitespace());

    // parsing residue number
    let residue_number = line[22..26]
        .trim()
        .parse::<usize>()
        .map_err(|_| ParsePdbError::ParseAtomLineErr(line.to_string()))?;

    // parsing position
    let mut curr = 30usize;
    let mut position = [0.0, 0.0, 0.0];
    for pos in &mut position {
        *pos = line[curr..curr + 8]
            .trim()
            .parse::<f32>()
            .map(|x| x / 10.0)
            .map_err(|_| ParsePdbError::ParseAtomLineErr(line.to_string()))?;

        curr += 8;
    }

    let atom = Atom::new(
        residue_number,
        &residue_name,
        atom_number,
        &atom_name,
        position.into(),
        Vector3D::default(),
        Vector3D::default(),
    );

    // add chain information, if available
    match chain {
        Some(x) => Ok(atom.with_chain(x)),
        None => Ok(atom),
    }
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
    for dim in &mut boxsize {
        *dim = line[curr..curr + 9]
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

/// Parse a single line as connectivity information.
///
/// ## Notes
/// - Parses a line starting with CONECT.
/// - Can read CONECT line of any length.
fn line_as_conect(
    temp_bonded: &mut [Vec<usize>],
    line: &str,
    number2index: &HashMap<usize, usize>,
) -> Result<(), ParsePdbConnectivityError> {
    if line.len() < 11 {
        return Err(ParsePdbConnectivityError::ParseConectLineErr(
            line.to_string(),
        ));
    }

    // get atom number of the target atom
    let atom_number = line[6..11]
        .trim()
        .parse::<usize>()
        .map_err(|_| ParsePdbConnectivityError::ParseConectLineErr(line.to_string()))?;

    // get index of the target atom
    let atom_index = match number2index.get(&atom_number) {
        Some(i) => i,
        None => {
            return Err(ParsePdbConnectivityError::AtomNotFound(
                atom_number,
                line.to_string(),
            ))
        }
    };

    // parse atom numbers on the line
    let mut iterator = 11usize;
    while iterator + 4 < line.len() {
        let trimmed = line[iterator..(iterator + 5)].trim();

        if !trimmed.is_empty() {
            let number = trimmed
                .parse::<usize>()
                .map_err(|_| ParsePdbConnectivityError::ParseConectLineErr(line.to_string()))?;

            // convert atom number of the bonded atom to its index in the system structure
            let index = match number2index.get(&number) {
                Some(i) => i,
                None => {
                    return Err(ParsePdbConnectivityError::AtomNotFound(
                        number,
                        line.to_string(),
                    ))
                }
            };

            // check that the target atom is not the same atom as the bonded atom
            if atom_index == index {
                return Err(ParsePdbConnectivityError::SelfBonding(atom_number));
            }

            // storage the index of the bonded atom
            temp_bonded[*atom_index].push(*index);
        }

        iterator += 5;
    }

    Ok(())
}

fn write_line<W: Write>(writer: &mut W, line: &str) -> Result<(), WritePdbError> {
    writeln!(writer, "{}", line).map_err(|_| WritePdbError::CouldNotWrite)
}

fn write<W: Write>(writer: &mut W, string: &str) -> Result<(), WritePdbError> {
    write!(writer, "{}", string).map_err(|_| WritePdbError::CouldNotWrite)
}

/// Write a header for a PDB file.
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

/// Write connectivity section for the system.
fn write_connectivity_section(
    system: &System,
    writer: &mut BufWriter<File>,
    group_name: &str,
) -> Result<(), WritePdbError> {
    for atom in system
        .group_iter(group_name)
        .expect("FATAL GROAN ERROR | pdb_io::write_connectivity_section (1) | Group should exist but it does not.")
    {
        // select atoms from bonded which are inside the printed group
        let bonded_in_group = match atom.get_bonded() {
            None => return Ok(()),
            Some(container) => {
                container.iter().filter(|index| system
                    .group_isin(group_name, *index)
                    .expect("FATAL GROAN ERROR | pdb_io::write_connectivity_section (2) | Group should exist but it does not.")
                )
                .collect::<Vec<usize>>()
            }
        };

        let n_bonded = bonded_in_group.len();
        if n_bonded == 0 {
            continue;
        }

        let atom_number = atom.get_atom_number();
        // connectivity section can not accommodate larger than 5-digit numbers
        if atom_number > 99_999 {
            return Err(WritePdbError::ConectInvalidNumber(atom_number));
        }

        for (i, bonded_index) in bonded_in_group.into_iter().enumerate() {

            // traditionally, PDB file only supports at most 4 bonds for a single atom on a single line
            // higher number of bonds must be separated into multiple CONECT lines
            if i % 4 == 0 {
                write(writer, &format!("CONECT{:>5}", atom_number))?;
            }

            let bonded_number = system
                .get_atom_as_ref(bonded_index)
                .expect("FATAL GROAN ERROR | pdb_io::write_connectivity_section | Invalid atom index.")
                .get_atom_number();

                // connectivity section can not accomodate larger than 5-digit numbers
                if bonded_number > 99_999 {
                    return Err(WritePdbError::ConectInvalidNumber(bonded_number));
                }

                // write bonded atom number to the output file
                write(writer, &format!("{:>5}", bonded_number))?;

            if i % 4 == 3 || i == n_bonded - 1 {
                write_line(writer, "")?;
            }
        }
    }

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
        assert_eq!(first.get_chain().unwrap(), 'A');

        assert!(approx_eq!(f32, first.get_position().x, 1.660));
        assert!(approx_eq!(f32, first.get_position().y, 2.061));
        assert!(approx_eq!(f32, first.get_position().z, 3.153));

        // check atom somewhere in the middle
        let middle = &atoms[24];
        assert_eq!(middle.get_residue_number(), 11);
        assert_eq!(middle.get_residue_name(), "LEU");
        assert_eq!(middle.get_atom_name(), "SC1");
        assert_eq!(middle.get_atom_number(), 25);
        assert_eq!(middle.get_chain().unwrap(), 'B');

        assert!(approx_eq!(f32, middle.get_position().x, 3.161));
        assert!(approx_eq!(f32, middle.get_position().y, 2.868));
        assert!(approx_eq!(f32, middle.get_position().z, 2.797));

        // check the last atom
        let last = &atoms[49];
        assert_eq!(last.get_residue_number(), 21);
        assert_eq!(last.get_residue_name(), "LYS");
        assert_eq!(last.get_atom_name(), "SC2");
        assert_eq!(last.get_atom_number(), 50);
        assert_eq!(last.get_chain().unwrap(), 'C');

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
    fn read_endmdl() {
        let system = read_pdb("test_files/example_endmdl.pdb").unwrap();

        assert_eq!(system.get_name(), "Buforin II peptide P11L");
        assert_eq!(system.get_n_atoms(), 17);

        assert_eq!(system.atoms_iter().nth(0).unwrap().get_atom_number(), 1);
        assert_eq!(system.atoms_iter().nth(16).unwrap().get_atom_number(), 17);
    }

    #[test]
    fn read_end() {
        let system = read_pdb("test_files/example_end.pdb").unwrap();

        assert_eq!(system.get_name(), "Buforin II peptide P11L");
        assert_eq!(system.get_n_atoms(), 17);

        assert_eq!(system.atoms_iter().nth(0).unwrap().get_atom_number(), 1);
        assert_eq!(system.atoms_iter().nth(16).unwrap().get_atom_number(), 17);
    }

    #[test]
    fn read_nochain() {
        let system_chain = read_pdb("test_files/example.pdb").unwrap();
        let system_nochain = read_pdb("test_files/example_nochain.pdb").unwrap();

        assert_eq!(system_chain.get_name(), system_nochain.get_name());
        assert_eq!(
            system_chain.get_box_as_ref().x,
            system_nochain.get_box_as_ref().x
        );
        assert_eq!(
            system_chain.get_box_as_ref().y,
            system_nochain.get_box_as_ref().y
        );
        assert_eq!(
            system_chain.get_box_as_ref().z,
            system_nochain.get_box_as_ref().z
        );

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

            assert_eq!(anc.get_chain(), None);
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
        assert_eq!(first.get_chain().unwrap(), 'A');

        assert!(approx_eq!(f32, first.get_position().x, 1.660));
        assert!(approx_eq!(f32, first.get_position().y, 2.061));
        assert!(approx_eq!(f32, first.get_position().z, 3.153));

        // check atom somewhere in the middle
        let middle = &atoms[24];
        assert_eq!(middle.get_residue_number(), 11);
        assert_eq!(middle.get_residue_name(), "LEU");
        assert_eq!(middle.get_atom_name(), "SC1");
        assert_eq!(middle.get_atom_number(), 25);
        assert_eq!(middle.get_chain().unwrap(), 'A');

        assert!(approx_eq!(f32, middle.get_position().x, 3.161));
        assert!(approx_eq!(f32, middle.get_position().y, 2.868));
        assert!(approx_eq!(f32, middle.get_position().z, 2.797));

        // check the last atom
        let last = &atoms[49];
        assert_eq!(last.get_residue_number(), 21);
        assert_eq!(last.get_residue_name(), "LYS");
        assert_eq!(last.get_atom_name(), "SC2");
        assert_eq!(last.get_atom_number(), 50);
        assert_eq!(last.get_chain().unwrap(), 'A');

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

    #[test]
    fn add_bonds_from_pdb() {
        let mut system = read_pdb("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        let expected_bonded: Vec<Vec<usize>> = vec![
            vec![3, 2],
            vec![1],
            vec![4, 1, 6],
            vec![3, 5],
            vec![4],
            vec![7, 3, 8],
            vec![6],
            vec![9, 6, 10],
            vec![8],
            vec![10, 12],
            vec![11],
            vec![10, 15, 14],
            vec![13],
            vec![13, 16],
            vec![17, 15, 18],
            vec![16],
            vec![16, 20, 19],
            vec![18],
            vec![21, 18, 24],
            vec![20, 22, 23],
            vec![21, 23],
            vec![21, 22],
            vec![25, 20, 26],
            vec![24],
            vec![24, 28, 27],
            vec![26],
            vec![26, 29],
            vec![11, 8, 13],
            vec![30, 28, 32, 36, 38, 42, 48],
            vec![29, 31],
            vec![30],
            vec![29, 34, 33],
            vec![32],
            vec![35, 32, 38],
            vec![34, 36, 37],
            vec![35, 37, 29],
            vec![35, 36],
            vec![39, 34, 41, 29],
            vec![38, 40],
            vec![39],
            vec![42, 38, 43],
            vec![41, 29],
            vec![44, 41, 45],
            vec![43],
            vec![46, 43, 48],
            vec![45, 47],
            vec![46],
            vec![49, 45, 29],
            vec![48],
            vec![],
        ];

        for (a, atom) in system.atoms_iter().enumerate() {
            let expected = expected_bonded.get(a).unwrap();
            for index in atom.get_bonded().unwrap().iter() {
                let bonded = system
                    .get_atoms_as_ref()
                    .get(index)
                    .unwrap()
                    .get_atom_number();
                assert!(expected.contains(&bonded));
            }
        }

        assert!(system.get_atom_as_ref(49).unwrap().get_bonded().is_some());
    }

    #[test]
    fn add_bonds_from_pdb_2() {
        let mut system1 = read_pdb("test_files/conect.pdb").unwrap();
        system1.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        let mut system2 = read_pdb("test_files/conect.pdb").unwrap();
        system2
            .add_bonds_from_pdb("test_files/bonds_for_example.pdb")
            .unwrap();

        for (atom1, atom2) in system1.atoms_iter().zip(system2.atoms_iter()) {
            assert_eq!(atom1.get_bonded(), atom2.get_bonded());
        }
    }

    #[test]
    fn add_bonds_empty_pdb() {
        let mut system = read_pdb("test_files/example.pdb").unwrap();

        match system.add_bonds_from_pdb("test_files/example.pdb") {
            Err(ParsePdbConnectivityError::NoBondsWarning(path)) => {
                assert_eq!(path, Box::from(Path::new("test_files/example.pdb")));
                assert!(system.has_bonds());
            }
            Ok(_) => {
                panic!("Parsing should have returned a warning but it succeeded without warning.")
            }
            Err(e) => panic!(
                "Parsing failed with an error `{:?}` instead of a warning.",
                e
            ),
        }
    }

    #[test]
    fn add_bonds_end() {
        let mut system = read_pdb("test_files/conect.pdb").unwrap();

        match system.add_bonds_from_pdb("test_files/conect_end.pdb") {
            Err(ParsePdbConnectivityError::NoBondsWarning(path)) => {
                assert_eq!(path, Box::from(Path::new("test_files/conect_end.pdb")));
                assert!(system.has_bonds());
            }
            Ok(_) => {
                panic!("Parsing should have returned a warning but it succeeded without warning.")
            }
            Err(e) => panic!(
                "Parsing failed with an error `{:?}` instead of a warning.",
                e
            ),
        }
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

    macro_rules! read_bonds_fails {
        ($name:ident, $file_struct:expr, $file_bonds:expr, $variant:path, $expected:expr) => {
            #[test]
            fn $name() {
                let mut system = read_pdb($file_struct).unwrap();

                match system.add_bonds_from_pdb($file_bonds) {
                    Err($variant(e)) => {
                        assert_eq!(e, $expected);
                        for atom in system.get_atoms_as_ref() {
                            assert!(atom.get_bonded().is_none());
                        }
                    }
                    Ok(_) => panic!("Parsing should have failed, but it succeeded."),
                    Err(e) => panic!("Parsing successfully failed but incorrect error type `{:?}` was returned.", e),
                }
            }
        };
    }

    read_bonds_fails!(
        pdb_bonds_nonexistent,
        "test_files/example.pdb",
        "test_files/nonexistent.pdb",
        ParsePdbConnectivityError::FileNotFound,
        Box::from(Path::new("test_files/nonexistent.pdb"))
    );

    read_bonds_fails!(
        pdb_bonds_parse_error_1,
        "test_files/example.pdb",
        "test_files/bonds_parse_error_1.pdb",
        ParsePdbConnectivityError::ParseConectLineErr,
        "CONECT"
    );

    read_bonds_fails!(
        pdb_bonds_parse_error_2,
        "test_files/example.pdb",
        "test_files/bonds_parse_error_2.pdb",
        ParsePdbConnectivityError::ParseConectLineErr,
        "CONECT   43   4A   41   45                                            "
    );

    #[test]
    fn pdb_bonds_invalid_index_1() {
        let mut system = read_pdb("test_files/example.pdb").unwrap();
        match system.add_bonds_from_pdb("test_files/bonds_invalid_index_1.pdb") {
            Err(ParsePdbConnectivityError::AtomNotFound(index, string)) => {
                assert_eq!(index, 51);
                assert_eq!(string, "CONECT    3    4    1   51");
            }
            Ok(_) => panic!("Parsing should have failed, but it succeeded."),
            Err(e) => panic!(
                "Parsing successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn pdb_bonds_invalid_index_2() {
        let mut system = read_pdb("test_files/example.pdb").unwrap();
        match system.add_bonds_from_pdb("test_files/bonds_invalid_index_2.pdb") {
            Err(ParsePdbConnectivityError::AtomNotFound(index, string)) => {
                assert_eq!(index, 55);
                assert_eq!(
                    string,
                    "CONECT   55   35   37   29                                            "
                );
            }
            Ok(_) => panic!("Parsing should have failed, but it succeeded."),
            Err(e) => panic!(
                "Parsing successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn pdb_bonds_inconsistent() {
        let mut system = read_pdb("test_files/example.pdb").unwrap();
        match system.add_bonds_from_pdb("test_files/bonds_inconsistency.pdb") {
            Err(ParsePdbConnectivityError::BondingInconsistency(atom1, atom2)) => {
                assert_eq!(atom1, 18);
                assert_eq!(atom2, 16);
            }
            Ok(_) => panic!("Parsing should have failed, but it succeeded."),
            Err(e) => panic!(
                "Parsing successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    read_bonds_fails!(
        pdb_bonds_selfbonding,
        "test_files/example.pdb",
        "test_files/bonds_selfbonding.pdb",
        ParsePdbConnectivityError::SelfBonding,
        44
    );

    #[test]
    fn pdb_bonds_duplicate_numbers() {
        let mut system = read_pdb("test_files/example.pdb").unwrap();
        system.get_atom_as_ref_mut(10).unwrap().set_atom_number(25);

        match system.add_bonds_from_pdb("test_files/bonds_inconsistency.pdb") {
            Err(ParsePdbConnectivityError::DuplicateAtomNumbers) => (),
            Ok(_) => panic!("Parsing should have failed, but it succeeded."),
            Err(e) => panic!(
                "Parsing successfully failed but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }
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

        if let Err(_) = system.write_pdb(path_to_output, false) {
            panic!("Writing pdb file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example_nochain.pdb").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_fails() {
        let system = System::from_file("test_files/example.gro").unwrap();

        match system.write_pdb("Xhfguiaghqueiowhd/nonexistent.ndx", false) {
            Err(WritePdbError::CouldNotCreate(e)) => {
                assert_eq!(e, Box::from(Path::new("Xhfguiaghqueiowhd/nonexistent.ndx")))
            }
            Ok(_) => panic!("Writing should have failed, but it did not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }
    }

    #[test]
    fn write_with_chains() {
        let system = System::from_file("test_files/example.pdb").unwrap();

        let pdb_output = NamedTempFile::new().unwrap();
        let path_to_output = pdb_output.path();

        if let Err(_) = system.write_pdb(path_to_output, false) {
            panic!("Writing pdb file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example.pdb").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
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

        if let Err(_) = system.write_pdb(path_to_output, false) {
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

        if let Err(_) = system.group_write_pdb("Protein", path_to_output, false) {
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

        match system.group_write_pdb("Protein", path_to_output, false) {
            Err(WritePdbError::GroupNotFound(e)) => assert_eq!(e, "Protein"),
            Ok(_) => panic!("Writing should have failed, but it did not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }
    }

    #[test]
    fn write_with_connectivity() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        let pdb_output = NamedTempFile::new().unwrap();
        let path_to_output = pdb_output.path();

        system.write_pdb(path_to_output, true).unwrap();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/expected_bonds.pdb").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_with_connectivity_no_bonds() {
        let mut system = System::from_file("test_files/example.pdb").unwrap();
        match system.add_bonds_from_pdb("test_files/example.pdb") {
            Ok(_) | Err(ParsePdbConnectivityError::NoBondsWarning(_)) => (),
            Err(e) => panic!("Could not read bonds from file: `{}`", e),
        }

        let pdb_output = NamedTempFile::new().unwrap();
        let path_to_output = pdb_output.path();

        system.write_pdb(path_to_output, true).unwrap();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example.pdb").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group_with_connectivity() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        system.group_create("Group", "serial 20 to 30").unwrap();

        let pdb_output = NamedTempFile::new().unwrap();
        let path_to_output = pdb_output.path();

        system
            .group_write_pdb("Group", path_to_output, true)
            .unwrap();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/group_expected_bonds.pdb").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_with_connectivity_fails_too_large() {
        let mut atoms = Vec::with_capacity(100_000);
        let atom = Atom::new(1, "RES", 1, "ATM", Vector3D::default(), Vector3D::default(), Vector3D::default());

        for _ in 0..100_000 {
            atoms.push(atom.clone());
        }

        let system = System::new("Test system", atoms, [10.0, 10.0, 10.0].into());

        let pdb_output = NamedTempFile::new().unwrap();
        let path_to_output = pdb_output.path();

        match system.write_pdb(path_to_output, true) {
                Ok(_) => panic!("Writing should have failed but it succeeded."),
                Err(WritePdbError::ConectTooLarge(e)) => assert_eq!(e, system.get_n_atoms()),
                Err(e) => panic!("Writing successfully failed but incorrect error type `{:?}` was returned.", e),
        }
    }

    #[test]
    fn write_with_connectivity_fails_duplicate() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        system.get_atom_as_ref_mut(10).unwrap().set_atom_number(4);

        let pdb_output = NamedTempFile::new().unwrap();
        let path_to_output = pdb_output.path();

        match system.write_pdb(path_to_output, true) {
                Ok(_) => panic!("Writing should have failed but it succeeded."),
                Err(WritePdbError::ConectDuplicateAtomNumbers) => (),
                Err(e) => panic!("Writing successfully failed but incorrect error type `{:?}` was returned.", e),
        }
    }

    #[test]
    fn write_with_connectivity_fails_number_too_high() {
        let mut system = System::from_file("test_files/conect.pdb").unwrap();
        system.add_bonds_from_pdb("test_files/conect.pdb").unwrap();

        system.get_atom_as_ref_mut(10).unwrap().set_atom_number(100_000);

        let pdb_output = NamedTempFile::new().unwrap();
        let path_to_output = pdb_output.path();

        match system.write_pdb(path_to_output, true) {
                Ok(_) => panic!("Writing should have failed but it succeeded."),
                Err(WritePdbError::ConectInvalidNumber(e)) => assert_eq!(e, 100_000),
                Err(e) => panic!("Writing successfully failed but incorrect error type `{:?}` was returned.", e),
        }
    }
}
