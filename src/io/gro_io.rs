// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of functions for reading and writing gro files.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::auxiliary::{GRO_MAX_COORDINATE, GRO_MIN_COORDINATE};
use crate::errors::{ParseGroError, WriteGroError};
use crate::structures::{atom::Atom, simbox::SimBox};
use crate::system::System;

use super::check_coordinate_sizes;

/// ## Methods for writing gro files.
impl System {
    /// Write all atoms of the `System` into a gro file with the given name.
    ///
    /// ## Returns
    /// `Ok` if writing has been successful. Otherwise `WriteGroError`.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let system = System::from_file("system.gro").unwrap();
    /// if let Err(e) = system.write_gro("system_copy.gro", true) {
    ///     eprintln!("{}", e);
    ///     return;
    /// }
    /// ```
    ///
    /// ## Notes
    /// - The function will write velocities for atoms only if `write_velocities == true`.
    /// - The function will write all 9 box coordinates only if necessary
    /// (any of the last 6 coordinates is non-zero). Otherwise, it assumes the box is
    /// orthogonal and writes out only 3 dimensions of the box.
    /// - If simulation box is undefined, it is written as a sequence of zeros.
    pub fn write_gro(
        &self,
        filename: impl AsRef<Path>,
        write_velocities: bool,
    ) -> Result<(), WriteGroError> {
        match self.group_write_gro("all", filename, write_velocities) {
            Ok(_) => Ok(()),
            Err(WriteGroError::GroupNotFound(_)) => {
                panic!(
                    "FATAL GROAN ERROR | System::write_gro | Default group 'all' does not exist."
                )
            }
            Err(e) => Err(e),
        }
    }

    /// Write atoms of the specified group into a gro file with the given name.
    ///
    /// ## Returns
    /// `Ok` if writing has been successful. Otherwise `WriteGroError`.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// if let Err(e) = system.group_write_gro("Protein", "protein.gro", true) {
    ///     eprintln!("{}", e);
    ///     return;
    /// }
    /// ```
    ///
    /// ## Notes
    /// - The function will write velocities for atoms only if `write_velocities == true`.
    /// - The function will write all 9 box coordinates only if necessary
    /// (any of the last 6 coordinates is non-zero). Otherwise, it assumes the box is
    /// orthogonal and writes out only 3 dimensions of the box.
    /// - If the simulation box is undefined, it is written as a sequence of zeros.
    pub fn group_write_gro(
        &self,
        group_name: &str,
        filename: impl AsRef<Path>,
        write_velocities: bool,
    ) -> Result<(), WriteGroError> {
        if !self.group_exists(group_name) {
            return Err(WriteGroError::GroupNotFound(group_name.to_string()));
        }

        // check that coordinates of the atoms are in the range supported by the data format
        // this has to be done before the file is even created
        if !check_coordinate_sizes(
            self.group_iter(group_name).unwrap(),
            GRO_MIN_COORDINATE,
            GRO_MAX_COORDINATE,
        ) {
            return Err(WriteGroError::CoordinateTooLarge);
        }

        let output = File::create(&filename)
            .map_err(|_| WriteGroError::CouldNotCreate(Box::from(filename.as_ref())))?;

        let mut writer = BufWriter::new(output);

        // write gro file header
        let title = match group_name {
            "all" => self.get_name().to_owned(),
            _ => format!("Group `{}` from {}", group_name, self.get_name()),
        };

        write_header(
            &mut writer,
            &title,
            self.group_get_n_atoms(group_name).unwrap(),
        )?;

        // write atoms
        for atom in self.group_iter(group_name).unwrap() {
            atom.write_gro(&mut writer, write_velocities)?;
        }

        // write simulation box
        self.write_box(&mut writer)?;

        writer.flush().map_err(|_| WriteGroError::CouldNotWrite)?;

        Ok(())
    }

    /// Write box dimensions into an open gro file.
    fn write_box(&self, writer: &mut BufWriter<File>) -> Result<(), WriteGroError> {
        match self.get_box_as_ref() {
            Some(simbox) if simbox.is_orthogonal() => {
                writeln!(
                    writer,
                    " {:9.5} {:9.5} {:9.5}",
                    simbox.x, simbox.y, simbox.z
                )
                .map_err(|_| WriteGroError::CouldNotWrite)?;
            }
            Some(simbox) => {
                writeln!(
                    writer,
                    " {:9.5} {:9.5} {:9.5} {:9.5} {:9.5} {:9.5} {:9.5} {:9.5} {:9.5}",
                    simbox.x,
                    simbox.y,
                    simbox.z,
                    simbox.v1y,
                    simbox.v1z,
                    simbox.v2x,
                    simbox.v2z,
                    simbox.v3x,
                    simbox.v3y
                )
                .map_err(|_| WriteGroError::CouldNotWrite)?;
            }
            None => {
                let x = 0.0;
                writeln!(writer, " {x:9.5} {x:9.5} {x:9.5}",)
                    .map_err(|_| WriteGroError::CouldNotWrite)?;
            }
        }
        Ok(())
    }
}

/// Read a gro file and construct a System structure.
pub fn read_gro(filename: impl AsRef<Path>) -> Result<System, ParseGroError> {
    let file = match File::open(filename.as_ref()) {
        Ok(x) => x,
        Err(_) => return Err(ParseGroError::FileNotFound(Box::from(filename.as_ref()))),
    };

    let mut buffer = BufReader::new(file);

    // get title and number of atoms
    let title = get_title(&mut buffer, filename.as_ref())?;
    let n_atoms = get_natoms(&mut buffer, filename.as_ref())?;
    let mut simulation_box = None;

    let mut atoms: Vec<Atom> = Vec::with_capacity(n_atoms);

    // parse all remaining lines
    for (gmx_index, raw_line) in buffer.lines().enumerate() {
        let line = match raw_line {
            Ok(x) => x,
            Err(_) => return Err(ParseGroError::LineNotFound(Box::from(filename.as_ref()))),
        };

        if gmx_index == n_atoms {
            simulation_box = Some(line_as_box(&line)?);
            if simulation_box.as_ref().unwrap().is_zero() {
                simulation_box = None;
            }
        } else {
            let atom = line_as_atom(&line)?;
            atoms.push(atom);
        }
    }

    if atoms.len() != n_atoms {
        return Err(ParseGroError::LineNotFound(Box::from(filename.as_ref())));
    }

    Ok(System::new(&title, atoms, simulation_box))
}

/// Read the next line in the provided buffer and parse it as a title.
fn get_title(
    buffer: &mut BufReader<File>,
    filename: impl AsRef<Path>,
) -> Result<String, ParseGroError> {
    let mut title = String::new();
    match buffer.read_line(&mut title) {
        Ok(0) | Err(_) => return Err(ParseGroError::LineNotFound(Box::from(filename.as_ref()))),
        Ok(_) => return Ok(title.trim().to_string()),
    };
}

/// Read the next line in the provided buffer and parse it as the number of atoms.
fn get_natoms(
    buffer: &mut BufReader<File>,
    filename: impl AsRef<Path>,
) -> Result<usize, ParseGroError> {
    let mut line = String::new();
    match buffer.read_line(&mut line) {
        Ok(0) | Err(_) => Err(ParseGroError::LineNotFound(Box::from(filename.as_ref()))),
        Ok(_) => match line.trim().parse::<usize>() {
            Ok(x) => Ok(x),
            Err(_) => Err(ParseGroError::ParseLineErr(line.trim().to_string())),
        },
    }
}

/// Parse a line as atom.
fn line_as_atom(line: &str) -> Result<Atom, ParseGroError> {
    if line.len() < 44 {
        return Err(ParseGroError::ParseAtomLineErr(line.to_string()));
    }

    // parse residue number
    let resid = line[0..5]
        .trim()
        .parse::<usize>()
        .map_err(|_| ParseGroError::ParseAtomLineErr(line.to_string()))?;

    // parse residue name
    let resname = line[5..10].trim().to_string();
    if resname.is_empty() {
        return Err(ParseGroError::ParseAtomLineErr(line.to_string()));
    }

    // parse atom name
    let atomname = line[10..15].trim().to_string();
    if atomname.is_empty() {
        return Err(ParseGroError::ParseAtomLineErr(line.to_string()));
    }

    // parse atom number
    let atomid = line[15..20]
        .trim()
        .parse::<usize>()
        .map_err(|_| ParseGroError::ParseAtomLineErr(line.to_string()))?;

    // parse position
    let mut position = [0.0; 3];
    for (i, item) in position.iter_mut().enumerate() {
        let curr = 20 + i * 8;
        *item = line[curr..curr + 8]
            .trim()
            .parse::<f32>()
            .map_err(|_| ParseGroError::ParseAtomLineErr(line.to_string()))?;
    }

    let atom = Atom::new(resid, &resname, atomid, &atomname).with_position(position.into());

    // parse velocity, if present
    if line.len() >= 68 {
        let mut velocity = [0.0; 3];

        for (i, item) in velocity.iter_mut().enumerate() {
            let curr = 44 + i * 8;
            *item = line[curr..curr + 8]
                .trim()
                .parse::<f32>()
                .map_err(|_| ParseGroError::ParseAtomLineErr(line.to_string()))?;
        }

        Ok(atom.with_velocity(velocity.into()))
    } else {
        Ok(atom)
    }
}

/// Parse a line as simulation box dimensions.
fn line_as_box(line: &str) -> Result<SimBox, ParseGroError> {
    let mut simulation_box = [0.0f32; 9];
    let mut i = 0usize;
    for split in line.split_whitespace() {
        simulation_box[i] = split
            .trim()
            .parse::<f32>()
            .map_err(|_| ParseGroError::ParseBoxLineErr(line.to_string()))?;
        i += 1;
    }

    if i != 3 && i != 9 {
        Err(ParseGroError::ParseBoxLineErr(line.to_string()))?;
    }

    // check that the simulation box is valid
    if simulation_box[3] != 0.0 || simulation_box[4] != 0.0 || simulation_box[6] != 0.0 {
        return Err(ParseGroError::UnsupportedBox(line.to_string()));
    }

    Ok(simulation_box.into())
}

/// Write gro file header into an open gro file.
fn write_header(
    writer: &mut BufWriter<File>,
    title: &str,
    n_atoms: usize,
) -> Result<(), WriteGroError> {
    writeln!(writer, "{}", title).map_err(|_| WriteGroError::CouldNotWrite)?;

    writeln!(writer, "{:>5}", n_atoms).map_err(|_| WriteGroError::CouldNotWrite)?;

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
    fn read() {
        let system = read_gro("test_files/example.gro").unwrap();

        // check the system's name
        assert_eq!(
            system.get_name(),
            "INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1"
        );

        // check the number of atoms
        assert_eq!(system.get_n_atoms(), 16844);

        // check box size
        let simbox = system.get_box_as_ref().unwrap();
        assert!(approx_eq!(f32, simbox.x, 13.01331));
        assert!(approx_eq!(f32, simbox.y, 13.01331));
        assert!(approx_eq!(f32, simbox.z, 11.25347));

        assert!(approx_eq!(f32, simbox.v1x, 13.01331));
        assert!(approx_eq!(f32, simbox.v2y, 13.01331));
        assert!(approx_eq!(f32, simbox.v3z, 11.25347));

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
        assert_eq!(first.get_residue_name(), "GLY");
        assert_eq!(first.get_atom_name(), "BB");
        assert_eq!(first.get_atom_number(), 1);

        assert!(approx_eq!(f32, first.get_position().unwrap().x, 9.497));
        assert!(approx_eq!(f32, first.get_position().unwrap().y, 1.989));
        assert!(approx_eq!(f32, first.get_position().unwrap().z, 7.498));

        assert!(approx_eq!(f32, first.get_velocity().unwrap().x, -0.0683));
        assert!(approx_eq!(f32, first.get_velocity().unwrap().y, 0.1133));
        assert!(approx_eq!(f32, first.get_velocity().unwrap().z, 0.0005));

        assert_eq!(first.get_force(), None);

        // check atom somewhere in the middle
        let middle = &atoms[4932];
        assert_eq!(middle.get_residue_number(), 435);
        assert_eq!(middle.get_residue_name(), "POPC");
        assert_eq!(middle.get_atom_name(), "C4B");
        assert_eq!(middle.get_atom_number(), 4933);

        assert!(approx_eq!(f32, middle.get_position().unwrap().x, 6.384));
        assert!(approx_eq!(f32, middle.get_position().unwrap().y, 11.908));
        assert!(approx_eq!(f32, middle.get_position().unwrap().z, 5.471));

        assert!(approx_eq!(f32, middle.get_velocity().unwrap().x, -0.2271));
        assert!(approx_eq!(f32, middle.get_velocity().unwrap().y, 0.1287));
        assert!(approx_eq!(f32, middle.get_velocity().unwrap().z, 0.1784));

        assert_eq!(middle.get_force(), None);

        // check the last atom
        let last = &atoms[16843];
        assert_eq!(last.get_residue_number(), 11180);
        assert_eq!(last.get_residue_name(), "ION");
        assert_eq!(last.get_atom_name(), "CL");
        assert_eq!(last.get_atom_number(), 16844);

        assert!(approx_eq!(f32, last.get_position().unwrap().x, 8.829));
        assert!(approx_eq!(f32, last.get_position().unwrap().y, 11.186));
        assert!(approx_eq!(f32, last.get_position().unwrap().z, 2.075));

        assert!(approx_eq!(f32, last.get_velocity().unwrap().x, 0.0712));
        assert!(approx_eq!(f32, last.get_velocity().unwrap().y, 0.2294));
        assert!(approx_eq!(f32, last.get_velocity().unwrap().z, -0.1673));

        assert_eq!(last.get_force(), None);
    }

    #[test]
    fn read_novelocities() {
        let system = read_gro("test_files/example_novelocities.gro").unwrap();

        assert_eq!(system.get_name(), "Buforin II peptide P11L");
        assert_eq!(system.get_n_atoms(), 50);

        // check box size
        let simbox = system.get_box_as_ref().unwrap();
        assert!(approx_eq!(f32, simbox.x, 6.08608));
        assert!(approx_eq!(f32, simbox.y, 6.08608));
        assert!(approx_eq!(f32, simbox.z, 6.08608));

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

        assert!(approx_eq!(f32, first.get_position().unwrap().x, 1.660));
        assert!(approx_eq!(f32, first.get_position().unwrap().y, 2.061));
        assert!(approx_eq!(f32, first.get_position().unwrap().z, 3.153));

        // check atom somewhere in the middle
        let middle = &atoms[24];
        assert_eq!(middle.get_residue_number(), 11);
        assert_eq!(middle.get_residue_name(), "LEU");
        assert_eq!(middle.get_atom_name(), "SC1");
        assert_eq!(middle.get_atom_number(), 25);

        assert!(approx_eq!(f32, middle.get_position().unwrap().x, 3.161));
        assert!(approx_eq!(f32, middle.get_position().unwrap().y, 2.868));
        assert!(approx_eq!(f32, middle.get_position().unwrap().z, 2.797));

        // check the last atom
        let last = &atoms[49];
        assert_eq!(last.get_residue_number(), 21);
        assert_eq!(last.get_residue_name(), "LYS");
        assert_eq!(last.get_atom_name(), "SC2");
        assert_eq!(last.get_atom_number(), 50);

        assert!(approx_eq!(f32, last.get_position().unwrap().x, 4.706));
        assert!(approx_eq!(f32, last.get_position().unwrap().y, 4.447));
        assert!(approx_eq!(f32, last.get_position().unwrap().z, 2.813));

        // check that the velocity and force of all atoms is zero
        for atom in atoms.iter() {
            assert_eq!(atom.get_velocity(), None);
            assert_eq!(atom.get_force(), None);
        }
    }

    #[test]
    fn read_box9() {
        let system = read_gro("test_files/example_box9.gro").unwrap();

        let simbox = system.get_box_as_ref().unwrap();
        assert!(approx_eq!(f32, simbox.x, 6.08608));
        assert!(approx_eq!(f32, simbox.y, 6.08608));
        assert!(approx_eq!(f32, simbox.z, 6.08608));

        assert!(approx_eq!(f32, simbox.v1y, 0.0));
        assert!(approx_eq!(f32, simbox.v1z, 0.0));
        assert!(approx_eq!(f32, simbox.v2x, 2.2));

        assert!(approx_eq!(f32, simbox.v2z, 0.0));
        assert!(approx_eq!(f32, simbox.v3x, 1.4));
        assert!(approx_eq!(f32, simbox.v3y, 3.856));
    }

    #[test]
    fn read_box_zero() {
        let system = read_gro("test_files/example_box_zero.gro").unwrap();
        assert!(!system.has_box());
    }

    #[test]
    fn read_nonexistent() {
        if read_gro("test_files/nonexistent.gro").is_ok() {
            panic!("Nonexistent file seems to exist.");
        }
    }

    macro_rules! read_gro_fails {
        ($name:ident, $file:expr, $variant:path, $expected:expr) => {
            #[test]
            fn $name() {
                let file = $file;
                match read_gro(file) {
                    Err($variant(e)) => assert_eq!(e, $expected),
                    Ok(_) => panic!("Parsing should have failed, but it succeeded."),
                    Err(e) => panic!("Parsing successfully failed but incorrect error type `{:?}` was returned.", e),
                }
            }
        };
    }

    read_gro_fails!(
        read_incomplete_line,
        "test_files/example_incomplete_line.gro",
        ParseGroError::ParseAtomLineErr,
        "   16HIS    SC1   35   3.458   3.653   "
    );

    read_gro_fails!(
        read_empty_file,
        "test_files/example_empty.gro",
        ParseGroError::LineNotFound,
        Box::from(Path::new("test_files/example_empty.gro"))
    );

    read_gro_fails!(
        read_only_title,
        "test_files/example_only_title.gro",
        ParseGroError::LineNotFound,
        Box::from(Path::new("test_files/example_only_title.gro"))
    );

    read_gro_fails!(
        read_missing_natoms,
        "test_files/example_missing_natoms.gro",
        ParseGroError::ParseLineErr,
        ""
    );

    read_gro_fails!(
        read_unparsable_natoms,
        "test_files/example_unparsable_natoms.gro",
        ParseGroError::ParseLineErr,
        "6A5F"
    );

    read_gro_fails!(
        read_missing_atom,
        "test_files/example_missing_atom.gro",
        ParseGroError::ParseAtomLineErr,
        "   6.08608   6.08608   6.08608"
    );

    read_gro_fails!(
        read_invalid_resid,
        "test_files/example_invalid_resid.gro",
        ParseGroError::ParseAtomLineErr,
        "   1APHE    SC2   22   2.519   3.025   3.387"
    );

    read_gro_fails!(
        read_invalid_atomid,
        "test_files/example_invalid_atomid.gro",
        ParseGroError::ParseAtomLineErr,
        "   21LYS     BB        4.362   4.008   3.161"
    );

    read_gro_fails!(
        read_invalid_position,
        "test_files/example_invalid_position.gro",
        ParseGroError::ParseAtomLineErr,
        "    2ARG    SC1    4   1.877   1. 73   3.023"
    );

    read_gro_fails!(
        read_invalid_velocity,
        "test_files/example_invalid_velocity.gro",
        ParseGroError::ParseAtomLineErr,
        "   15LEU    SC1   31   9.638   2.052   5.595  0.0685  O.0634  0.1453"
    );

    read_gro_fails!(
        read_shifted_line,
        "test_files/example_shifted_line.gro",
        ParseGroError::ParseAtomLineErr,
        "    20ARG     BB   45   4.265   3.832   2.925"
    );

    read_gro_fails!(
        read_empty_box_line,
        "test_files/example_empty_box_line.gro",
        ParseGroError::ParseBoxLineErr,
        ""
    );

    read_gro_fails!(
        read_short_box,
        "test_files/example_short_box.gro",
        ParseGroError::ParseBoxLineErr,
        "   6.08608   6.08608"
    );

    read_gro_fails!(
        read_long_box,
        "test_files/example_long_box.gro",
        ParseGroError::ParseBoxLineErr,
        "   6.08608   6.08608   6.08608   6.08608   6.08608"
    );

    read_gro_fails!(
        read_unparsable_box,
        "test_files/example_unparsable_box.gro",
        ParseGroError::ParseBoxLineErr,
        "   6.08608   6.08608   6,08608"
    );

    read_gro_fails!(
        read_unsupported_box,
        "test_files/example_unsupported_box.gro",
        ParseGroError::UnsupportedBox,
        "   6.08608   6.08608   6.08608   0.00000   0.60000   2.20000   0.00000   1.40000   3.85600"
    );

    #[test]
    fn from_file() {
        let system = System::from_file("test_files/example.gro").expect("File not found.");

        assert_eq!(
            system.get_name(),
            "INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1"
        );
        assert_eq!(system.get_n_atoms(), 16844);

        assert!(system.group_exists("all"));
    }

    #[test]
    fn from_file_fails() {
        if System::from_file("test_files/example_invalid_position.gro").is_ok() {
            panic!("Parsing should have failed, but it succeeded.")
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
        let system = System::from_file("test_files/example.gro").unwrap();

        let gro_output = NamedTempFile::new().unwrap();
        let path_to_output = gro_output.path();

        if system.write_gro(path_to_output, true).is_err() {
            panic!("Writing gro file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example.gro").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_no_velocities() {
        let system = System::from_file("test_files/example_novelocities.gro").unwrap();

        let gro_output = NamedTempFile::new().unwrap();
        let path_to_output = gro_output.path();

        if system.write_gro(path_to_output, false).is_err() {
            panic!("Writing gro file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example_novelocities.gro").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_fails() {
        let system = System::from_file("test_files/example.gro").unwrap();

        match system.write_gro("Xhfguiedhqueiowhd/nonexistent.ndx", true) {
            Err(WriteGroError::CouldNotCreate(e)) => {
                assert_eq!(e, Box::from(Path::new("Xhfguiedhqueiowhd/nonexistent.ndx")))
            }
            Ok(_) => panic!("Writing should have failed, but it did not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }
    }

    #[test]
    fn write_gro9() {
        let system = System::from_file("test_files/example_box9.gro").unwrap();

        let gro_output = NamedTempFile::new().unwrap();
        let path_to_output = gro_output.path();

        if system.write_gro(path_to_output, false).is_err() {
            panic!("Writing gro file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example_box9.gro").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_wrap() {
        let atom1 = Atom::new(158, "THR", 1, "BBBBBT");

        let atom2 = Atom::new(158, "THR", 99999, "SC1");

        let atom3 = Atom::new(100003, "ARG", 100000, "BB");

        let atom4 = Atom::new(100003, "ARGGGT", 200001, "SC1");

        let atom5 = Atom::new(100003, "ARG", 200005, "SC2");

        let atoms = vec![atom1, atom2, atom3, atom4, atom5];
        let simbox = SimBox::from([1.0, 1.0, 1.0]);

        let system = System::new("Expected atom and residue wrapping", atoms, Some(simbox));

        let gro_output = NamedTempFile::new().unwrap();
        let path_to_output = gro_output.path();

        if system.write_gro(path_to_output, false).is_err() {
            panic!("Writing gro file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/wrapping_expected.gro").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.read_ndx("test_files/index.ndx").unwrap();

        let gro_output = NamedTempFile::new().unwrap();
        let path_to_output = gro_output.path();

        if system
            .group_write_gro("Protein", path_to_output, true)
            .is_err()
        {
            panic!("Writing gro file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/protein.gro").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_group_fails() {
        let system = System::from_file("test_files/example.gro").unwrap();

        let gro_output = NamedTempFile::new().unwrap();
        let path_to_output = gro_output.path();

        match system.group_write_gro("Protein", path_to_output, true) {
            Err(WriteGroError::GroupNotFound(e)) => assert_eq!(e, "Protein"),
            Ok(_) => panic!("Writing should have failed, but it did not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }
    }

    #[test]
    fn write_no_box() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        system.reset_box();

        let gro_output = NamedTempFile::new().unwrap();
        let path_to_output = gro_output.path();

        system.write_gro(path_to_output, false).unwrap();

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/example_box_zero.gro").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_too_large_coordinate_1() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.get_atom_as_mut(16).unwrap().set_position_z(10000.0);

        match system.write_gro("will_not_be_created_1.gro", false) {
            Err(WriteGroError::CoordinateTooLarge) => (),
            Ok(_) => panic!("Writing should have failed, but it did not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }
    }

    #[test]
    fn write_too_large_coordinate_2() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.get_atom_as_mut(16).unwrap().set_position_x(-9999.0);

        match system.group_write_gro("all", "will_not_be_created_2.gro", false) {
            Err(WriteGroError::CoordinateTooLarge) => (),
            Ok(_) => panic!("Writing should have failed, but it did not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }
    }
}
