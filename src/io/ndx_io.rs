// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Implementation of functions for reading and writing ndx files.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use indexmap::IndexMap;
use std::collections::HashSet;

use crate::errors::{GroupError, ParseNdxError, WriteNdxError};
use crate::structures::atom::Atom;
use crate::system::System;

/// ## Methods for reading and writing ndx files.
impl System {
    /// Read an ndx file and create atom Groups in the System structure.
    ///
    /// ## Returns
    /// - `Ok` if the parsing is successful.
    /// - `ParseNdxError::InvalidNamesWarning` if any of the groups has an invalid name.
    ///    Has priority over `DuplicateGroupsWarning`.
    /// - `ParseNdxError::DuplicateGroupsWarning` if any of the groups already exists in the system.
    /// - Other `ParseNdxError` errors if the file does not exist or parsing failed.
    ///
    /// ## Notes
    /// - Overwrites all groups with the same names in the system, returning a warning.
    /// - In case duplicate groups are present in the ndx file, the last one
    ///   is input into the system
    /// - The indices in an ndx file do not correspond to atom numbers
    ///   from a gro file, but to actual atom numbers as used by gromacs.
    /// - In case an error other than `ParseNdxError::DuplicateGroupsWarning` or
    ///   `ParseNdxError::InvalidNamesWarning` occurs, the system is not changed.
    /// - Atom numbers can be in any order and will be properly reordered.
    /// - Duplicate atom numbers are ignored.
    /// - Empty lines are skipped.
    pub fn read_ndx(&mut self, filename: impl AsRef<Path>) -> Result<(), ParseNdxError> {
        let file = match File::open(filename.as_ref()) {
            Ok(x) => x,
            Err(_) => return Err(ParseNdxError::FileNotFound(Box::from(filename.as_ref()))),
        };

        let buffer = BufReader::new(file);
        let mut groups: IndexMap<String, Vec<usize>> = IndexMap::new();

        let mut current_name = "".to_string();
        let mut atom_indices = Vec::new();

        let mut duplicate_names: HashSet<String> = HashSet::new();

        for raw_line in buffer.lines() {
            let line = match raw_line {
                Ok(x) => x,
                Err(_) => return Err(ParseNdxError::LineNotFound(Box::from(filename.as_ref()))),
            };

            // skip empty lines
            if line.trim().is_empty() {
                continue;
            }

            // read ndx group name
            if line.contains('[') && line.contains(']') {
                // store previously loaded group
                if !current_name.is_empty()
                    && groups.insert(current_name.clone(), atom_indices).is_some()
                {
                    duplicate_names.insert(current_name);
                }
                atom_indices = Vec::new();

                // read next group name
                current_name = parse_group_name(&line)?;

            // read standard line
            } else {
                atom_indices.extend(parse_ndx_line(&line, self.get_atoms())?);
            }
        }

        // load the last group
        if !current_name.is_empty() && groups.insert(current_name.clone(), atom_indices).is_some() {
            duplicate_names.insert(current_name);
        }

        let mut invalid_names: HashSet<String> = HashSet::new();

        // create groups
        for (groupname, atoms) in groups.into_iter() {
            match self.group_create_from_indices(&groupname, atoms) {
                Ok(_) => (),
                Err(GroupError::AlreadyExistsWarning(e)) => {
                    duplicate_names.insert(e.to_string());
                }
                Err(GroupError::InvalidName(e)) => {
                    invalid_names.insert(e.to_string());
                }
                Err(_) => panic!(
                    "FATAL GROAN ERROR | System::read_ndx | Unexpected error returned from `group_create_from_indices`."
                ),
            }
        }

        if !invalid_names.is_empty() {
            Err(ParseNdxError::InvalidNamesWarning(Box::new(invalid_names)))
        } else if !duplicate_names.is_empty() {
            Err(ParseNdxError::DuplicateGroupsWarning(Box::new(
                duplicate_names,
            )))
        } else {
            Ok(())
        }
    }

    /// Open and write an ndx file using Groups from System as ndx groups.
    ///
    /// ## Returns
    /// `Ok` if writing is successful, else `WriteNdxError`.
    ///
    /// ## Example
    /// Creating groups for residue names and writing them into an ndx file.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// let (_, _residues) = system.atoms_split_by_resname();
    /// if let Err(e) = system.write_ndx("output.ndx") {
    ///     eprintln!("{}", e);
    ///     return;
    /// }
    /// ```
    ///
    /// ## Notes
    /// - Overwrites the contents of any previously existing file with the same `filename`.
    /// - Default System groups such as `all` and `All` are not written out unless a new
    ///   group with such name has been created.
    /// - Groups are written out in the same order as in which they were added into the System,
    ///   unless further manipulated.
    pub fn write_ndx(&self, filename: impl AsRef<Path>) -> Result<(), WriteNdxError> {
        let output = match File::create(&filename) {
            Ok(x) => x,
            Err(_) => return Err(WriteNdxError::CouldNotCreate(Box::from(filename.as_ref()))),
        };

        let mut writer = BufWriter::new(output);

        for (name, group) in self.get_groups().iter() {
            // skip default groups
            if group.print_ndx {
                group.write_ndx(&mut writer, name)?
            };
        }

        writer.flush().map_err(|_| WriteNdxError::CouldNotWrite)?;

        Ok(())
    }
}

/// Parse a line of an ndx file as a group name.
fn parse_group_name(line: &str) -> Result<String, ParseNdxError> {
    let name = line.replace(['[', ']'], "").trim().to_string();

    if name.is_empty() {
        Err(ParseNdxError::ParseGroupNameErr(line.to_string()))
    } else {
        Ok(name)
    }
}

/// Parse a line of an ndx file as gmx atom numbers for an atom Group.
fn parse_ndx_line(line: &str, atoms: &[Atom]) -> Result<Vec<usize>, ParseNdxError> {
    let mut indices = Vec::new();

    for raw_id in line.split_whitespace() {
        let id = match raw_id.parse::<usize>() {
            Ok(x) => x,
            Err(_) => return Err(ParseNdxError::ParseLineErr(line.to_string())),
        };

        if id == 0 || id > atoms.len() {
            return Err(ParseNdxError::InvalidAtomIndex(id));
        }

        indices.push(id - 1);
    }

    Ok(indices)
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests_read_ndx {
    use super::*;

    #[test]
    fn read() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        assert_eq!(system.get_n_groups(), 23);

        // assert that the groups were created
        assert!(system.group_exists("System"));
        assert!(system.group_exists("Protein"));
        assert!(system.group_exists("Protein-H"));
        assert!(system.group_exists("C-alpha"));
        assert!(system.group_exists("Backbone"));
        assert!(system.group_exists("MainChain"));
        assert!(system.group_exists("MainChain+Cb"));
        assert!(system.group_exists("MainChain+H"));
        assert!(system.group_exists("SideChain"));
        assert!(system.group_exists("SideChain-H"));
        assert!(system.group_exists("Prot-Masses"));
        assert!(system.group_exists("non-Protein"));
        assert!(system.group_exists("Other"));
        assert!(system.group_exists("POPC"));
        assert!(system.group_exists("W"));
        assert!(system.group_exists("ION"));
        assert!(system.group_exists("Transmembrane_all"));
        assert!(system.group_exists("Transmembrane"));
        assert!(system.group_exists("Membrane"));
        assert!(system.group_exists("Protein_Membrane"));
        assert!(system.group_exists("W_ION"));

        // assert that the groups have the correct number of atoms
        assert_eq!(system.group_get_n_atoms("System").unwrap(), 16844);
        assert_eq!(system.group_get_n_atoms("Protein").unwrap(), 61);
        assert_eq!(system.group_get_n_atoms("Protein-H").unwrap(), 61);
        assert_eq!(system.group_get_n_atoms("C-alpha").unwrap(), 0);
        assert_eq!(system.group_get_n_atoms("Backbone").unwrap(), 0);
        assert_eq!(system.group_get_n_atoms("MainChain").unwrap(), 0);
        assert_eq!(system.group_get_n_atoms("MainChain+Cb").unwrap(), 0);
        assert_eq!(system.group_get_n_atoms("MainChain+H").unwrap(), 0);
        assert_eq!(system.group_get_n_atoms("SideChain").unwrap(), 61);
        assert_eq!(system.group_get_n_atoms("SideChain-H").unwrap(), 61);
        assert_eq!(system.group_get_n_atoms("Prot-Masses").unwrap(), 61);
        assert_eq!(system.group_get_n_atoms("non-Protein").unwrap(), 16783);
        assert_eq!(system.group_get_n_atoms("Other").unwrap(), 16783);
        assert_eq!(system.group_get_n_atoms("POPC").unwrap(), 6144);
        assert_eq!(system.group_get_n_atoms("W").unwrap(), 10399);
        assert_eq!(system.group_get_n_atoms("ION").unwrap(), 240);
        assert_eq!(system.group_get_n_atoms("Transmembrane_all").unwrap(), 61);
        assert_eq!(system.group_get_n_atoms("Transmembrane").unwrap(), 29);
        assert_eq!(system.group_get_n_atoms("Membrane").unwrap(), 6144);
        assert_eq!(system.group_get_n_atoms("Protein_Membrane").unwrap(), 6205);
        assert_eq!(system.group_get_n_atoms("W_ION").unwrap(), 10639);

        // assert that the groups contain the correct atoms
        for (group_atom, system_atom) in system
            .group_iter("System")
            .unwrap()
            .zip(system.atoms_iter())
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }

        for (group_atom, system_atom) in system
            .group_iter("Protein")
            .unwrap()
            .zip(system.atoms_iter().take(61))
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }

        for (group_atom, system_atom) in system
            .group_iter("Transmembrane_all")
            .unwrap()
            .zip(system.atoms_iter().take(61))
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }

        for (group_atom, system_atom) in system
            .group_iter("W_ION")
            .unwrap()
            .zip(system.atoms_iter().skip(6205))
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }

        for (group_atom, system_atom) in system
            .group_iter("Membrane")
            .unwrap()
            .zip(system.atoms_iter().skip(61).take(6144))
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }
    }

    #[test]
    fn read_small() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        system.read_ndx("test_files/index_small.ndx").unwrap();

        assert_eq!(system.get_n_groups(), 4);

        assert!(system.group_exists("System"));
        assert!(system.group_exists("Protein"));

        assert_eq!(system.group_get_n_atoms("System").unwrap(), 50);
        assert_eq!(system.group_get_n_atoms("Protein").unwrap(), 50);

        // assert that the groups contain the correct atoms
        for (group_atom, system_atom) in system
            .group_iter("System")
            .unwrap()
            .zip(system.atoms_iter())
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }

        for (group_atom, system_atom) in system
            .group_iter("Protein")
            .unwrap()
            .zip(system.atoms_iter())
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }
    }

    #[test]
    fn read_shuffled() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        system.read_ndx("test_files/index_shuffled.ndx").unwrap();

        assert_eq!(system.get_n_groups(), 4);

        assert!(system.group_exists("System"));
        assert!(system.group_exists("Protein"));

        assert_eq!(system.group_get_n_atoms("System").unwrap(), 50);
        assert_eq!(system.group_get_n_atoms("Protein").unwrap(), 50);

        // assert that the groups contain the correct atoms
        for (group_atom, system_atom) in system
            .group_iter("System")
            .unwrap()
            .zip(system.atoms_iter())
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }

        for (group_atom, system_atom) in system
            .group_iter("Protein")
            .unwrap()
            .zip(system.atoms_iter())
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }
    }

    #[test]
    fn red_duplicate_atoms() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        system.read_ndx("test_files/index_duplicate.ndx").unwrap();

        assert_eq!(system.get_n_groups(), 4);

        assert!(system.group_exists("System"));
        assert!(system.group_exists("Protein"));

        assert_eq!(system.group_get_n_atoms("System").unwrap(), 50);
        assert_eq!(system.group_get_n_atoms("Protein").unwrap(), 50);

        // assert that the groups contain the correct atoms
        for (group_atom, system_atom) in system
            .group_iter("System")
            .unwrap()
            .zip(system.atoms_iter())
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }

        for (group_atom, system_atom) in system
            .group_iter("Protein")
            .unwrap()
            .zip(system.atoms_iter())
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }
    }

    #[test]
    fn read_empty() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        system.read_ndx("test_files/index_empty.ndx").unwrap();

        assert_eq!(system.get_n_groups(), 2);

        assert!(!system.group_exists("System"));
        assert!(!system.group_exists("Protein"));
        assert!(system.group_exists("all"));
        assert!(system.group_exists("All"));
    }

    #[test]
    fn read_empy_lines() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        system.read_ndx("test_files/index_empty_lines.ndx").unwrap();

        assert_eq!(system.get_n_groups(), 4);

        assert!(system.group_exists("System"));
        assert!(system.group_exists("Protein"));

        assert_eq!(system.group_get_n_atoms("System").unwrap(), 50);
        assert_eq!(system.group_get_n_atoms("Protein").unwrap(), 50);

        // assert that the groups contain the correct atoms
        for (group_atom, system_atom) in system
            .group_iter("System")
            .unwrap()
            .zip(system.atoms_iter())
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }

        for (group_atom, system_atom) in system
            .group_iter("Protein")
            .unwrap()
            .zip(system.atoms_iter())
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }
    }

    #[test]
    fn read_multiword_group() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        system
            .read_ndx("test_files/index_multiword_group.ndx")
            .unwrap();

        assert_eq!(system.get_n_groups(), 4);

        assert!(system.group_exists("System"));
        assert!(system.group_exists("Protein Named Buforin II P11L"));

        assert_eq!(system.group_get_n_atoms("System").unwrap(), 50);
        assert_eq!(
            system
                .group_get_n_atoms("Protein Named Buforin II P11L")
                .unwrap(),
            50
        );
    }

    macro_rules! read_ndx_fails {
        ($name:ident, $file:expr, $variant:path, $expected:expr) => {
            #[test]
            fn $name() {
                let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
                match system.read_ndx($file) {
                    Err($variant(e)) => assert_eq!(e, $expected),
                    Ok(_) => panic!("Parsing should have failed, but it succeeded."),
                    Err(e) => panic!("Parsing successfully failed but incorrect error type `{:?}` was returned.", e),
                }

                assert!(!system.group_exists("System"));
                assert!(!system.group_exists("Protein"));
                assert!(system.group_exists("all"));
                assert!(system.group_exists("All"));
            }
        };
    }

    read_ndx_fails!(
        read_nonexistent,
        "nonexistent.ndx",
        ParseNdxError::FileNotFound,
        Box::from(Path::new("nonexistent.ndx"))
    );

    read_ndx_fails!(
        read_name_invalid,
        "test_files/index_invalid_name.ndx",
        ParseNdxError::ParseGroupNameErr,
        "[   ] "
    );

    read_ndx_fails!(
        read_unfinished_name,
        "test_files/index_unfinished_name.ndx",
        ParseNdxError::ParseLineErr,
        "[ Protein "
    );

    read_ndx_fails!(
        read_invalid_line,
        "test_files/index_invalid_line.ndx",
        ParseNdxError::ParseLineErr,
        "  16   17   18   19   20   21   -22   23   24   25   26   27   28   29   30"
    );

    read_ndx_fails!(
        read_invalid_index,
        "test_files/index_invalid_index1.ndx",
        ParseNdxError::InvalidAtomIndex,
        0
    );

    read_ndx_fails!(
        read_invalid_index2,
        "test_files/index_invalid_index2.ndx",
        ParseNdxError::InvalidAtomIndex,
        51
    );

    #[test]
    fn read_duplicate_groups() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        match system.read_ndx("test_files/index_duplicate_groups.ndx") {
            Err(ParseNdxError::DuplicateGroupsWarning(e)) => {
                assert_eq!(e, Box::new(HashSet::from(["Protein".to_string()])))
            }
            Ok(_) => panic!("Warning should have been returned, but it was not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }

        assert_eq!(system.get_n_groups(), 4);

        assert!(system.group_exists("System"));
        assert!(system.group_exists("Protein"));

        assert_eq!(system.group_get_n_atoms("System").unwrap(), 50);
        assert_eq!(system.group_get_n_atoms("Protein").unwrap(), 32);
    }

    #[test]
    fn read_duplicate_groups2() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        match system.read_ndx("test_files/index_duplicate_groups2.ndx") {
            Err(ParseNdxError::DuplicateGroupsWarning(e)) => {
                assert_eq!(e, Box::new(HashSet::from(["Protein".to_string()])))
            }
            Ok(_) => panic!("Warning should have been returned, but it was not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }

        assert_eq!(system.get_n_groups(), 4);

        assert!(system.group_exists("System"));
        assert!(system.group_exists("Protein"));

        assert_eq!(system.group_get_n_atoms("System").unwrap(), 50);
        assert_eq!(system.group_get_n_atoms("Protein").unwrap(), 15);
    }

    #[test]
    fn read_group_exists() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        match system.read_ndx("test_files/index_group_exists.ndx") {
            Err(ParseNdxError::DuplicateGroupsWarning(e)) => {
                assert_eq!(e, Box::new(HashSet::from(["All".to_string()])))
            }
            Ok(_) => panic!("Warning should have been returned, but it was not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }

        assert_eq!(system.get_n_groups(), 4);

        assert!(system.group_exists("System"));
        assert!(system.group_exists("Protein"));
        assert!(system.group_exists("All"));

        assert_eq!(system.group_get_n_atoms("System").unwrap(), 50);
        assert_eq!(system.group_get_n_atoms("Protein").unwrap(), 50);
        assert_eq!(system.group_get_n_atoms("All").unwrap(), 35);
    }

    #[test]
    fn read_groups_exist() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        match system.read_ndx("test_files/index_groups_exist.ndx") {
            Err(ParseNdxError::DuplicateGroupsWarning(e)) => assert_eq!(
                e,
                Box::new(HashSet::from(["All".to_string(), "Protein".to_string()]))
            ),
            Ok(_) => panic!("Warning should have been returned, but it was not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }

        assert_eq!(system.get_n_groups(), 4);

        assert!(system.group_exists("System"));
        assert!(system.group_exists("Protein"));

        assert_eq!(system.group_get_n_atoms("System").unwrap(), 50);
        assert_eq!(system.group_get_n_atoms("Protein").unwrap(), 15);
        assert_eq!(system.group_get_n_atoms("All").unwrap(), 35);
    }

    #[test]
    fn read_invalid_names() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        match system.read_ndx("test_files/index_invalid_names.ndx") {
            Err(ParseNdxError::InvalidNamesWarning(e)) => assert_eq!(
                e,
                Box::new(HashSet::from([
                    "inval@id".to_string(),
                    "&also_invalid".to_string(),
                    "(parentheses are invalid)".to_string()
                ]))
            ),
            Ok(_) => panic!("Warning should have been returned, but it was not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }

        assert_eq!(system.get_n_groups(), 4);

        assert!(system.group_exists("System"));
        assert!(system.group_exists("Valid Name"));

        assert_eq!(system.group_get_n_atoms("System").unwrap(), 50);
        assert_eq!(system.group_get_n_atoms("Valid Name").unwrap(), 50);
        assert_eq!(system.group_get_n_atoms("All").unwrap(), 50);
    }
}

#[cfg(test)]
mod tests_write_ndx {
    use super::*;
    use file_diff;
    use tempfile::NamedTempFile;

    #[test]
    fn write() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let ndx_output = NamedTempFile::new().unwrap();
        let path_to_output = ndx_output.path();

        if system.write_ndx(path_to_output).is_err() {
            panic!("Writing ndx file failed.");
        }

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/index.ndx").unwrap();

        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn write_fails() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        match system.write_ndx("Xhfguiedhqueiowhd/nonexistent.ndx") {
            Err(WriteNdxError::CouldNotCreate(e)) => {
                assert_eq!(e, Box::from(Path::new("Xhfguiedhqueiowhd/nonexistent.ndx")))
            }
            Ok(_) => panic!("Writing should have failed, but it did not."),
            Err(e) => panic!("Incorrect error type `{:?}` was returned.", e),
        }
    }
}
