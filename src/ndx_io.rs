// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of functions for reading and writing ndx files.

use std::path::Path;
use std::fs::File;
use std::io::{BufReader, BufRead};

use crate::errors::{ParseNdxError, GroupError};
use crate::group::Group;
use crate::system::System;
use crate::atom::Atom;

impl System {

    /// Read an ndx file and create atom Groups in the System structure.
    /// ## Returns
    /// Ok if the parsing is successful or ParseNdxError if parsing fails 
    /// or any of the groups in an ndx file already exists for the system.
    /// ## Notes
    /// - In case an error occurs, the system is not changed.
    /// - Atom numbers can be in any order and will be properly reordered.
    /// - Duplicate atom numbers are ignored.
    /// - Empty lines are skipped.
    pub fn read_ndx(&mut self, filename: impl AsRef<Path>) -> Result<(), ParseNdxError> {

        let file = match File::open(filename.as_ref()) {
            Ok(x) => x,
            Err(_) => return Err(ParseNdxError::FileNotFound(Box::from(filename.as_ref()))),
        };

        let buffer = BufReader::new(file);
        let mut names: Vec<String> = Vec::new();
        let mut groups: Vec<Group> = Vec::new();

        let mut current_name = "".to_string();
        let mut atom_indices = Vec::new();

        for raw_line in buffer.lines() {
            
            let line = match raw_line {
                Ok(x) => x,
                Err(_) => return Err(ParseNdxError::LineNotFound(Box::from(filename.as_ref()))),
            };
            
            // skip empty lines
            if line.trim().is_empty() { continue; }

            // read ndx group name
            if line.contains("[") && line.contains("]") {
                // create previously loaded group
                add_group(current_name, atom_indices, &mut names, &mut groups, self.get_n_atoms());
                atom_indices = Vec::new();
                
                // read next group name
                current_name = parse_group_name(&line)?;
                // check that the group does not already exist
                if self.group_exists(&current_name) {
                    return Err(ParseNdxError::GroupAlreadyExists(current_name));
                }

                // check that there are no duplicate group names in the ndx file
                if let Some(_) = names.iter().position(|x| *x == current_name) {
                    return Err(ParseNdxError::GroupsShareName(current_name));
                }
                

            // read standard line
            } else {
                atom_indices.extend(parse_ndx_line(&line, self.get_atoms_as_ref())?);
            }
        }

        add_group(current_name, atom_indices, &mut names, &mut groups, self.get_n_atoms());
    
        // transfer all groups into the system
        for (name, group) in names.into_iter().zip(groups.into_iter()) {
            if let Err(_) = self.group_add(&name, group) {
                panic!("Internal Error. Group {} already exists, which should be impossible.", name);
            }
        }

        Ok(())
    }


    /// Adds an already constructed group into the system. 
    /// ## Returns
    /// Ok if the group was added or GroupError in case a group with the same name already exists.
    fn group_add(&mut self, name: &str, group: Group) -> Result<(), GroupError> {
        unsafe {
            match self.get_groups_as_ref_mut().insert(name.to_string(), group) {
                Some(_) => Err(GroupError::AlreadyExists(name.to_string())),
                None => Ok(())
            }
        }
    }
}

/// Parse a line of an ndx file as a group name.
fn parse_group_name(line: &str) -> Result<String, ParseNdxError> {
    let name = line.replace("[", "").replace("]", "").trim().to_string();
    
    if name.is_empty() {
        Err(ParseNdxError::ParseGroupNameErr(line.to_string()))
    } else {
        Ok(name)
    }
}

/// Parse a line of an ndx file as gmx atom numbers for an atom Group.
/// ## Panics
/// - If the atoms in the `atoms` variable are not sorted by their gmx atom numbers. 
/// In `atoms`, the gmx atom numbers must be continuously increasing.
fn parse_ndx_line(line: &str, atoms: &Vec<Atom>) -> Result<Vec<u64>, ParseNdxError> {
    let mut indices = Vec::new();
    
    for raw_id in line.split_whitespace() {
        let id = match raw_id.parse::<u64>() {
            Ok(x) => x,
            Err(_) => return Err(ParseNdxError::ParseLineErr(line.to_string())),
        };

        if id == 0 {
            return Err(ParseNdxError::InvalidAtomIndex(id));
        }

        // check that the index matches the expected gmx atom number
        match atoms.get(id as usize - 1) {
            Some(atom) => {
                if atom.get_gmx_atom_number() == id {
                    indices.push(id - 1);
                } else {
                    panic!("Internal Error. Expected gmx atom number {} was not found for atom at index {}.", id, id - 1);
                }
            },
            // index is out of range, return Error
            None => return Err(ParseNdxError::InvalidAtomIndex(id)),
        }
    }

    Ok(indices)
}

/// Create a group and add it and its name to vectors `groups` and `names`, respectively.
fn add_group(
    name: String, 
    atom_indices: Vec<u64>, 
    names: &mut Vec<String>, 
    groups: &mut Vec<Group>,
    n_atoms: u64) {
    
    if !name.is_empty() {
        names.push(name);
        groups.push(Group::from_indices(atom_indices, n_atoms));
    }
}




#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_ndx() { 

        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

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
        assert_eq!(system.group_n_atoms("System").unwrap(), 16844);
        assert_eq!(system.group_n_atoms("Protein").unwrap(), 61);
        assert_eq!(system.group_n_atoms("Protein-H").unwrap(), 61);
        assert_eq!(system.group_n_atoms("C-alpha").unwrap(), 0);
        assert_eq!(system.group_n_atoms("Backbone").unwrap(), 0);
        assert_eq!(system.group_n_atoms("MainChain").unwrap(), 0);
        assert_eq!(system.group_n_atoms("MainChain+Cb").unwrap(), 0);
        assert_eq!(system.group_n_atoms("MainChain+H").unwrap(), 0);
        assert_eq!(system.group_n_atoms("SideChain").unwrap(), 61);
        assert_eq!(system.group_n_atoms("SideChain-H").unwrap(), 61);
        assert_eq!(system.group_n_atoms("Prot-Masses").unwrap(), 61);
        assert_eq!(system.group_n_atoms("non-Protein").unwrap(), 16783);
        assert_eq!(system.group_n_atoms("Other").unwrap(), 16783);
        assert_eq!(system.group_n_atoms("POPC").unwrap(), 6144);
        assert_eq!(system.group_n_atoms("W").unwrap(), 10399);
        assert_eq!(system.group_n_atoms("ION").unwrap(), 240);
        assert_eq!(system.group_n_atoms("Transmembrane_all").unwrap(), 61);
        assert_eq!(system.group_n_atoms("Transmembrane").unwrap(), 29);
        assert_eq!(system.group_n_atoms("Membrane").unwrap(), 6144);
        assert_eq!(system.group_n_atoms("Protein_Membrane").unwrap(), 6205);
        assert_eq!(system.group_n_atoms("W_ION").unwrap(), 10639);

        // assert that the groups contain the correct atoms
        for (group_atom, system_atom) in system.group_iter("System").unwrap().zip(system.atoms_iter()) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
        }

        for (group_atom, system_atom) in system.group_iter("Protein").unwrap().zip(system.atoms_iter().take(61)) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
        }

        for (group_atom, system_atom) in system.group_iter("Transmembrane_all").unwrap().zip(system.atoms_iter().take(61)) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
        }

        for (group_atom, system_atom) in system.group_iter("W_ION").unwrap().zip(system.atoms_iter().skip(6205)) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
        }

        for (group_atom, system_atom) in system.group_iter("Membrane").unwrap().zip(system.atoms_iter().skip(61).take(6144)) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
        }
    }

    #[test]
    fn test_read_ndx_small() { 

        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        system.read_ndx("test_files/index_small.ndx").unwrap();

        assert!(system.group_exists("System"));
        assert!(system.group_exists("Protein"));

        assert_eq!(system.group_n_atoms("System").unwrap(), 50);
        assert_eq!(system.group_n_atoms("Protein").unwrap(), 50);

        // assert that the groups contain the correct atoms
        for (group_atom, system_atom) in system.group_iter("System").unwrap().zip(system.atoms_iter()) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
        }

        for (group_atom, system_atom) in system.group_iter("Protein").unwrap().zip(system.atoms_iter()) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
        }
    }

    #[test]
    fn test_read_ndx_shuffled() { 

        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        system.read_ndx("test_files/index_shuffled.ndx").unwrap();

        assert!(system.group_exists("System"));
        assert!(system.group_exists("Protein"));

        assert_eq!(system.group_n_atoms("System").unwrap(), 50);
        assert_eq!(system.group_n_atoms("Protein").unwrap(), 50);

        // assert that the groups contain the correct atoms
        for (group_atom, system_atom) in system.group_iter("System").unwrap().zip(system.atoms_iter()) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
        }

        for (group_atom, system_atom) in system.group_iter("Protein").unwrap().zip(system.atoms_iter()) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
        }
    }

    #[test]
    fn test_read_ndx_duplicate_atoms() { 

        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        system.read_ndx("test_files/index_duplicate.ndx").unwrap();

        assert!(system.group_exists("System"));
        assert!(system.group_exists("Protein"));

        assert_eq!(system.group_n_atoms("System").unwrap(), 50);
        assert_eq!(system.group_n_atoms("Protein").unwrap(), 50);

        // assert that the groups contain the correct atoms
        for (group_atom, system_atom) in system.group_iter("System").unwrap().zip(system.atoms_iter()) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
        }

        for (group_atom, system_atom) in system.group_iter("Protein").unwrap().zip(system.atoms_iter()) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
        }
    }

    #[test]
    fn test_read_ndx_empty() { 

        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        system.read_ndx("test_files/index_empty.ndx").unwrap();

        assert!(!system.group_exists("System"));
        assert!(!system.group_exists("Protein"));
        assert!(system.group_exists("all"));
        assert!(system.group_exists("All"));
    }

    #[test]
    fn test_read_ndx_empy_lines() { 

        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();
        system.read_ndx("test_files/index_empty_lines.ndx").unwrap();

        assert!(system.group_exists("System"));
        assert!(system.group_exists("Protein"));

        assert_eq!(system.group_n_atoms("System").unwrap(), 50);
        assert_eq!(system.group_n_atoms("Protein").unwrap(), 50);

        // assert that the groups contain the correct atoms
        for (group_atom, system_atom) in system.group_iter("System").unwrap().zip(system.atoms_iter()) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
        }

        for (group_atom, system_atom) in system.group_iter("Protein").unwrap().zip(system.atoms_iter()) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
        }
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

    read_ndx_fails!(test_read_ndx_nonexistent, "nonexistent.ndx", 
    ParseNdxError::FileNotFound, Box::from(Path::new("nonexistent.ndx")));

    read_ndx_fails!(test_read_ndx_duplicate_groups, "test_files/index_duplicate_groups.ndx", 
    ParseNdxError::GroupsShareName, "Protein");

    read_ndx_fails!(test_read_ndx_group_exists, "test_files/index_group_exists.ndx", 
    ParseNdxError::GroupAlreadyExists, "All");

    read_ndx_fails!(test_read_ndx_name_invalid, "test_files/index_invalid_name.ndx", 
    ParseNdxError::ParseGroupNameErr, "[   ] ");

    read_ndx_fails!(test_read_ndx_unfinished_name, "test_files/index_unfinished_name.ndx", 
    ParseNdxError::ParseLineErr, "[ Protein ");

    read_ndx_fails!(test_read_ndx_invalid_line, "test_files/index_invalid_line.ndx", 
    ParseNdxError::ParseLineErr, "  16   17   18   19   20   21   -22   23   24   25   26   27   28   29   30");

    read_ndx_fails!(test_read_ndx_invalid_index, "test_files/index_invalid_index1.ndx", 
    ParseNdxError::InvalidAtomIndex, 0);

    read_ndx_fails!(test_read_ndx_invalid_index2, "test_files/index_invalid_index2.ndx", 
    ParseNdxError::InvalidAtomIndex, 51);

}