// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of the Group structure and its methods.

use std::io::Write;

use crate::errors::{SelectError, WriteNdxError};
use crate::selections::select::{self, Select};
use crate::structures::shape::Shape;
use crate::system::general::System;

/******************************/
/*       GROUP STRUCTURE      */
/******************************/

/// Group of atoms in target system.
#[derive(Debug, Clone)]
pub struct Group {
    pub atom_ranges: Vec<(usize, usize)>,
    pub print_ndx: bool,
}

impl Group {
    /// Get immutable reference to the atom ranges of the Group.
    pub fn get_atom_ranges(&self) -> &Vec<(usize, usize)> {
        &self.atom_ranges
    }

    /// Create a new valid `Group` structure from query in Groan Selection Language.
    pub fn from_query(query: &str, system: &System) -> Result<Self, SelectError> {
        // parse groan selection language query into binary selection tree
        let select = select::parse_query(query)?;
        // apply the selection tree to the system
        Group::from_select(select, system)
    }

    /// Create a new valid `Group` structure from query in Groan Selection Language and from geometry constraint.
    pub fn from_query_and_geometry(
        query: &str,
        geometry: impl Shape,
        system: &System,
    ) -> Result<Self, SelectError> {
        // create group from query
        let group = Group::from_query(query, system)?;
        // now apply geometry to the group
        Ok(group.apply_geometry(geometry, system))
    }

    /// Create a new valid Group structure from atom ranges.
    ///
    /// ## Parameters
    /// Expects a vector of atom ranges (start, end) and the total number of atoms in the system.
    /// The individual "atom ranges" in the atom ranges can be overlapping.
    /// In the final Group structure, they will be merged together.
    pub fn from_ranges(atom_ranges: Vec<(usize, usize)>, n_atoms: usize) -> Self {
        let merged_ranges = Group::fix_atom_ranges(atom_ranges, n_atoms);
        Group {
            atom_ranges: merged_ranges,
            print_ndx: true,
        }
    }

    /// Create a new valid Group structure from atom indices.
    ///
    /// ## Parameters
    /// Expects a vector of atom indices and the total number of atoms in the system.
    /// There can be duplicate atoms in the 'atom indices'. In the final Group structure, they will be removed.
    pub fn from_indices(atom_indices: Vec<usize>, n_atoms: usize) -> Self {
        let ranges = Group::make_atom_ranges(atom_indices, n_atoms);
        Group {
            atom_ranges: ranges,
            print_ndx: true,
        }
    }

    /// Create a new valid Group structure using Select tree.
    fn from_select(select: Box<Select>, system: &System) -> Result<Self, SelectError> {
        let mut indices = Vec::new();

        for i in 0usize..system.get_n_atoms() {
            if Group::matches_select(i, &select, system)? {
                indices.push(i);
            }
        }

        let ranges = Group::make_atom_ranges(indices, system.get_n_atoms());
        Ok(Group {
            atom_ranges: ranges,
            print_ndx: true,
        })
    }

    /// Consumes self returning a new `Group` which only contains atoms fullfilling the geometry condition.
    fn apply_geometry(self, geometry: impl Shape, system: &System) -> Self {
        let mut indices = Vec::new();

        let atoms = system.get_atoms_as_ref();
        let simbox = system.get_box_as_ref();

        // iterate through the atoms of the old group and select atoms fullfilling the geometry condition
        for (min, max) in self.atom_ranges {
            for (i, atom) in atoms.iter().enumerate().take(max + 1).skip(min) {
                if geometry.inside(atom.get_position(), simbox) {
                    indices.push(i);
                }
            }
        }

        let ranges = Group::make_atom_ranges(indices, system.get_n_atoms());

        Group {
            atom_ranges: ranges,
            print_ndx: true,
        }
    }

    /// Check whether properties of target atom match conditions prescribed by target Select tree.
    fn matches_select(
        atom_index: usize,
        select: &Select,
        system: &System,
    ) -> Result<bool, SelectError> {
        match select {
            Select::ResidueName(names) => Ok(names
                .iter()
                .any(|name| name == system.get_atoms_as_ref()[atom_index].get_residue_name())),

            Select::AtomName(names) => Ok(names
                .iter()
                .any(|name| name == system.get_atoms_as_ref()[atom_index].get_atom_name())),

            Select::ResidueNumber(numbers) => {
                let resnum = system.get_atoms_as_ref()[atom_index].get_residue_number();
                Ok(numbers
                    .iter()
                    .any(|&(start, end)| resnum >= start && resnum <= end))
            }

            Select::GmxAtomNumber(numbers) => Ok(numbers
                .iter()
                .any(|&(start, end)| atom_index + 1 >= start && atom_index < end)),

            Select::AtomNumber(numbers) => {
                let atomnum = system.get_atoms_as_ref()[atom_index].get_atom_number();
                Ok(numbers
                    .iter()
                    .any(|&(start, end)| atomnum >= start && atomnum <= end))
            }

            Select::Chain(identifiers) => {
                let chain = match system.get_atoms_as_ref()[atom_index].get_chain() {
                    None => return Ok(false),
                    Some(x) => x,
                };

                Ok(identifiers.iter().any(|&target| target == chain))
            }

            Select::GroupName(names) => {
                for name in names.iter() {
                    match system.group_isin(name.to_str(), atom_index) {
                        Ok(true) => return Ok(true),
                        Ok(false) => (),
                        // if the group does not exist, return an error
                        Err(_) => return Err(SelectError::GroupNotFound(name.to_str().to_owned())),
                    }
                }

                Ok(false)
            }

            Select::And(left, right) => Ok(Group::matches_select(atom_index, left, system)?
                && Group::matches_select(atom_index, right, system)?),

            Select::Or(left, right) => Ok(Group::matches_select(atom_index, left, system)?
                || Group::matches_select(atom_index, right, system)?),

            Select::Not(operand) => Ok(!Group::matches_select(atom_index, operand, system)?),
        }
    }

    /// Fix overlaps in atom ranges.
    /// Makes sure that the atom ranges are valid and not overflowing the number of atoms in the system.
    pub fn fix_atom_ranges(
        mut atom_ranges: Vec<(usize, usize)>,
        n_atoms: usize,
    ) -> Vec<(usize, usize)> {
        if atom_ranges.is_empty() {
            return atom_ranges;
        }

        // sort the atom ranges in ascending order
        atom_ranges.sort_unstable();

        let mut merged_indices = Vec::new();
        let mut current_start = std::usize::MAX;
        let mut current_end = 0usize;

        for (start, end) in &atom_ranges {
            // start must be smaller than n_atoms and not larger than end
            if *start >= n_atoms || *start > *end {
                continue;
            }

            // current range does not overlap with the previous one nor is adjacent to it
            if *start > current_end + 1 || (current_end == 0usize && current_start != 0usize) {
                if current_start != std::usize::MAX {
                    merged_indices.push((current_start, current_end.min(n_atoms - 1)));
                }
                current_start = *start;
                current_end = *end;
            // current range overlaps
            } else if *end > current_end {
                current_end = *end;
            }
        }

        if current_start != std::usize::MAX {
            // add the last merged range to the result if it exists
            merged_indices.push((current_start, current_end.min(n_atoms - 1)));
        }

        merged_indices
    }

    /// Create valid atom ranges from atom indices.
    /// Makes sure that the atom ranges are valid and not overflowing the number of atoms in the system.
    fn make_atom_ranges(mut atom_indices: Vec<usize>, n_atoms: usize) -> Vec<(usize, usize)> {
        let mut atom_ranges = Vec::new();

        if atom_indices.is_empty() {
            return atom_ranges;
        }

        // sort the atom indices in ascending order
        atom_indices.sort_unstable();

        let mut start = atom_indices[0];
        let mut end = atom_indices[0];

        for index in atom_indices.iter().skip(1) {
            // the range can't exceed the number of atoms
            if *index >= n_atoms {
                end = n_atoms - 1;
                break;
            }

            // remove duplicate atoms
            if *index == end {
                continue;
            }

            if *index == end + 1 {
                // extend the current range
                end = *index;
            } else {
                // add the completed range to the atom_indices vector and start a new range
                atom_ranges.push((start, end));
                start = *index;
                end = *index;
            }
        }

        // Add the last completed range to the atom_indices vector
        atom_ranges.push((start, end));

        atom_ranges
    }

    /// Get the number of atoms in the group.
    pub fn get_n_atoms(&self) -> usize {
        self.atom_ranges
            .iter()
            .fold(0, |acc, (start, end)| acc + (end - start + 1))
    }

    /// Writes an ndx file for a group
    pub fn write_ndx(&self, stream: &mut impl Write, name: &str) -> Result<(), WriteNdxError> {
        writeln!(stream, "[ {} ]", name).map_err(|_| WriteNdxError::CouldNotWrite)?;

        let last_index = match self.atom_ranges.last() {
            Some(x) => x.1,
            None => return Ok(()),
        };

        let mut iterator = 0usize;
        for (start, end) in self.atom_ranges.iter() {
            for index in *start..=*end {
                iterator += 1;
                if iterator % 15 == 0 || index == last_index {
                    writeln!(stream, "{:4}", index + 1)
                        .map_err(|_| WriteNdxError::CouldNotWrite)?;
                } else {
                    write!(stream, "{:4} ", index + 1).map_err(|_| WriteNdxError::CouldNotWrite)?;
                }
            }
        }

        Ok(())
    }

    /// Check whether the name for the group is a valid group name.
    /// Characters '"&|!@()<>= are not allowed. Names containing whitespace only are also not allowed.
    pub fn name_is_valid(string: &str) -> bool {
        if string.trim().is_empty() {
            return false;
        }

        let forbidden_chars = "'\"&|!@()<>=";

        for c in string.chars() {
            if forbidden_chars.contains(c) {
                return false;
            }
        }

        true
    }
}

/******************************/
/*    UNIT TESTS FOR GROUP    */
/******************************/

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_simple() {
        let atoms = vec![(20, 32)];

        let group = Group::from_ranges(atoms, 33);

        assert!(group.print_ndx);
        assert_eq!(group.atom_ranges[0], (20, 32));
    }

    #[test]
    fn new_ranges_empty() {
        let atoms = vec![];
        let group = Group::from_ranges(atoms, 1028);
        assert!(group.atom_ranges.is_empty());
    }

    #[test]
    fn new_multiple_nooverlap() {
        let atoms = vec![(20, 32), (64, 64), (84, 143)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (20, 32));
        assert_eq!(group.atom_ranges[1], (64, 64));
        assert_eq!(group.atom_ranges[2], (84, 143));
    }

    #[test]
    fn new_adjacent() {
        let atoms = vec![(20, 32), (33, 42)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (20, 42));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn new_overlap_simple() {
        let atoms = vec![(20, 32), (24, 42)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (20, 42));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn new_overlap_full() {
        let atoms = vec![(20, 32), (28, 30)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (20, 32));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn new_overlap_full_flipped() {
        let atoms = vec![(28, 30), (20, 32)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (20, 32));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn new_barely_overlaps() {
        let atoms = vec![(20, 32), (32, 42)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (20, 42));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn new_overlap_complex() {
        let atoms = vec![
            (64, 128),
            (5, 32),
            (1, 25),
            (129, 133),
            (133, 200),
            (35, 78),
            (10, 15),
        ];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (1, 32));
        assert_eq!(group.atom_ranges[1], (35, 200));
        assert_eq!(group.atom_ranges.len(), 2);
    }

    #[test]
    fn new_overlaps_with_zeros() {
        let atoms = vec![(1, 25), (0, 1), (0, 0), (0, 34)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (0, 34));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn new_start_larger_than_end() {
        let atoms = vec![(32, 25), (14, 17)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (14, 17));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn new_end_larger_than_natoms() {
        let atoms = vec![(543, 1020), (1000, 1432)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (543, 1027));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn new_start_larger_than_natoms() {
        let atoms = vec![(543, 1020), (1043, 1432)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (543, 1020));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn new_end_larger_than_natoms_nooverlap() {
        let atoms = vec![(0, 43), (1006, 1432)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (0, 43));
        assert_eq!(group.atom_ranges[1], (1006, 1027));
        assert_eq!(group.atom_ranges.len(), 2);
    }

    #[test]
    fn new_single_atoms() {
        let atoms = vec![
            (5, 5),
            (4, 4),
            (11, 11),
            (12, 12),
            (0, 0),
            (1, 1),
            (6, 6),
            (2, 2),
            (3, 3),
            (7, 7),
            (13, 13),
            (10, 10),
            (8, 8),
            (9, 9),
        ];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (0, 13));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn new_from_indices_basic() {
        let atom_indices = vec![6, 2, 13, 1, 10, 8, 3, 12, 7, 14, 15];
        let group = Group::from_indices(atom_indices, 16);

        assert!(group.print_ndx);

        assert_eq!(group.atom_ranges[0], (1, 3));
        assert_eq!(group.atom_ranges[1], (6, 8));
        assert_eq!(group.atom_ranges[2], (10, 10));
        assert_eq!(group.atom_ranges[3], (12, 15));
        assert_eq!(group.atom_ranges.len(), 4);
    }

    #[test]
    fn new_from_indices_empty() {
        let atoms = vec![];
        let group = Group::from_indices(atoms, 1028);
        assert!(group.atom_ranges.is_empty());
    }

    #[test]
    fn new_from_indices_duplicate() {
        let atom_indices = vec![1, 6, 3, 2, 13, 1, 10, 8, 3, 12, 7, 14, 15, 10];
        let group = Group::from_indices(atom_indices, 20);

        assert_eq!(group.atom_ranges[0], (1, 3));
        assert_eq!(group.atom_ranges[1], (6, 8));
        assert_eq!(group.atom_ranges[2], (10, 10));
        assert_eq!(group.atom_ranges[3], (12, 15));
        assert_eq!(group.atom_ranges.len(), 4);
    }

    #[test]
    fn new_from_indices_larger_than_natoms() {
        let atom_indices = vec![1, 6, 3, 2, 13, 1, 10, 8, 3, 12, 7, 14, 15, 10];
        let group = Group::from_indices(atom_indices, 15);

        assert_eq!(group.atom_ranges[0], (1, 3));
        assert_eq!(group.atom_ranges[1], (6, 8));
        assert_eq!(group.atom_ranges[2], (10, 10));
        assert_eq!(group.atom_ranges[3], (12, 14));
        assert_eq!(group.atom_ranges.len(), 4);
    }

    #[test]
    fn get_n_atoms_empty() {
        let atom_ranges = vec![];
        let group = Group::from_ranges(atom_ranges, 100);

        assert_eq!(group.get_n_atoms(), 0);
    }

    #[test]
    fn get_n_atoms_basic() {
        let atom_ranges = vec![
            (64, 128),
            (5, 32),
            (1, 25),
            (129, 133),
            (133, 200),
            (35, 78),
            (10, 15),
        ];
        let group = Group::from_ranges(atom_ranges, 1028);

        assert_eq!(group.get_n_atoms(), 198);
    }
}
