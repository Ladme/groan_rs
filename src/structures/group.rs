// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Implementation of the Group structure and its methods.

use std::io::Write;

use crate::errors::{SelectError, WriteNdxError};
use crate::select::Select;
use crate::structures::container::AtomContainer;
use crate::structures::shape::Shape;
use crate::system::iterating::get_molecule_indices;
use crate::system::System;

/******************************/
/*       GROUP STRUCTURE      */
/******************************/

/// Group of atoms in target system.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(deny_unknown_fields))]
pub struct Group {
    atoms: AtomContainer,
    pub print_ndx: bool,
}

impl Group {
    /// Get immutable reference to the `AtomContainer` of the Group.
    pub fn get_atoms(&self) -> &AtomContainer {
        &self.atoms
    }

    /// Create a new valid `Group` structure from query in Groan Selection Language.
    pub fn from_query(query: &str, system: &System) -> Result<Self, SelectError> {
        // parse groan selection language query into binary selection tree
        let select = Select::parse_query(query)?;
        // apply the selection tree to the system
        Group::from_select(*select, system)
    }

    /// Create a new valid `Group` structure from query in Groan Selection Language and from geometry constraint.
    pub fn from_query_and_geometry(
        query: &str,
        geometry: impl Shape,
        system: &System,
    ) -> Result<Self, SelectError> {
        // create group from query
        let group = Group::from_query(query, system)?;
        // now apply geometry constraint to the group
        Ok(group.apply_geometry(geometry, system))
    }

    /// Create a new valid `Group` structure from query in Groan Selection Language
    /// and from an arbitrary number of geometry constraints.
    pub fn from_query_and_geometries(
        query: &str,
        geometries: Vec<Box<dyn Shape>>,
        system: &System,
    ) -> Result<Self, SelectError> {
        // create group from query
        let group = Group::from_query(query, system)?;

        // now apply geometry constraints to the group
        Ok(group.apply_geometries(geometries, system))
    }

    /// Create a new valid Group structure from atom ranges.
    ///
    /// ## Parameters
    /// Expects a vector of atom ranges (start, end) and the total number of atoms in the system.
    /// The individual "atom ranges" in the atom ranges can be overlapping.
    /// In the final Group structure, they will be merged together.
    pub fn from_ranges(atom_ranges: Vec<(usize, usize)>, n_atoms: usize) -> Self {
        Group {
            atoms: AtomContainer::from_ranges(atom_ranges, n_atoms),
            print_ndx: true,
        }
    }

    /// Create a new valid Group structure from atom indices.
    ///
    /// ## Parameters
    /// Expects a vector of atom indices and the total number of atoms in the system.
    /// There can be duplicate atoms in the 'atom indices'. In the final Group structure, they will be removed.
    pub fn from_indices(atom_indices: Vec<usize>, n_atoms: usize) -> Self {
        Group {
            atoms: AtomContainer::from_indices(atom_indices, n_atoms),
            print_ndx: true,
        }
    }

    /// Create a new valid Group structure using Select tree.
    pub(crate) fn from_select(select: Select, system: &System) -> Result<Self, SelectError> {
        // expand regex group names and labels
        let select = Box::new(select.expand_regex_group_label(system)?);

        let indices: Vec<usize> = (0usize..system.get_n_atoms())
            .filter_map(|i| {
                if Group::matches_select(i, &select, system).ok()? {
                    Some(i)
                } else {
                    None
                }
            })
            .collect();

        Ok(Group {
            atoms: AtomContainer::from_indices(indices, system.get_n_atoms()),
            print_ndx: true,
        })
    }

    /// Consumes self returning a new `Group` which only contains atoms fullfilling the geometry condition.
    /// Atoms that have no positions are never inside any geometric shape.
    ///
    /// ## Panics
    /// Panics if the system has no simulation box.
    fn apply_geometry(self, geometry: impl Shape, system: &System) -> Self {
        let mut indices = Vec::new();

        let atoms = system.get_atoms();
        let simbox = system.get_box();

        // iterate through the atoms of the old group and select atoms fullfilling the geometry condition
        for index in self.atoms.iter() {
            // safety: `AtomContainerIterator` guarantees that we will only iterate through valid atoms
            let atom = unsafe { atoms.get_unchecked(index) };

            // atoms that have no positions are not inside the shape
            if let Some(pos) = atom.get_position() {
                if geometry.inside(pos, simbox.expect("FATAL GROAN ERROR | Group::apply_geometry | System does not have a simulation box.")) {
                    indices.push(index);
                }
            }
        }

        Group {
            atoms: AtomContainer::from_indices(indices, system.get_n_atoms()),
            print_ndx: true,
        }
    }

    /// Consumes self returning a new `Group` which only contains atoms fulfilling the specified geometry conditions.
    /// Atoms that have no positions are never inside any geometric shape.
    ///
    /// ## Panics
    /// Panics if the system has no simulation box.
    fn apply_geometries(self, geometries: Vec<Box<dyn Shape>>, system: &System) -> Self {
        let mut indices = Vec::new();

        let atoms = system.get_atoms();
        let simbox = system.get_box();

        // iterate through the atoms of the old group and select atoms fulfilling all geometry conditions
        for index in self.atoms.iter() {
            // safety: `AtomContainerIterator` guarantees that we will only iterate through valid atoms
            let atom = unsafe { atoms.get_unchecked(index) };
            let mut inside = true;

            for geom in &geometries {
                if let Some(pos) = atom.get_position() {
                    if !geom.inside(pos, simbox.expect("FATAL GROAN ERROR | Group::apply_geometries | System does not have a simulation box.")) {
                        inside = false;
                        break;
                    }
                // atoms that have no positions are not inside the shape
                } else {
                    inside = false;
                    break;
                }
            }

            if inside {
                indices.push(index);
            }
        }

        Group {
            atoms: AtomContainer::from_indices(indices, system.get_n_atoms()),
            print_ndx: true,
        }
    }

    /// Check whether properties of target atom match conditions prescribed by target Select tree.
    pub(super) fn matches_select(
        atom_index: usize,
        select: &Select,
        system: &System,
    ) -> Result<bool, SelectError> {
        match select {
            Select::ResidueName(names) => Ok(names
                .iter()
                .any(|name| name == system.get_atoms()[atom_index].get_residue_name())),

            Select::AtomName(names) => Ok(names
                .iter()
                .any(|name| name == system.get_atoms()[atom_index].get_atom_name())),

            Select::ResidueNumber(numbers) => {
                let resnum = system.get_atoms()[atom_index].get_residue_number();
                Ok(numbers
                    .iter()
                    .any(|&(start, end)| resnum >= start && resnum <= end))
            }

            Select::GmxAtomNumber(numbers) => Ok(numbers
                .iter()
                .any(|&(start, end)| atom_index + 1 >= start && atom_index < end)),

            Select::AtomNumber(numbers) => {
                let atomnum = system.get_atoms()[atom_index].get_atom_number();
                Ok(numbers
                    .iter()
                    .any(|&(start, end)| atomnum >= start && atomnum <= end))
            }

            Select::Chain(identifiers) => {
                let chain = match system.get_atoms()[atom_index].get_chain() {
                    None => return Ok(false),
                    Some(x) => x,
                };

                Ok(identifiers.iter().any(|&target| target == chain))
            }

            Select::GroupName(names) => {
                for name in names.iter() {
                    match name.match_groups(system, atom_index) {
                        Ok(true) => return Ok(true),
                        Ok(false) => (),
                        Err(e) => return Err(e),
                    }
                }

                Ok(false)
            }

            Select::LabeledAtom(names) => {
                for name in names.iter() {
                    match name.match_labels(system, atom_index) {
                        Ok(true) => return Ok(true),
                        Ok(false) => (),
                        Err(e) => return Err(e),
                    }
                }

                Ok(false)
            }

            Select::ElementName(names) => {
                let elname = match system.get_atoms()[atom_index].get_element_name() {
                    None => return Ok(false),
                    Some(x) => x,
                };

                Ok(names.iter().any(|target| target == elname))
            }

            Select::ElementSymbol(symbols) => {
                let elsymbol = match system.get_atoms()[atom_index].get_element_symbol() {
                    None => return Ok(false),
                    Some(x) => x,
                };

                Ok(symbols.iter().any(|target| target == elsymbol))
            }

            Select::Molecule(operand) => {
                let molecule_indices = get_molecule_indices(system, atom_index)
                    .expect("FATAL GROAN ERROR | Group::matches_select | Atom index should exist.");

                for index in molecule_indices {
                    if Group::matches_select(index, operand, system)? {
                        return Ok(true);
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

    /// Get the number of atoms in the group.
    pub fn get_n_atoms(&self) -> usize {
        self.atoms.get_n_atoms()
    }

    /// Writes an ndx file for a group
    pub fn write_ndx(&self, stream: &mut impl Write, name: &str) -> Result<(), WriteNdxError> {
        writeln!(stream, "[ {} ]", name).map_err(|_| WriteNdxError::CouldNotWrite)?;

        let last_index = match self.atoms.last() {
            Some(x) => x,
            None => return Ok(()),
        };

        for (iterator, index) in self.atoms.iter().enumerate() {
            if (iterator + 1) % 15 == 0 || index == last_index {
                writeln!(stream, "{:4}", index + 1).map_err(|_| WriteNdxError::CouldNotWrite)?;
            } else {
                write!(stream, "{:4} ", index + 1).map_err(|_| WriteNdxError::CouldNotWrite)?;
            }
        }

        Ok(())
    }

    /// Create a new valid `Group` as a union of two other groups.
    pub fn union(group1: &Group, group2: &Group) -> Group {
        let container = AtomContainer::union(group1.get_atoms(), group2.get_atoms());

        Group {
            atoms: container,
            print_ndx: true,
        }
    }

    /// Create a new valid `Group` as an interaction of two other groups.
    pub fn intersection(group1: &Group, group2: &Group) -> Group {
        let container = AtomContainer::intersection(group1.get_atoms(), group2.get_atoms());

        Group {
            atoms: container,
            print_ndx: true,
        }
    }
}

/******************************/
/*    UNIT TESTS FOR GROUP    */
/******************************/

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_ranges() {
        let ranges = vec![
            (10, 15),
            (17, 25),
            (11, 11),
            (7, 3),
            (9, 10),
            (15, 15),
            (16, 18),
            (2, 5),
            (10, 15),
        ];

        let group = Group::from_ranges(ranges, 20);

        assert_eq!(group.get_n_atoms(), 15);
        let expected_indices: Vec<usize> =
            vec![2, 3, 4, 5, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19];
        let observed_indices: Vec<usize> = group.atoms.iter().collect();
        assert_eq!(expected_indices, observed_indices);
    }

    #[test]
    fn from_indices() {
        let indices = vec![11, 1, 2, 3, 20, 5, 0, 5, 4, 18, 6, 19, 1, 13, 20, 27];
        let group = Group::from_indices(indices, 20);

        assert_eq!(group.get_n_atoms(), 11);
        let expected_indices: Vec<usize> = vec![0, 1, 2, 3, 4, 5, 6, 11, 13, 18, 19];
        let observed_indices: Vec<usize> = group.atoms.iter().collect();
        assert_eq!(expected_indices, observed_indices);
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

    #[test]
    fn union() {
        let indices1 = vec![11, 1, 2, 3, 20, 5, 0, 5, 4, 18, 6, 19, 1, 13, 20, 27];
        let group1 = Group::from_indices(indices1, 20);

        let indices2 = vec![13, 1, 2, 7, 5, 19, 21, 1, 9, 10, 11];
        let group2 = Group::from_indices(indices2, 15);

        let union = Group::union(&group1, &group2);
        assert_eq!(
            union,
            Group::from_indices(vec![0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 13, 14, 18, 19], 20)
        );
    }

    #[test]
    fn intersection() {
        let indices1 = vec![11, 1, 2, 3, 20, 5, 0, 5, 4, 18, 6, 19, 1, 13, 20, 27];
        let group1 = Group::from_indices(indices1, 20);

        let indices2 = vec![13, 1, 2, 7, 5, 19, 21, 1, 9, 10, 11];
        let group2 = Group::from_indices(indices2, 15);

        let union = Group::intersection(&group1, &group2);
        assert_eq!(union, Group::from_indices(vec![1, 2, 5, 11, 13], 20));
    }
}

#[cfg(test)]
#[cfg(feature = "serde")]
mod serde_tests {
    use super::*;

    #[test]
    fn group_to_yaml() {
        let ranges = vec![(20, 32), (64, 64), (84, 143)];
        let group = Group::from_ranges(ranges, 1028);

        let string = serde_yaml::to_string(&group).unwrap();
        let expected = "atoms:
  atom_blocks:
  - start: 20
    end: 32
  - start: 64
    end: 64
  - start: 84
    end: 143
print_ndx: true
";
        assert_eq!(string, expected);
    }

    #[test]
    fn group_from_yaml() {
        let string = "atoms:
  atom_blocks:
  - start: 20
    end: 32
  - start: 64
    end: 64
  - start: 84
    end: 143
print_ndx: true
";

        let group: Group = serde_yaml::from_str(string).unwrap();

        let ranges = vec![(20, 32), (64, 64), (84, 143)];
        let expected = Group::from_ranges(ranges, 1028);

        assert_eq!(group.atoms, expected.atoms);
        assert_eq!(group.print_ndx, expected.print_ndx);
    }
}
