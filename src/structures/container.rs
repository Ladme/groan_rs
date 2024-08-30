// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of the AtomContainer structure and its methods.

use std::cmp;

/// Structure describing a group of atoms.
/// Guaranteed to only contain valid atom indices.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(deny_unknown_fields))]
pub struct AtomContainer {
    /// Vector of atom blocks capturing the indices of atoms included in the container.
    atom_blocks: Vec<AtomBlock>,
}

/// Block of atoms described by atom indices in the System structure.
/// Guaranteed to only contain valid atom indices.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(feature = "serde", serde(deny_unknown_fields))]
struct AtomBlock {
    /// Index of the first atom in the `AtomBlock`.
    /// Atoms are indexed starting from 0.
    start: usize,
    /// Index of the last atom in the `AtomBlock`.
    /// Atoms are indexed starting from 0.
    /// The atom with `max` index will also be included in the `AtomBlock`.
    end: usize,
}

impl AtomContainer {
    /// Create an empty `AtomContainer`.
    pub fn empty() -> Self {
        AtomContainer {
            atom_blocks: Vec::new(),
        }
    }

    /// Create a new valid `AtomContainer` structure from individual indices of the atoms.
    ///
    /// ## Parameters
    /// - `atom_indices`: vector of atom indices from the `System` structure. Atoms are indexed from 0.
    /// - `n_atoms`: the total number of atoms in the `System` structure
    ///
    /// ## Notes
    /// - `atom_indices` structure can contain duplicate indices.
    /// - In case an index in the `atom_indices` is out of range (>= `n_atoms`),
    ///   the index of the last atom in the `System` is used instead.
    pub fn from_indices(mut atom_indices: Vec<usize>, n_atoms: usize) -> Self {
        let mut blocks: Vec<AtomBlock> = Vec::new();

        if atom_indices.is_empty() {
            return AtomContainer {
                atom_blocks: blocks,
            };
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
                // create the block
                unsafe {
                    // safety: we previously make sure that `end` is lower than `n_atoms`
                    // start can not be higher than end from logic
                    blocks.push(AtomBlock::new_unchecked(start, end));
                }
                // start new block
                start = *index;
                end = *index;
            }
        }

        // add the last completed block to the `blocks` vector
        unsafe {
            // safety: we previously make sure that `end` is lower than `n_atoms`
            // start can not be higher than end from logic
            blocks.push(AtomBlock::new_unchecked(start, end));
        }

        AtomContainer {
            atom_blocks: blocks,
        }
    }

    /// Create a new valid `AtomContainer` structure from atom ranges.
    ///
    /// ## Parameters
    /// - `atom_ranges`: vector of atom ranges from the `System` structure. Atoms are indexed from 0.
    ///   Atom range is a tuple representation of an `AtomBlock`. The first element of the tuple
    ///   corresponds to the index of the first atom in the block. The second element corresponds
    ///   to the index of the last atom in the block.
    ///
    /// - `n_atoms`: the total number of atoms in the `System` structure
    ///
    /// ## Notes
    /// - `atom_ranges` may be overlapping and multiple identical atom ranges may be provided.
    /// - The `start` value of the ranges may be higher than the `end` value. Such ranges are ignored.
    /// - `atom_ranges` may contain atom indices that are out of range (do not exist in the `System`).
    ///   In such cases, the index of the last atom in the `System` is used instead as the upper bound
    ///   of the range.
    pub fn from_ranges(atom_ranges: Vec<(usize, usize)>, n_atoms: usize) -> Self {
        // create `AtomBlocks` from the ranges
        let mut blocks = Vec::new();

        if n_atoms == 0 {
            return AtomContainer {
                atom_blocks: blocks,
            };
        }

        for range in atom_ranges {
            let start = range.0;
            let end = if range.1 < n_atoms {
                range.1
            } else {
                n_atoms - 1
            };

            // ignore ranges in which the start is higher than end
            if start > end {
                continue;
            }

            unsafe {
                // safety: we make sure that the `AtomBlock` is valid in previous steps
                blocks.push(AtomBlock::new_unchecked(start, end));
            }
        }

        // then create the AtomContainer
        AtomContainer::from_blocks(blocks)
    }

    /// Get the number of atoms in the `AtomContainer`.
    ///
    /// ## Warning
    /// This is not an O(1) operation!
    /// Time complexity of this operation depends on the complexity of the `AtomContainer`,
    /// i.e. on the number of `AtomBlocks` forming the `AtomContainer`.
    pub fn get_n_atoms(&self) -> usize {
        self.atom_blocks
            .iter()
            .fold(0, |acc, block| acc + block.get_n_atoms())
    }

    /// Create a new valid `AtomContainer` from vector of `AtomBlocks`.
    /// The `AtomBlocks` may be overlapping.
    fn from_blocks(mut atom_blocks: Vec<AtomBlock>) -> Self {
        if atom_blocks.is_empty() {
            return AtomContainer {
                atom_blocks: Vec::new(),
            };
        }

        // sort the blocks in ascending order
        atom_blocks.sort_unstable();

        let mut new_blocks = Vec::new();
        let mut current_start = usize::MAX;
        let mut current_end = 0usize;

        for block in atom_blocks {
            let start = block.start;
            let end = block.end;

            // current block does not overlap with the previous one nor is adjacent to it
            if block.start > current_end + 1 || (current_end == 0 && current_start != 0) {
                if current_start != usize::MAX {
                    unsafe {
                        // safety: each `block` is guaranteed to be valid; therefore, new `AtomBlock` must be valid as well
                        new_blocks.push(AtomBlock::new_unchecked(current_start, current_end));
                    }
                }

                current_start = start;
                current_end = end;
            // current range overlaps
            } else if end > current_end {
                current_end = end;
            }
        }

        // add the last merged block to the result if it exists
        if current_start != usize::MAX {
            unsafe {
                // safety: each `block` is guaranteed to be valid; therefore, new `AtomBlock` must be valid as well
                new_blocks.push(AtomBlock::new_unchecked(current_start, current_end));
            }
        }

        AtomContainer {
            atom_blocks: new_blocks,
        }
    }

    /// Iterate over the atom indices of the `AtomContainer`.
    pub fn iter(&self) -> AtomContainerIterator {
        AtomContainerIterator {
            container: self,
            current_block_index: 0,
            current_atom_index: 0,
        }
    }

    /// Convert `AtomContainer` into `OwnedAtomContainerIterator`.
    pub fn into_iter(self) -> OwnedAtomContainerIterator {
        OwnedAtomContainerIterator {
            container: self,
            current_block_index: 0,
            current_atom_index: 0,
        }
    }

    /// Get the index of the last atom in the `AtomContainer`.
    /// Returns `None` if the `AtomContainer` is empty.
    pub fn last(&self) -> Option<usize> {
        self.atom_blocks.last().map(|block| block.end)
    }

    /// Get the index of the first atom in the `AtomContainer`.
    /// Returns `None` if the `AtomContainer` is empty.
    pub fn first(&self) -> Option<usize> {
        self.atom_blocks.first().map(|block| block.start)
    }

    /// Returns `true` if the `AtomContainer` contains no atoms. Else returns `false`.
    pub fn is_empty(&self) -> bool {
        self.atom_blocks.is_empty()
    }

    /// Returns `true` if the `AtomContainer` contains the atom with target index.
    pub fn isin(&self, index: usize) -> bool {
        for block in &self.atom_blocks {
            // blocks are ordered; if index is smaller than the start of the current block,
            // it cannot be present in the container
            if index < block.start {
                return false;
            }

            // we only need to test for block.end because we already test for block.start before
            if index <= block.end {
                return true;
            }
        }

        false
    }

    /// Create a new `AtomContainer` that is the union of the provided `AtomContainers`.
    pub fn union(container1: &AtomContainer, container2: &AtomContainer) -> AtomContainer {
        let blocks = container1
            .atom_blocks
            .iter()
            .cloned()
            .chain(container2.atom_blocks.iter().cloned())
            .collect();
        AtomContainer::from_blocks(blocks)
    }

    /// Create a new `AtomContainer` that is the intersection of the provided `AtomContainers`.
    pub fn intersection(container1: &AtomContainer, container2: &AtomContainer) -> AtomContainer {
        let indices: Vec<usize> = container1
            .iter()
            .filter(|&index| container2.isin(index))
            .collect();

        let n_atoms = match indices.iter().max() {
            Some(x) => x + 1,
            None => return AtomContainer::empty(),
        };

        AtomContainer::from_indices(indices, n_atoms)
    }

    /// Add index to the `AtomContainer`.
    /// Maintains the validity of the container.
    /// If `index` is higher than or equal to `n_atoms`, it is not added to the container.
    pub fn add(&mut self, index: usize, n_atoms: usize) {
        if index >= n_atoms {
            return;
        }

        let mut new_blocks = self.atom_blocks.clone();
        // safety: we check that index is smaller than `n_atoms` in a previous step
        new_blocks.push(unsafe { AtomBlock::new_unchecked(index, index) });

        self.atom_blocks = AtomContainer::from_blocks(new_blocks).atom_blocks;
    }
}

impl AtomBlock {
    /// Create a new valid `AtomBlock` structure.
    ///
    /// ## Parameters
    /// - `start`: index of the first atom in the `AtomBlock`. Atoms are indexed from 0.
    /// - `end`: index of the last atom in the `AtomBlock`. Atoms are indexed from 0.
    /// - `n_atoms` the total number of atoms in the `System` structure.
    ///
    /// ## Panics
    /// - Panics if `end` is higher than or equal to `n_atoms`.
    /// - Panics if `start` is `higher` than `end`.
    #[allow(dead_code)]
    fn new(start: usize, end: usize, n_atoms: usize) -> Self {
        if end >= n_atoms {
            panic!(
                "FATAL GROAN ERROR | AtomBlock::new | `AtomBlock` could not be constructed: `end` is out of range."
            );
        }

        if start > end {
            panic!(
                "FATAL GROAN ERROR | AtomBlock::new | `AtomBlock` could not be constructed: `start` is higher than `end`."
            );
        }

        AtomBlock { start, end }
    }

    /// Create a new `AtomBlock` structure without checking its validity.
    ///
    /// ## Parameters
    /// - `start`: index of the first atom in the `AtomBlock`. Atoms are indexed from 0.
    /// - `end`: index of the last atom in the `AtomBlock`. Atoms are indexed from 0.
    ///
    /// ## Safety
    /// This function is only safe to use if:
    /// - a) `end` is lower than the number of atoms in the `System` and
    /// - b) `start` is lower than or equal to `end`.
    unsafe fn new_unchecked(start: usize, end: usize) -> Self {
        AtomBlock { start, end }
    }

    /// Get the number of atoms in the `AtomBlock`.
    fn get_n_atoms(&self) -> usize {
        self.end - self.start + 1
    }
}

/// Implementation of the `Ord` trait for `AtomBlock`.
impl cmp::Ord for AtomBlock {
    /// Comparison of `AtomBlock` is performed in the same way as for 2-member tuples.
    fn cmp(&self, other: &Self) -> cmp::Ordering {
        match self.start.cmp(&other.start) {
            cmp::Ordering::Equal => self.end.cmp(&other.end),
            ordering => ordering,
        }
    }
}

/// Implementation of the `PartialOrd` trait for `AtomBlock`.
impl cmp::PartialOrd for AtomBlock {
    fn partial_cmp(&self, other: &Self) -> Option<cmp::Ordering> {
        Some(self.cmp(other))
    }
}

// Private helper function to handle the iteration logic
#[inline(always)]
fn next_index(
    container: &AtomContainer,
    current_block_index: &mut usize,
    current_atom_index: &mut usize,
) -> Option<usize> {
    while *current_block_index < container.atom_blocks.len() {
        let block = unsafe { &container.atom_blocks.get_unchecked(*current_block_index) };

        // move current_atom_index to the start of the block if it's not there yet
        if *current_atom_index < block.start {
            // immediately move to the start + 1 position
            *current_atom_index = block.start + 1;
            // and return the start position
            return Some(block.start);
        }

        // check if current_atom_index is within the block range
        if *current_atom_index <= block.end {
            let index_to_return = *current_atom_index;
            *current_atom_index += 1;

            // return the current index
            return Some(index_to_return);
        }

        // move to the next block
        *current_block_index += 1;
    }

    None
}

/// Iterator over the atom indices in the `AtomContainer`.
#[derive(Debug, Clone)]
pub struct AtomContainerIterator<'a> {
    container: &'a AtomContainer,
    current_block_index: usize,
    current_atom_index: usize,
}

/// Implementation of the iteration over indices of the `AtomContainer`.
impl<'a> Iterator for AtomContainerIterator<'a> {
    type Item = usize;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        next_index(
            self.container,
            &mut self.current_block_index,
            &mut self.current_atom_index,
        )
    }
}

/// Iterator over the atom indices in the `AtomContainer`
/// which owns the `AtomContainer`.
#[derive(Debug, Clone)]
pub struct OwnedAtomContainerIterator {
    container: AtomContainer,
    current_block_index: usize,
    current_atom_index: usize,
}

/// Implementation of the iteration over the indicies of the `AtomContainer`.
impl Iterator for OwnedAtomContainerIterator {
    type Item = usize;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        next_index(
            &self.container,
            &mut self.current_block_index,
            &mut self.current_atom_index,
        )
    }
}

#[cfg(test)]
mod tests_block {
    use super::*;

    #[test]
    fn new() {
        let block = AtomBlock::new(14, 22, 30);

        assert_eq!(block.start, 14);
        assert_eq!(block.end, 22);

        let block = AtomBlock::new(14, 22, 23);

        assert_eq!(block.start, 14);
        assert_eq!(block.end, 22);

        unsafe {
            let block = AtomBlock::new_unchecked(14, 22);

            assert_eq!(block.start, 14);
            assert_eq!(block.end, 22);
        }
    }

    #[test]
    #[should_panic(
        expected = "FATAL GROAN ERROR | AtomBlock::new | `AtomBlock` could not be constructed: `end` is out of range."
    )]
    fn new_panic_1() {
        let _block = AtomBlock::new(14, 22, 22);
    }

    #[test]
    #[should_panic(
        expected = "FATAL GROAN ERROR | AtomBlock::new | `AtomBlock` could not be constructed: `start` is higher than `end`."
    )]
    fn new_panic_2() {
        let _block = AtomBlock::new(22, 14, 30);
    }

    #[test]
    fn get_n_atoms() {
        let block = AtomBlock::new(5, 5, 20);
        assert_eq!(block.get_n_atoms(), 1);

        let block = AtomBlock::new(0, 11, 20);
        assert_eq!(block.get_n_atoms(), 12);
    }
}

#[cfg(test)]
mod tests_container {
    use super::*;

    fn cmp_block_tuple(block: &AtomBlock, tuple: (usize, usize)) {
        assert_eq!(block.start, tuple.0);
        assert_eq!(block.end, tuple.1);
    }

    #[test]
    fn from_indices_basic() {
        let indices = vec![6, 2, 13, 1, 10, 8, 3, 12, 7, 14, 15];
        let container = AtomContainer::from_indices(indices, 16);

        assert_eq!(container.atom_blocks.len(), 4);
        cmp_block_tuple(&container.atom_blocks[0], (1, 3));
        cmp_block_tuple(&container.atom_blocks[1], (6, 8));
        cmp_block_tuple(&container.atom_blocks[2], (10, 10));
        cmp_block_tuple(&container.atom_blocks[3], (12, 15));
    }

    #[test]
    fn from_indices_empty() {
        let indices = vec![];
        let container = AtomContainer::from_indices(indices, 16);

        assert!(container.is_empty());
    }

    #[test]
    fn from_indices_duplicate() {
        let indices = vec![1, 6, 3, 2, 13, 1, 10, 8, 3, 12, 7, 14, 15, 10];
        let container = AtomContainer::from_indices(indices, 16);

        assert_eq!(container.atom_blocks.len(), 4);
        cmp_block_tuple(&container.atom_blocks[0], (1, 3));
        cmp_block_tuple(&container.atom_blocks[1], (6, 8));
        cmp_block_tuple(&container.atom_blocks[2], (10, 10));
        cmp_block_tuple(&container.atom_blocks[3], (12, 15));
    }

    #[test]
    fn from_indices_larger_than_natoms() {
        let indices = vec![1, 6, 3, 2, 13, 1, 10, 8, 3, 12, 7, 14, 15, 10];
        let container = AtomContainer::from_indices(indices, 15);

        assert_eq!(container.atom_blocks.len(), 4);
        cmp_block_tuple(&container.atom_blocks[0], (1, 3));
        cmp_block_tuple(&container.atom_blocks[1], (6, 8));
        cmp_block_tuple(&container.atom_blocks[2], (10, 10));
        cmp_block_tuple(&container.atom_blocks[3], (12, 14));
    }

    #[test]
    fn from_indices_complex() {
        let indices = vec![11, 1, 2, 3, 20, 5, 0, 5, 4, 18, 6, 19, 1, 13, 20, 27];
        let container = AtomContainer::from_indices(indices, 20);

        assert_eq!(container.atom_blocks.len(), 4);
        cmp_block_tuple(&container.atom_blocks[0], (0, 6));
        cmp_block_tuple(&container.atom_blocks[1], (11, 11));
        cmp_block_tuple(&container.atom_blocks[2], (13, 13));
        cmp_block_tuple(&container.atom_blocks[3], (18, 19));
    }

    #[test]
    fn from_blocks() {
        let blocks = vec![
            AtomBlock::new(10, 15, 20),
            AtomBlock::new(17, 19, 20),
            AtomBlock::new(11, 11, 20),
            AtomBlock::new(9, 10, 20),
            AtomBlock::new(15, 15, 20),
            AtomBlock::new(16, 18, 20),
            AtomBlock::new(2, 5, 20),
            AtomBlock::new(10, 15, 20),
        ];

        let container = AtomContainer::from_blocks(blocks);

        assert_eq!(container.atom_blocks.len(), 2);
        cmp_block_tuple(&container.atom_blocks[0], (2, 5));
        cmp_block_tuple(&container.atom_blocks[1], (9, 19));
    }

    #[test]
    fn from_ranges_simple() {
        let ranges = vec![(20, 32)];
        let container = AtomContainer::from_ranges(ranges, 33);

        assert_eq!(container.atom_blocks.len(), 1);
        cmp_block_tuple(&container.atom_blocks[0], (20, 32));
    }

    #[test]
    fn from_ranges_empty() {
        let ranges = vec![];
        let container = AtomContainer::from_ranges(ranges, 33);

        assert!(container.is_empty());
    }

    #[test]
    fn from_ranges_multiple_nooverlap() {
        let ranges = vec![(20, 32), (64, 64), (84, 143)];
        let container = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container.atom_blocks.len(), 3);
        cmp_block_tuple(&container.atom_blocks[0], (20, 32));
        cmp_block_tuple(&container.atom_blocks[1], (64, 64));
        cmp_block_tuple(&container.atom_blocks[2], (84, 143));
    }

    #[test]
    fn from_ranges_adjacent() {
        let ranges = vec![(20, 32), (33, 42)];
        let container = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container.atom_blocks.len(), 1);
        cmp_block_tuple(&container.atom_blocks[0], (20, 42));
    }

    #[test]
    fn from_ranges_overlap_simple() {
        let ranges = vec![(20, 32), (24, 42)];
        let container = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container.atom_blocks.len(), 1);
        cmp_block_tuple(&container.atom_blocks[0], (20, 42));
    }

    #[test]
    fn from_ranges_overlap_alt() {
        let ranges = vec![(20, 35), (20, 32)];
        let container = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container.atom_blocks.len(), 1);
        cmp_block_tuple(&container.atom_blocks[0], (20, 35));
    }

    #[test]
    fn from_ranges_overlap_full() {
        let ranges = vec![(20, 32), (28, 30)];
        let container = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container.atom_blocks.len(), 1);
        cmp_block_tuple(&container.atom_blocks[0], (20, 32));
    }

    #[test]
    fn from_ranges_overlap_full_flipped() {
        let ranges = vec![(28, 30), (20, 32)];
        let container = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container.atom_blocks.len(), 1);
        cmp_block_tuple(&container.atom_blocks[0], (20, 32));
    }

    #[test]
    fn from_ranges_barely_overlaps() {
        let ranges = vec![(20, 32), (32, 42)];
        let container = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container.atom_blocks.len(), 1);
        cmp_block_tuple(&container.atom_blocks[0], (20, 42));
    }

    #[test]
    fn from_ranges_overlap_complex() {
        let ranges = vec![
            (64, 128),
            (5, 32),
            (1, 25),
            (129, 133),
            (133, 200),
            (35, 78),
            (10, 15),
            (1033, 1055),
        ];

        let container = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container.atom_blocks.len(), 2);
        cmp_block_tuple(&container.atom_blocks[0], (1, 32));
        cmp_block_tuple(&container.atom_blocks[1], (35, 200));
    }

    #[test]
    fn from_ranges_overlaps_with_zeros() {
        let ranges = vec![(1, 25), (0, 1), (0, 0), (0, 34)];
        let container = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container.atom_blocks.len(), 1);
        cmp_block_tuple(&container.atom_blocks[0], (0, 34));
    }

    #[test]
    fn from_ranges_start_larger_than_end() {
        let ranges = vec![(32, 25), (14, 17)];
        let container = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container.atom_blocks.len(), 1);
        cmp_block_tuple(&container.atom_blocks[0], (14, 17));
    }

    #[test]
    fn from_ranges_end_larger_than_natoms() {
        let ranges = vec![(543, 1020), (1000, 1432)];
        let container = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container.atom_blocks.len(), 1);
        cmp_block_tuple(&container.atom_blocks[0], (543, 1027));
    }

    #[test]
    fn from_ranges_start_larger_than_natoms() {
        let ranges = vec![(543, 1020), (1043, 1432)];
        let container = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container.atom_blocks.len(), 1);
        cmp_block_tuple(&container.atom_blocks[0], (543, 1020));
    }

    #[test]
    fn from_ranges_end_larger_than_natoms_nooverlap() {
        let ranges = vec![(0, 43), (1006, 1432)];
        let container = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container.atom_blocks.len(), 2);
        cmp_block_tuple(&container.atom_blocks[0], (0, 43));
        cmp_block_tuple(&container.atom_blocks[1], (1006, 1027));
    }

    #[test]
    fn from_ranges_single_atoms() {
        let ranges = vec![
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

        let container = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container.atom_blocks.len(), 1);
        cmp_block_tuple(&container.atom_blocks[0], (0, 13));
    }

    #[test]
    fn from_ranges_complex() {
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

        let container = AtomContainer::from_ranges(ranges, 20);

        assert_eq!(container.atom_blocks.len(), 2);
        cmp_block_tuple(&container.atom_blocks[0], (2, 5));
        cmp_block_tuple(&container.atom_blocks[1], (9, 19));
    }

    #[test]
    fn get_n_atoms() {
        let indices = vec![11, 1, 2, 3, 20, 5, 0, 5, 4, 18, 6, 19, 1, 13, 20, 27];
        let container = AtomContainer::from_indices(indices, 20);

        assert_eq!(container.get_n_atoms(), 11);
    }

    #[test]
    fn iterate() {
        let indices = vec![11, 1, 2, 3, 20, 5, 0, 5, 4, 18, 6, 19, 1, 13, 20, 27];
        let container = AtomContainer::from_indices(indices, 20);

        let expected_indices = vec![0, 1, 2, 3, 4, 5, 6, 11, 13, 18, 19];
        let observed_indices: Vec<usize> = container.iter().collect();

        assert_eq!(expected_indices, observed_indices);
    }

    #[test]
    fn into_iter() {
        let indices = vec![11, 1, 2, 3, 20, 5, 0, 5, 4, 18, 6, 19, 1, 13, 20, 27];
        let container = AtomContainer::from_indices(indices, 20);

        let expected_indices = vec![0, 1, 2, 3, 4, 5, 6, 11, 13, 18, 19];
        let observed_indices: Vec<usize> = container.into_iter().collect();

        assert_eq!(expected_indices, observed_indices);
    }

    #[test]
    fn last() {
        let indices = vec![11, 1, 2, 3, 20, 5, 0, 5, 4, 18, 6, 19, 1, 13, 20, 27];
        let container = AtomContainer::from_indices(indices, 20);

        assert_eq!(container.last().unwrap(), 19);

        let indices = vec![5, 5];
        let container = AtomContainer::from_indices(indices, 20);

        assert_eq!(container.last().unwrap(), 5);

        let indices: Vec<usize> = Vec::new();
        let container = AtomContainer::from_indices(indices, 20);

        assert_eq!(container.last(), None);
    }

    #[test]
    fn first() {
        let indices = vec![11, 1, 2, 3, 20, 5, 0, 5, 4, 18, 6, 19, 1, 13, 20, 27];
        let container = AtomContainer::from_indices(indices, 20);

        assert_eq!(container.first().unwrap(), 0);

        let indices = vec![5, 5];
        let container = AtomContainer::from_indices(indices, 20);

        assert_eq!(container.first().unwrap(), 5);

        let indices: Vec<usize> = Vec::new();
        let container = AtomContainer::from_indices(indices, 20);

        assert_eq!(container.first(), None);
    }

    #[test]
    fn isin() {
        let indices = vec![11, 1, 2, 3, 20, 5, 0, 5, 4, 18, 6, 19, 1, 13, 20, 27];
        let container = AtomContainer::from_indices(indices, 20);

        assert!(container.isin(5));
        assert!(!container.isin(12));
        assert!(!container.isin(73));
    }

    #[test]
    fn isin2() {
        let indices = vec![11, 20, 5, 5, 4, 18, 6, 19, 13, 20, 27];
        let container = AtomContainer::from_indices(indices, 20);

        assert!(!container.isin(1));
        assert!(!container.isin(3));
        assert!(container.isin(5));
    }

    #[test]
    fn union() {
        let indices = vec![11, 1, 2, 3, 20, 5, 0, 5, 4, 18, 6, 19, 1, 13, 20, 27];
        let container1 = AtomContainer::from_indices(indices, 20);

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
        let container2 = AtomContainer::from_ranges(ranges, 20);

        let union = AtomContainer::union(&container1, &container2);

        //vec![0, 1, 2, 3, 4, 5, 6, 11, 13, 18, 19];
        assert_eq!(union.atom_blocks.len(), 2);
        cmp_block_tuple(&union.atom_blocks[0], (0, 6));
        cmp_block_tuple(&union.atom_blocks[1], (9, 19));
    }

    #[test]
    fn intersection() {
        let indices = vec![11, 1, 2, 3, 20, 5, 0, 5, 4, 18, 6, 19, 1, 13, 20, 27];
        let container1 = AtomContainer::from_indices(indices, 20);
        let container2 =
            AtomContainer::from_indices(vec![13, 1, 2, 7, 5, 19, 21, 1, 9, 10, 11], 15);

        let intersect1 = AtomContainer::intersection(&container1, &container2);
        let intersect2 = AtomContainer::intersection(&container2, &container1);

        assert_eq!(intersect1.atom_blocks.len(), 4);
        cmp_block_tuple(&intersect1.atom_blocks[0], (1, 2));
        cmp_block_tuple(&intersect1.atom_blocks[1], (5, 5));
        cmp_block_tuple(&intersect1.atom_blocks[2], (11, 11));
        cmp_block_tuple(&intersect1.atom_blocks[3], (13, 13));
        assert_eq!(intersect1, intersect2);
    }

    #[test]
    fn intersection_empty() {
        let indices = vec![11, 1, 2, 3, 20, 5, 0, 5, 4, 18, 6, 19, 1, 13, 20, 27];
        let container1 = AtomContainer::from_indices(indices, 20);
        let container2 = AtomContainer::empty();

        let intersect1 = AtomContainer::intersection(&container1, &container2);
        let intersect2 = AtomContainer::intersection(&container2, &container1);

        assert_eq!(intersect1, AtomContainer::empty());
        assert_eq!(intersect1, intersect2);
    }

    #[test]
    fn intersection_disjunct() {
        let indices = vec![11, 1, 2, 3, 20, 5, 0, 5, 4, 18, 6, 19, 1, 13, 20, 27];
        let container1 = AtomContainer::from_indices(indices, 20);
        let container2 = AtomContainer::from_indices(vec![7, 8, 9, 10, 12, 14], 15);

        let intersect1 = AtomContainer::intersection(&container1, &container2);
        let intersect2 = AtomContainer::intersection(&container2, &container1);

        assert_eq!(intersect1, AtomContainer::empty());
        assert_eq!(intersect1, intersect2);
    }

    #[test]
    fn add() {
        let indices = vec![15, 2, 3, 5, 6, 7, 10, 11, 12, 15];
        let mut container = AtomContainer::from_indices(indices, 20);

        container.add(18, 20);
        assert_eq!(container.atom_blocks.len(), 5);
        cmp_block_tuple(&container.atom_blocks[0], (2, 3));
        cmp_block_tuple(&container.atom_blocks[1], (5, 7));
        cmp_block_tuple(&container.atom_blocks[2], (10, 12));
        cmp_block_tuple(&container.atom_blocks[3], (15, 15));
        cmp_block_tuple(&container.atom_blocks[4], (18, 18));

        container.add(1, 20);
        assert_eq!(container.atom_blocks.len(), 5);
        cmp_block_tuple(&container.atom_blocks[0], (1, 3));
        cmp_block_tuple(&container.atom_blocks[1], (5, 7));
        cmp_block_tuple(&container.atom_blocks[2], (10, 12));
        cmp_block_tuple(&container.atom_blocks[3], (15, 15));
        cmp_block_tuple(&container.atom_blocks[4], (18, 18));

        container.add(8, 20);
        assert_eq!(container.atom_blocks.len(), 5);
        cmp_block_tuple(&container.atom_blocks[0], (1, 3));
        cmp_block_tuple(&container.atom_blocks[1], (5, 8));
        cmp_block_tuple(&container.atom_blocks[2], (10, 12));
        cmp_block_tuple(&container.atom_blocks[3], (15, 15));
        cmp_block_tuple(&container.atom_blocks[4], (18, 18));

        container.add(9, 20);
        assert_eq!(container.atom_blocks.len(), 4);
        cmp_block_tuple(&container.atom_blocks[0], (1, 3));
        cmp_block_tuple(&container.atom_blocks[1], (5, 12));
        cmp_block_tuple(&container.atom_blocks[2], (15, 15));
        cmp_block_tuple(&container.atom_blocks[3], (18, 18));

        container.add(19, 20);
        assert_eq!(container.atom_blocks.len(), 4);
        cmp_block_tuple(&container.atom_blocks[0], (1, 3));
        cmp_block_tuple(&container.atom_blocks[1], (5, 12));
        cmp_block_tuple(&container.atom_blocks[2], (15, 15));
        cmp_block_tuple(&container.atom_blocks[3], (18, 19));

        container.add(16, 20);
        container.add(14, 20);
        assert_eq!(container.atom_blocks.len(), 4);
        cmp_block_tuple(&container.atom_blocks[0], (1, 3));
        cmp_block_tuple(&container.atom_blocks[1], (5, 12));
        cmp_block_tuple(&container.atom_blocks[2], (14, 16));
        cmp_block_tuple(&container.atom_blocks[3], (18, 19));

        container.add(4, 20);
        assert_eq!(container.atom_blocks.len(), 3);
        cmp_block_tuple(&container.atom_blocks[0], (1, 12));
        cmp_block_tuple(&container.atom_blocks[1], (14, 16));
        cmp_block_tuple(&container.atom_blocks[2], (18, 19));

        container.add(20, 20);
        assert_eq!(container.atom_blocks.len(), 3);
        cmp_block_tuple(&container.atom_blocks[0], (1, 12));
        cmp_block_tuple(&container.atom_blocks[1], (14, 16));
        cmp_block_tuple(&container.atom_blocks[2], (18, 19));

        container.add(25, 20);
        assert_eq!(container.atom_blocks.len(), 3);
        cmp_block_tuple(&container.atom_blocks[0], (1, 12));
        cmp_block_tuple(&container.atom_blocks[1], (14, 16));
        cmp_block_tuple(&container.atom_blocks[2], (18, 19));

        container.add(7, 20);
        assert_eq!(container.atom_blocks.len(), 3);
        cmp_block_tuple(&container.atom_blocks[0], (1, 12));
        cmp_block_tuple(&container.atom_blocks[1], (14, 16));
        cmp_block_tuple(&container.atom_blocks[2], (18, 19));
    }
}

#[cfg(test)]
#[cfg(feature = "serde")]
mod serde_tests {
    use super::*;

    #[test]
    fn atomblock_to_yaml() {
        let block = AtomBlock::new(14, 22, 30);

        let string = serde_yaml::to_string(&block).unwrap();

        assert_eq!(string, "start: 14\nend: 22\n");
    }

    #[test]
    fn atomblock_from_yaml() {
        let string = "{ start: 14, end: 22}";

        let block: AtomBlock = serde_yaml::from_str(&string).unwrap();

        assert_eq!(block.start, 14);
        assert_eq!(block.end, 22);
    }

    #[test]
    fn atomcontainer_to_yaml() {
        let ranges = vec![(20, 32), (64, 64), (84, 143)];
        let container = AtomContainer::from_ranges(ranges, 1028);

        let string = serde_yaml::to_string(&container).unwrap();

        assert_eq!(
            string,
            "atom_blocks:
- start: 20
  end: 32
- start: 64
  end: 64
- start: 84
  end: 143
"
        )
    }

    #[test]
    fn atomcontainer_from_yaml() {
        let string = "atom_blocks: [ { start: 20, end: 32 }, { start: 64, end: 64 }, { start: 84, end: 143 } ]";

        let container: AtomContainer = serde_yaml::from_str(&string).unwrap();

        let ranges = vec![(20, 32), (64, 64), (84, 143)];
        let expected = AtomContainer::from_ranges(ranges, 1028);

        assert_eq!(container, expected);
    }
}
