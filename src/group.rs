// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of the Group structure and associated methods.

/// Group of atoms in target system.
#[derive(Debug, Clone)]
pub struct Group {
    pub atom_ranges: Vec<(u64, u64)>,
}

impl Group {

    /// Create a new valid Group structure from atom ranges.
    /// ## Parameters
    /// Expects a vector of atom ranges (start, end) and the total number of atoms in the system.
    /// The individual "atom ranges" in the atom ranges can be overlapping. 
    /// In the final Group structure, they will be merged together.
    pub fn from_ranges(atom_ranges: Vec<(u64, u64)>, n_atoms: u64) -> Self {

        let merged_ranges = Group::fix_atom_ranges(atom_ranges, n_atoms);
        Group { atom_ranges: merged_ranges }
    }

    /// Create a new valid Group structure from atom indices.
    /// ## Parameters
    /// Expects a vector of atom indices and the total number of atoms in the system.
    /// There can be duplicate atoms in the "atom indices". In the final Group structure, they will be removed.
    pub fn from_indices(atom_indices: Vec<u64>, n_atoms: u64) -> Self {

        let ranges = Group::make_atom_ranges(atom_indices, n_atoms);
        Group { atom_ranges: ranges }
    }


    /// Fix overlaps in atom ranges. 
    /// Makes sure that the atom ranges are valid and not overflowing the number of atoms in the system.
    fn fix_atom_ranges(mut atom_ranges: Vec<(u64, u64)>, n_atoms: u64) -> Vec<(u64, u64)> {
        if atom_ranges.is_empty() {
            return atom_ranges;
        }

        // sort the atom ranges in ascending order
        atom_ranges.sort_unstable(); 
        
        let mut merged_indices = Vec::new();
        let mut current_start = std::u64::MAX;
        let mut current_end = 0u64;
        
        for (start, end) in &atom_ranges {
            // start must be smaller than n_atoms and not larger than end
            if *start >= n_atoms || *start > *end  {
                continue;
            }

            // current range does not overlap with the previous one nor is adjacent to it
            if *start > current_end + 1 || (current_end == 0u64 && current_start != 0u64) {
                if current_start != std::u64::MAX {
                    merged_indices.push((current_start, current_end.min(n_atoms - 1)));
                }
                current_start = *start;
                current_end = *end;
            // current range overlaps
            } else if *end > current_end {
                current_end = *end;
            }
        }
        
        if current_start != std::u64::MAX {
            // add the last merged range to the result if it exists
            merged_indices.push((current_start, current_end.min(n_atoms - 1)));
        }
        
        merged_indices
    }

    /// Create valid atom ranges from atom indices. 
    /// Makes sure that the atom ranges are valid and not overflowing the number of atoms in the system.
    fn make_atom_ranges(mut atom_indices: Vec<u64>, n_atoms: u64) -> Vec<(u64, u64)> {
        let mut atom_ranges = Vec::new();
        
        if atom_indices.is_empty() {
            return atom_ranges;
        }

        // sort the atom indices in ascending order
        atom_indices.sort_unstable();

        let mut start = atom_indices[0];
        let mut end = atom_indices[0];

        for i in 1..atom_indices.len() {
            // the range can't exceed the number of atoms
            if atom_indices[i] >= n_atoms {
                end = n_atoms - 1;
                break;
            }

            // remove duplicate atoms
            if atom_indices[i] == end {
                continue;
            }

            if atom_indices[i] == end + 1 {
                // extend the current range
                end = atom_indices[i];
            } else {
                // add the completed range to the atom_indices vector and start a new range
                atom_ranges.push((start, end));
                start = atom_indices[i];
                end = atom_indices[i];
            }
        }

        // Add the last completed range to the atom_indices vector
        atom_ranges.push((start, end));

        atom_ranges
    }


    /// Get the number of atoms in the group.
    pub fn get_n_atoms(&self) -> u64 {
        self.atom_ranges
            .iter()
            .fold(0, |acc, (start, end)| acc + (end - start + 1))
    }


}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_group_simple() {

        let atoms = vec![(20, 32)];

        let group = Group::from_ranges(atoms, 33);

        assert_eq!(group.atom_ranges[0], (20, 32));
    }

    #[test]
    fn test_new_group_ranges_empty() {

        let atoms = vec![];
        let group = Group::from_ranges(atoms, 1028);
        assert!(group.atom_ranges.is_empty());
    }

    #[test]
    fn test_new_group_multiple_nooverlap() {

        let atoms = vec![(20, 32), (64, 64), (84, 143)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (20, 32));
        assert_eq!(group.atom_ranges[1], (64, 64));
        assert_eq!(group.atom_ranges[2], (84, 143));
    }

    #[test]
    fn test_new_group_adjacent() {
        
        let atoms = vec![(20, 32), (33, 42)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (20, 42));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn test_new_group_overlap_simple() {
        
        let atoms = vec![(20, 32), (24, 42)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (20, 42));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn test_new_group_overlap_full() {
        
        let atoms = vec![(20, 32), (28, 30)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (20, 32));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn test_new_group_overlap_full_flipped() {
        
        let atoms = vec![(28, 30), (20, 32)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (20, 32));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn test_new_group_barely_overlaps() {

        let atoms = vec![(20, 32), (32, 42)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (20, 42));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn test_new_group_overlap_complex() {
        
        let atoms = vec![(64, 128), (5, 32), (1, 25), (129, 133), (133, 200), (35, 78), (10, 15)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (1, 32));
        assert_eq!(group.atom_ranges[1], (35, 200));
        assert_eq!(group.atom_ranges.len(), 2);
    }

    #[test]
    fn test_new_group_overlaps_with_zeros() {
        
        let atoms = vec![(1, 25), (0, 1), (0, 0), (0, 34)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (0, 34));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn test_new_group_start_larger_than_end() {
        
        let atoms = vec![(32, 25), (14, 17)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (14, 17));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn test_new_group_end_larger_than_natoms() {
        
        let atoms = vec![(543, 1020), (1000, 1432)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (543, 1027));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn test_new_group_start_larger_than_natoms() {
        
        let atoms = vec![(543, 1020), (1043, 1432)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (543, 1020));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn test_new_group_end_larger_than_natoms_nooverlap() {
        
        let atoms = vec![(0, 43), (1006, 1432)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (0, 43));
        assert_eq!(group.atom_ranges[1], (1006, 1027));
        assert_eq!(group.atom_ranges.len(), 2);
    }

    #[test]
    fn test_new_group_single_atoms() {
        
        let atoms = vec![(5, 5), (4, 4), (11, 11), (12, 12), (0, 0), (1, 1), (6, 6), (2, 2), (3, 3), (7, 7), (13, 13), (10, 10), (8, 8), (9, 9)];

        let group = Group::from_ranges(atoms, 1028);

        assert_eq!(group.atom_ranges[0], (0, 13));
        assert_eq!(group.atom_ranges.len(), 1);
    }

    #[test]
    fn test_new_group_from_indices_basic() {
        let atom_indices = vec![6, 2, 13, 1, 10, 8, 3, 12, 7, 14, 15];
        let group = Group::from_indices(atom_indices, 16);

        assert_eq!(group.atom_ranges[0], (1, 3));
        assert_eq!(group.atom_ranges[1], (6, 8));
        assert_eq!(group.atom_ranges[2], (10, 10));
        assert_eq!(group.atom_ranges[3], (12, 15));
        assert_eq!(group.atom_ranges.len(), 4);
    }

    #[test]
    fn test_new_group_from_indices_empty() {

        let atoms = vec![];
        let group = Group::from_indices(atoms, 1028);
        assert!(group.atom_ranges.is_empty());
    }

    #[test]
    fn test_new_group_from_indices_duplicate() {
        let atom_indices = vec![1, 6, 3, 2, 13, 1, 10, 8, 3, 12, 7, 14, 15, 10];
        let group = Group::from_indices(atom_indices, 20);

        assert_eq!(group.atom_ranges[0], (1, 3));
        assert_eq!(group.atom_ranges[1], (6, 8));
        assert_eq!(group.atom_ranges[2], (10, 10));
        assert_eq!(group.atom_ranges[3], (12, 15));
        assert_eq!(group.atom_ranges.len(), 4);
    }

    #[test]
    fn test_new_group_from_indices_larger_than_natoms() {
        let atom_indices = vec![1, 6, 3, 2, 13, 1, 10, 8, 3, 12, 7, 14, 15, 10];
        let group = Group::from_indices(atom_indices, 15);

        assert_eq!(group.atom_ranges[0], (1, 3));
        assert_eq!(group.atom_ranges[1], (6, 8));
        assert_eq!(group.atom_ranges[2], (10, 10));
        assert_eq!(group.atom_ranges[3], (12, 14));
        assert_eq!(group.atom_ranges.len(), 4);
    }

    #[test]
    fn test_get_n_atoms_empty() {

        let atom_ranges = vec![];
        let group = Group::from_ranges(atom_ranges, 100);

        assert_eq!(group.get_n_atoms(), 0);
    }

    #[test]
    fn test_get_n_atoms_basic() {

        let atom_ranges = vec![(64, 128), (5, 32), (1, 25), (129, 133), (133, 200), (35, 78), (10, 15)];
        let group = Group::from_ranges(atom_ranges, 1028);

        assert_eq!(group.get_n_atoms(), 198);
    }


}