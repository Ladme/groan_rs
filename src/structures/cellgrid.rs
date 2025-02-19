// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Implementation of CellGrid for speeding up distance calculations.

// TODO! Write documentation for methods.

use ndarray::Array3;

use crate::{
    errors::{AtomError, CellGridError, PositionError, SimBoxError},
    prelude::Vector3D,
    system::System,
};

use super::{atom::Atom, simbox::SimBox};

/// Only supports orthogonal simulation boxes.
#[derive(Debug, Clone)]
pub struct CellGrid {
    grid: Array3<Vec<usize>>,
    cell_dim: Vector3D,
}

impl CellGrid {
    /// If `cell_size` is larger than any dimension of the simulation box, the smallest dimension of the box is used instead.
    pub fn new(system: &System, group: &str, cell_size: f32) -> Result<CellGrid, CellGridError> {
        let simbox = system
            .get_box()
            .ok_or_else(|| CellGridError::SimBoxError(SimBoxError::DoesNotExist))?;

        if simbox.is_zero() {
            return Err(CellGridError::SimBoxError(SimBoxError::AllDimensionsZero));
        }

        if !simbox.is_orthogonal() {
            return Err(CellGridError::SimBoxError(SimBoxError::NotOrthogonal));
        }

        let xcells = Self::n_cells(simbox.x, cell_size);
        let ycells = Self::n_cells(simbox.y, cell_size);
        let zcells = Self::n_cells(simbox.z, cell_size);

        let cell_dim = Vector3D::new(
            simbox.x / xcells as f32,
            simbox.y / ycells as f32,
            simbox.z / zcells as f32,
        );

        let mut grid: Array3<Vec<usize>> = Array3::from_elem((xcells, ycells, zcells), Vec::new());

        // assign atoms to the grid
        for atom in system
            .group_iter(group)
            .map_err(CellGridError::GroupError)?
        {
            let index = Self::atom2cell(atom, &cell_dim, simbox)?;

            match grid.get_mut(index) {
                Some(x) => {
                    x.push(atom.get_index());
                    //assignment.insert(atom.get_index(), index);
                }
                None => panic!(
                    "FATAL GROAN ERROR | CellGrid::new | Bin index `{:?}` is out of range.",
                    index
                ),
            }
        }

        Ok(CellGrid { grid, cell_dim })
    }

    /// Iterate through all atoms that are located inside the same cell as the `reference` and in the neighbouring cells.
    /// The `reference` atom IS ALSO included in the iterator. The order of iteration is undefined.
    /// The `atoms` must be all the atoms of the system.
    pub fn iter_neighbours<'a>(
        &self,
        reference: usize,
        atoms: &'a [Atom],
        simbox: &SimBox,
    ) -> Result<impl Iterator<Item = &'a Atom> + use<'a, '_>, CellGridError> {
        let (xcells, ycells, zcells) = self.grid.dim();

        let atom = atoms
            .get(reference)
            .ok_or_else(|| CellGridError::AtomError(AtomError::OutOfRange(reference)))?;
        let [x, y, z] = Self::atom2cell(atom, &self.cell_dim, simbox)?;

        Ok([-1isize, 0, 1]
            .iter()
            .flat_map(move |&dz| {
                [-1isize, 0, 1].iter().flat_map(move |&dy| {
                    [-1isize, 0, 1].iter().flat_map(move |&dx| {
                        let nx = ((x as isize + dx).rem_euclid(xcells as isize)) as usize;
                        let ny = ((y as isize + dy).rem_euclid(ycells as isize)) as usize;
                        let nz = ((z as isize + dz).rem_euclid(zcells as isize)) as usize;
                        self.grid
                            .get([nx, ny, nz])
                            .expect("FATAL GROAN ERROR | CellGrid::iter_neighbours | Nonexistent cell visited.")
                            .iter()
                    })
                })
            })
            .map(|index| atoms
                .get(*index)
                .expect("FATAL GROAN ERROR | CellGrid::iter_neighbours | Atom should be present in the vector of atoms.")
            )
        )
    }

    /// Calculate the number of bins along a dimension of the simulation box.
    ///
    /// ## Panics
    /// - Panics if `box_len` or `bin_size` are not positive.
    #[inline]
    fn n_cells(box_len: f32, cell_size: f32) -> usize {
        assert!(
            box_len > 0.0,
            "FATAL GROAN ERROR | CellGrid::n_bins | `box_len` is not positive."
        );
        assert!(
            cell_size > 0.0,
            "FATAL GROAN ERROR | CellGrid::n_bins | `cell_size` is not positive."
        );

        (box_len / cell_size).floor().max(1.0) as usize
    }

    /// Get the cell of the grid this atom is located in.
    fn atom2cell(
        atom: &Atom,
        cell_dim: &Vector3D,
        simbox: &SimBox,
    ) -> Result<[usize; 3], CellGridError> {
        let mut pos = atom
            .get_position()
            .ok_or_else(|| {
                CellGridError::AtomError(AtomError::InvalidPosition(PositionError::NoPosition(
                    atom.get_index(),
                )))
            })?
            .clone();
        // wrap the atom into the simulation box
        pos.wrap(simbox);

        Ok(Self::pos2index(&pos, cell_dim))
    }

    /// Convert a position to index of the cell of the grid.
    #[inline(always)]
    fn pos2index(pos: &Vector3D, cell_dim: &Vector3D) -> [usize; 3] {
        [
            (pos.x / cell_dim.x).floor() as usize,
            (pos.y / cell_dim.y).floor() as usize,
            (pos.z / cell_dim.z).floor() as usize,
        ]
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::Dimension;

    use super::*;

    // TODO: tests no simple, model systems

    #[test]
    fn test_real() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let distances = system
            .group_all_distances("Protein", "Protein", Dimension::XYZ)
            .unwrap();

        let nearby: Vec<Vec<usize>> = distances
            .outer_iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter_map(
                        |(col_idx, &value)| {
                            if value < 0.5 {
                                Some(col_idx)
                            } else {
                                None
                            }
                        },
                    )
                    .collect()
            })
            .collect();

        println!("{:?}", nearby);

        let cellgrid = CellGrid::new(&system, "Protein", 0.5).unwrap();

        let mut nearby_cellgrid: Vec<Vec<usize>> = Vec::new();
        for atom in system.group_iter("Protein").unwrap() {
            nearby_cellgrid.push(
                cellgrid
                    .iter_neighbours(
                        atom.get_index(),
                        system.get_atoms(),
                        system.get_box().unwrap(),
                    )
                    .unwrap()
                    .filter_map(|a| {
                        let a1_index = atom.get_index();
                        let a2_index = a.get_index();

                        if system
                            .atoms_distance(a1_index, a2_index, Dimension::XYZ)
                            .unwrap()
                            < 0.5
                        {
                            Some(a2_index)
                        } else {
                            None
                        }
                    })
                    .collect::<Vec<usize>>(),
            )
        }

        println!("{:?}", nearby_cellgrid);
    }
}
