// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Implementation of CellGrid for speeding up distance calculations.

// TODO! Write documentation for methods.

use std::ops::Range;

use itertools::iproduct;
use ndarray::Array3;

use crate::{
    errors::{AtomError, CellGridError, PositionError, SimBoxError},
    prelude::Vector3D,
    system::System,
};

use super::{atom::Atom, iterators::UnorderedAtomIterator, simbox::SimBox};

/// Only supports orthogonal simulation boxes.
#[derive(Debug, Clone)]
pub struct CellGrid<'a> {
    /// Grid of cells storing atom indices.
    grid: Array3<Vec<usize>>,
    /// Dimensions of each cell of the grid.
    cell_size: Vector3D,
    /// List of all atoms in the system.
    atoms: &'a [Atom],
    /// Simulation box.
    simbox: &'a SimBox,
}

/// Specifies the range of cells, relative to the reference cell, which should
/// be visited in each dimension when searching for neighbours in a CellGrid.
/// When this structure is used in `neighbors_iter` method, it guarantess that no
/// cell is visited multiple times.
#[derive(Debug, Clone)]
pub struct CellNeighbors {
    /// Range of neighboring cells along the x-axis.
    x: Range<isize>,
    /// Range of neighboring cells along the y-axis.
    y: Range<isize>,
    /// Range of neighboring cells along the z-axis.
    z: Range<isize>,
}

impl Default for CellNeighbors {
    /// Construct a default `CellNeighbors` structure.
    /// Use only if your cutoff is isotropic and the size of
    /// the cells of your grid is equal to or larger than the cutoff.
    fn default() -> Self {
        CellNeighbors {
            x: -1..2,
            y: -1..2,
            z: -1..2,
        }
    }
}

impl CellNeighbors {
    /// Set the range of cells, relative to the reference cell,
    /// which should be visited in each dimension when searching for
    /// neighbours in a CellGrid.
    ///
    /// ## Examples
    /// The examples assume that the size of each cell in your grid is 1×1×1 nm.
    /// - If your cutoff is 0.8 nm, use a range of `-1..=1` for each dimension.
    ///   (In this scenario, you can also simply call `CellNeighbors::default`.)
    /// - If your cutoff is 2.9 nm, use a range of `-3..=3` for each dimension.
    /// - If your cutoff is 0.9 nm in the xy-plane and 1.3 nm along the z-axis,
    ///   use a range of `-1..=1` for x- and y-dimensions and `-2..=2` for the z-dimension.
    /// - If your cutoff is 1.3 nm in the z-dimension and there is no cutoff for the
    ///   x- and y-dimensions, use a range of `-2..=2` for the z-dimension and
    ///   `isize::MIN..=isize::MAX` for the x- and y-dimensions.
    ///   *The range will be automatically "trimmed" to include all cells along these dimensions exactly once.*
    /// - If your cutoff is 1.8 nm, but only in the positive direction of each dimension
    ///   (i.e., you are not interested what is "behind" and "below" the reference point),
    ///   use a range of `0..=2` for each dimension.
    pub fn new(x: Range<isize>, y: Range<isize>, z: Range<isize>) -> Self {
        CellNeighbors { x, y, z }
    }

    /// Trim the ranges so that each cell is visited at most once.
    fn trim(&mut self, xcells: usize, ycells: usize, zcells: usize) {
        if self.x.len() > xcells {
            self.x = 0..xcells as isize;
        }

        if self.y.len() > ycells {
            self.y = 0..ycells as isize;
        }

        if self.z.len() > zcells {
            self.z = 0..zcells as isize;
        }
    }
}

impl<'a> CellGrid<'a> {
    /// If `cell_size` is larger than any dimension of the simulation box, the smallest dimension of the box is used instead.
    /// TODO: error for negative cell size
    pub fn new(
        system: &'a System,
        group: &str,
        cell_size: f32,
    ) -> Result<CellGrid<'a>, CellGridError> {
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

        let cell_size = Vector3D::new(
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
            let index = Self::atom2cell(atom, &cell_size, simbox)?;

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

        Ok(CellGrid {
            grid,
            cell_size,
            atoms: system.get_atoms(),
            simbox,
        })
    }

    /// Returns an iterator over all atoms that have been assigned into the `CellGrid`
    /// and that are located in the same cell as the `reference` point and its neighboring cells.
    ///
    /// ## Params
    /// - `reference` - coordinates of the reference point
    /// - `ranges` - specifies what cells should be considered neighbors,
    ///    see [`CellNeighbors`] for more information
    ///
    /// ## Notes
    /// - Assumes periodic boundary conditions in all dimensions.
    /// - The order in which the neighboring atoms are visited is **undefined**,
    ///   but each atom will be visited at most once.
    /// - The `reference` point does not need to be wrapped into the simulation box.
    ///   It is wrapped automatically inside the method.
    /// - Note that **no** optimization in the form of removing distant diagonal cells is performed.
    ///
    pub fn neighbors_iter(
        &'a self,
        mut reference: Vector3D,
        mut ranges: CellNeighbors,
    ) -> Result<UnorderedAtomIterator<'a, impl Iterator<Item = usize> + 'a>, CellGridError> {
        reference.wrap(&self.simbox);
        let [x, y, z] = Self::pos2index(&reference, &self.cell_size);
        let (xcells, ycells, zcells) = self.grid.dim();

        // make sure that no cell is visited multiple times
        ranges.trim(xcells, ycells, zcells);

        let iterator = iproduct!(ranges.z, ranges.y, ranges.x).flat_map(move |(dz, dy, dx)| {
            let nx = ((x as isize + dx).rem_euclid(xcells as isize)) as usize;
            let ny = ((y as isize + dy).rem_euclid(ycells as isize)) as usize;
            let nz = ((z as isize + dz).rem_euclid(zcells as isize)) as usize;

            self.grid
                .get([nx, ny, nz])
                .expect("FATAL GROAN ERROR | CellGrid::iter_neighbours | Nonexistent cell visited.")
                .iter()
                .cloned()
        });

        Ok(UnorderedAtomIterator::new(
            &self.atoms,
            iterator,
            Some(&self.simbox),
        ))
    }

    /// Calculate the number of cells along a dimension of the simulation box.
    ///
    /// ## Panics
    /// - Panics if `box_len` or `cell_size` are not positive.
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

    /// Get the index of the cell this atom is located in.
    fn atom2cell(
        atom: &Atom,
        cell_size: &Vector3D,
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

        Ok(Self::pos2index(&pos, cell_size))
    }

    /// Convert a position to index of the cell of the grid.
    #[inline(always)]
    fn pos2index(pos: &Vector3D, cell_size: &Vector3D) -> [usize; 3] {
        [
            (pos.x / cell_size.x).floor() as usize,
            (pos.y / cell_size.y).floor() as usize,
            (pos.z / cell_size.z).floor() as usize,
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

        println!("{:?}", cellgrid.cell_size);
        println!("{:?}", cellgrid.grid.dim());
        println!("{:?}", cellgrid.grid);

        let mut nearby_cellgrid: Vec<Vec<usize>> = Vec::new();
        for atom in system.group_iter("Protein").unwrap() {
            nearby_cellgrid.push(
                cellgrid
                    .neighbors_iter(
                        atom.get_position().unwrap().clone(),
                        CellNeighbors::default(),
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
