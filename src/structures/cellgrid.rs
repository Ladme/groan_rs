// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Implementation of CellGrid for speeding up distance calculations.

use std::ops::{Range, RangeFull, RangeInclusive};

use itertools::iproduct;
use ndarray::Array3;

use crate::{
    errors::{AtomError, CellGridError, PositionError, SimBoxError},
    prelude::Vector3D,
    system::System,
};

use super::{atom::Atom, iterators::UnorderedAtomIterator, simbox::SimBox};

/// A structure for efficient searches within a cutoff.
///
/// Also commonly known as **cell lists**. See [Wikipedia](https://en.wikipedia.org/wiki/Cell_lists) for more information.
///
/// **`CellGrid` can only be used in systems with orthogonal simulation boxes!**
///
/// Theoretically, `CellGrid` reduces the time complexity of pairwise calculations from O(n^2) to O(n),
/// where **n** is the number of atoms. In practice, the usefulness of this structure increases with system size.
///
/// Construct a `CellGrid` using [`CellGrid::new`].
///
/// ## Example usage
///
/// Suppose you want to calculate the local membrane center for each lipid molecule in your system.
/// One way to achieve this is by performing geometry filtering using a cylinder of a specified radius
/// for each lipid molecule:
///
/// ```no_run
/// # use groan_rs::prelude::*;
/// #
/// // Calculate the local membrane center for each lipid molecule in a molecular system.
/// fn calc_local_membrane_centers() -> Result<Vec<Vector3D>, Box<dyn std::error::Error + Send + Sync>>
/// {
///     // read the system structure from a GRO file
///     let mut system = System::from_file("system.gro")?;
///     // create a group for membrane atoms...
///     system.group_create("Membrane", "@membrane")?;
///     // ...and a group for lipid 'heads'; assuming there is only one such atom per lipid molecule
///     system.group_create("Heads", "Membrane and name P")?;
///
///     let mut local_centers = Vec::new();
///     // iterate through lipid heads
///     for head in system.group_iter("Heads")? {
///         // construct a cylinder for selecting membrane atoms in the vicinity of the head
///         let cylinder = Cylinder::new(
///             // the center of the cylinder base is the lipid head
///             head.get_position().unwrap().clone(),
///             // the cylinder has a radius of 2 nm
///             2.0,
///             // the cylinder has infinite height
///             f32::INFINITY,
///             // the cylinder is oriented along the z-axis since the membrane is in the xy-plane
///             Dimension::Z,
///         );
///
///         let center = system
///             // iterate through ALL membrane atoms
///             .group_iter("Membrane")?
///             // select only the atoms inside the cylinder
///             .filter_geometry(cylinder)
///             // compute the geometric center of the selected atoms
///             .get_center()?;
///
///         // store the geometric center for this lipid
///         local_centers.push(center);
///     }
///
///     Ok(local_centers)
/// }
/// ```
///
/// However, the above approach can be computationally expensive, especially in large systems.
/// The time complexity is quadratic (`O(n * m)`), where `n` is the number of membrane atoms
/// and `m` is the number of lipid molecules.
///
/// In large systems, it is much more efficient to construct a `CellGrid`
/// and use it to efficiently select nearby atoms:
///
/// ```no_run
/// # use groan_rs::prelude::*;
/// #
/// // Calculate the local membrane center for each lipid molecule in a molecular system using a CellGrid.
/// fn calc_local_membrane_centers_cellgrid(
/// ) -> Result<Vec<Vector3D>, Box<dyn std::error::Error + Send + Sync>> {
///     // read the system structure and construct the necessary groups
///     let mut system = System::from_file("system.gro")?;
///     system.group_create("Membrane", "@membrane")?;
///     system.group_create("Heads", "Membrane and name P")?;
///
///     // construct the cell grid
///     // only the atoms of the `Membrane` group will be included
///     // the size of each cell inside the grid is 2 nm, which corresponds to the cylinder's radius
///     let cellgrid = CellGrid::new(&system, "Membrane", 2.0)?;
///
///     let mut local_centers = Vec::new();
///     // iterate through lipid heads
///     for head in system.group_iter("Heads")? {
///         // construct a cylinder for geometry selection as before
///         let cylinder = Cylinder::new(
///             head.get_position().unwrap().clone(),
///             2.0,
///             f32::INFINITY,
///             Dimension::Z,
///         );
///
///         // iterate only through atoms in neighboring cells
///         let center = cellgrid
///             .neighbors_iter(
///                 // reference point: position of the lipid head
///                 head.get_position().unwrap().clone(),
///                 // define neighboring cells
///                 // selects the cell containing the reference point and all adjacent cells in the xy-plane
///                 // additionally, it considers all cells above and below along the z-axis, since the
///                 // cylinder has an infinite height in the z-direction
///                 CellNeighbors::new(-1..=1, -1..=1, ..),
///             )
///             // select only the neighboring atoms that are inside the cylinder
///             .filter_geometry(cylinder)
///             // compute the geometric center of these atoms
///             .get_center()?;
///
///         local_centers.push(center);
///     }
///
///     Ok(local_centers)
/// }
/// ```
///
/// **Important note:**  
/// When comparing local membrane centers calculated using the "naive" approach and the `CellGrid` method,
/// you may notice that the results are not *exactly* identical. This is because [`System::group_iter`]
/// and [`CellGrid::neighbors_iter`] do not always return atoms in the same order. In fact, the order
/// of atoms returned by `CellGrid::neighbors_iter` is undefined.
/// Since the geometric center is computed by averaging atom positions, and floating-point summation is **not**
/// a commutative operation (due to rounding errors), small differences may arise.
///
/// However, there is no definitive way to determine which approach is more "correct," so this difference
/// should **not** be considered a drawback of the `CellGrid` method.
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
///
/// See [`CellNeighbors::new`] for more information.
#[derive(Debug, Clone)]
pub struct CellNeighbors {
    /// Range of neighboring cells along the x-axis.
    x: NeighborsRange,
    /// Range of neighboring cells along the y-axis.
    y: NeighborsRange,
    /// Range of neighboring cells along the z-axis.
    z: NeighborsRange,
}

/// Helper enum for specifying the range for neighbors selection.
#[derive(Debug, Clone)]
pub enum NeighborsRange {
    Exclusive(Range<isize>),
    Inclusive(RangeInclusive<isize>),
    Full(RangeFull),
}

impl Default for NeighborsRange {
    fn default() -> Self {
        Self::Exclusive(-1..2)
    }
}

impl From<Range<isize>> for NeighborsRange {
    fn from(value: Range<isize>) -> Self {
        Self::Exclusive(value)
    }
}

impl From<RangeInclusive<isize>> for NeighborsRange {
    fn from(value: RangeInclusive<isize>) -> Self {
        Self::Inclusive(value)
    }
}

impl From<RangeFull> for NeighborsRange {
    fn from(value: RangeFull) -> Self {
        Self::Full(value)
    }
}

impl NeighborsRange {
    /// Trim the provided `NeighborsRange` so that it does not iterate
    /// through the same cell multiple times due to PBC wrapping,
    /// and convert it to exclusive range for iteration.
    fn convert(self, n: usize) -> Range<isize> {
        match self {
            Self::Exclusive(x) => {
                if x.len() > n {
                    0..n as isize
                } else {
                    x
                }
            }
            Self::Inclusive(x) => {
                if x.end() - x.start() >= n as isize {
                    0..n as isize
                } else {
                    (*x.start())..(*(x.end()) + 1)
                }
            }
            Self::Full(_) => 0..n as isize,
        }
    }
}

impl Default for CellNeighbors {
    /// Construct a default `CellNeighbors` structure.
    /// Use only if your cutoff is isotropic and the size of
    /// the cells of your grid is equal to or larger than the cutoff.
    fn default() -> Self {
        CellNeighbors {
            x: NeighborsRange::default(),
            y: NeighborsRange::default(),
            z: NeighborsRange::default(),
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
    ///   `FullRange` (or `..`) for the x- and y-dimensions.
    /// - If your cutoff is 1.8 nm, but only in the positive direction of each dimension
    ///   (i.e., you are not interested what is "behind" and "below" the reference point),
    ///   use a range of `0..=2` for each dimension.
    pub fn new(
        x: impl Into<NeighborsRange>,
        y: impl Into<NeighborsRange>,
        z: impl Into<NeighborsRange>,
    ) -> Self {
        CellNeighbors {
            x: x.into(),
            y: y.into(),
            z: z.into(),
        }
    }

    /// Convert `NeighborsRanges` into iterable exclusive ranges.
    /// Makes sure that no cell is visited multiple times.
    #[inline(always)]
    fn convert(
        self,
        xcells: usize,
        ycells: usize,
        zcells: usize,
    ) -> (Range<isize>, Range<isize>, Range<isize>) {
        (
            self.x.convert(xcells),
            self.y.convert(ycells),
            self.z.convert(zcells),
        )
    }
}

impl<'a> CellGrid<'a> {
    /// Create a new [`CellGrid`]. The time complexity of this operation is O(N), where N is the number
    /// of atoms in the specified group.
    ///
    /// ## Notes
    /// - Only atoms belonging to the specified `group` will be assigned to the `CellGrid`.
    /// - The simulation box must be defined and orthogonal; otherwise, an error will be raised.
    /// - If `cell_size` is larger than any dimension of the simulation box, the smallest dimension of the box will be used instead.
    /// - The `cell_size` must be positive; otherwise, an error will be raised.
    pub fn new(
        system: &'a System,
        group: &str,
        cell_size: f32,
    ) -> Result<CellGrid<'a>, CellGridError> {
        if cell_size <= 0.0 {
            return Err(CellGridError::InvalidCellSize(cell_size.to_string()));
        }

        let simbox = Self::check_box(system)?;

        let xcells = Self::n_cells(simbox.x, cell_size);
        let ycells = Self::n_cells(simbox.y, cell_size);
        let zcells = Self::n_cells(simbox.z, cell_size);
        let cells = [xcells, ycells, zcells];

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
            let index = Self::atom2cell(atom, &cell_size, &cells, simbox)?;

            match grid.get_mut(index) {
                Some(x) => {
                    x.push(atom.get_index());
                }
                None => panic!(
                    "FATAL GROAN ERROR | CellGrid::new | Cell index `{:?}` is out of range.",
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
    /// and that are located in cells neighboring the `reference` point.
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
    ///
    pub fn neighbors_iter(
        &'a self,
        mut reference: Vector3D,
        ranges: CellNeighbors,
    ) -> UnorderedAtomIterator<'a, impl Iterator<Item = usize> + 'a> {
        reference.wrap(&self.simbox);
        let cells = self.grid.dim();
        let [x, y, z] = Self::pos2index(&reference, &self.cell_size, &[cells.0, cells.1, cells.2]);
        let (xcells, ycells, zcells) = self.grid.dim();

        // make sure that no cell is visited multiple times
        let (xrange, yrange, zrange) = ranges.convert(xcells, ycells, zcells);

        let iterator = iproduct!(xrange, yrange, zrange).flat_map(move |(dx, dy, dz)| {
            let nx = ((x as isize + dx).rem_euclid(xcells as isize)) as usize;
            let ny = ((y as isize + dy).rem_euclid(ycells as isize)) as usize;
            let nz = ((z as isize + dz).rem_euclid(zcells as isize)) as usize;

            self.grid
                .get([nx, ny, nz])
                .expect("FATAL GROAN ERROR | CellGrid::neighbors_iter | Nonexistent cell visited.")
                .iter()
                .cloned()
        });

        UnorderedAtomIterator::new(&self.atoms, iterator, Some(&self.simbox))
    }

    /// Check that the simulation box exists, is non-zero and orthogonal.
    /// Returns the simulation box, if it is valid. Otherwise returns an error.
    fn check_box(system: &'a System) -> Result<&'a SimBox, CellGridError> {
        let simbox = system
            .get_box()
            .ok_or_else(|| CellGridError::SimBoxError(SimBoxError::DoesNotExist))?;

        if simbox.is_zero() {
            return Err(CellGridError::SimBoxError(SimBoxError::AllDimensionsZero));
        }

        if !simbox.is_orthogonal() {
            return Err(CellGridError::SimBoxError(SimBoxError::NotOrthogonal));
        }

        Ok(simbox)
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
        ncells: &[usize; 3],
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

        Ok(Self::pos2index(&pos, cell_size, ncells))
    }

    /// Convert a position to index of the cell of the grid.
    #[inline(always)]
    fn pos2index(pos: &Vector3D, cell_size: &Vector3D, ncells: &[usize; 3]) -> [usize; 3] {
        [
            ((pos.x / cell_size.x).floor() as usize).rem_euclid(ncells[0]),
            ((pos.y / cell_size.y).floor() as usize).rem_euclid(ncells[1]),
            ((pos.z / cell_size.z).floor() as usize).rem_euclid(ncells[2]),
        ]
    }
}

#[cfg(test)]
mod tests {
    use float_cmp::assert_approx_eq;
    use itertools::izip;
    use rand::{seq::SliceRandom, Rng, SeedableRng};

    use crate::prelude::{AtomIteratorWithBox, Cylinder, Dimension, Rectangular, Shape, Sphere};

    use super::*;

    #[test]
    fn test_new_cellgrid_pass_too_large_cell() {
        let system = System::from_file("test_files/example.gro").unwrap();
        let cellgrid = CellGrid::new(&system, "all", 20.0).unwrap();
        let dim = cellgrid.grid.dim();
        assert_eq!(dim.0, 1);
        assert_eq!(dim.1, 1);
        assert_eq!(dim.2, 1);

        let simbox = system.get_box().unwrap();
        assert_approx_eq!(f32, cellgrid.cell_size.x, simbox.x);
        assert_approx_eq!(f32, cellgrid.cell_size.y, simbox.y);
        assert_approx_eq!(f32, cellgrid.cell_size.z, simbox.z);
    }

    #[test]
    fn test_new_cellgrid_fail_zero_cell() {
        let system = System::from_file("test_files/example.gro").unwrap();
        match CellGrid::new(&system, "all", 0.0) {
            Ok(_) => panic!("Function should have failed."),
            Err(CellGridError::InvalidCellSize(_)) => (),
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_new_cellgrid_fail_negative_cell() {
        let system = System::from_file("test_files/example.gro").unwrap();
        match CellGrid::new(&system, "all", -1.5) {
            Ok(_) => panic!("Function should have failed."),
            Err(CellGridError::InvalidCellSize(_)) => (),
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_new_cellgrid_fail_no_simbox() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.reset_box();
        match CellGrid::new(&system, "all", 1.0) {
            Ok(_) => panic!("Function should have failed."),
            Err(CellGridError::SimBoxError(_)) => (),
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_new_cellgrid_fail_simbox_not_orthogonal() {
        let system = System::from_file("test_files/dodecahedron.gro").unwrap();
        match CellGrid::new(&system, "all", 1.0) {
            Ok(_) => panic!("Function should have failed."),
            Err(CellGridError::SimBoxError(_)) => (),
            Err(e) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn test_new_cellgrid_at_edges() {
        let atoms = iproduct!(1..11, 1..10, 1..14)
            .map(|(x, y, z)| {
                Atom::new(1, "LYS", 1, "BB")
                    .with_position(Vector3D::new(x as f32, y as f32, z as f32))
            })
            .collect::<Vec<Atom>>();

        let system = System::new(
            "Artificial system for testing cell grid",
            atoms,
            Some([10.0, 9.0, 13.0].into()),
        );

        let cellgrid = CellGrid::new(&system, "all", 1.0).unwrap();
        for cell in cellgrid.grid.iter() {
            assert_eq!(cell.len(), 1);
        }
    }

    #[test]
    fn test_artificial_neighbors_iter() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let atoms = iproduct!(0..10, 0..9, 0..13)
            .map(|(x, y, z)| {
                let (dx, dy, dz): (f32, f32, f32) = (rng.gen(), rng.gen(), rng.gen());
                Atom::new(1, "LYS", 1, "BB").with_position(Vector3D::new(
                    x as f32 + dx,
                    y as f32 + dy,
                    z as f32 + dz,
                ))
            })
            .collect::<Vec<Atom>>();

        let system = System::new(
            "Artificial system for testing cell grid",
            atoms,
            Some([10.0, 9.0, 13.0].into()),
        );

        let cellgrid = CellGrid::new(&system, "all", 1.0).unwrap();
        // sanity check cell grid
        for cell in cellgrid.grid.iter() {
            assert_eq!(cell.len(), 1);
        }
        assert_approx_eq!(f32, cellgrid.cell_size.x, 1.0);
        assert_approx_eq!(f32, cellgrid.cell_size.y, 1.0);
        assert_approx_eq!(f32, cellgrid.cell_size.z, 1.0);

        for (neighbors, expected) in [
            (CellNeighbors::default(), 27),
            (
                CellNeighbors::new(-1..2, NeighborsRange::default(), -2..=2),
                45,
            ),
            (CellNeighbors::new(-2..=2, -1..=1, -1..2), 45),
            (
                CellNeighbors::new(-1..=1, -2..3, NeighborsRange::default()),
                45,
            ),
            (CellNeighbors::new(-1..=1, -1..2, RangeFull), 117),
            (CellNeighbors::new(-1..2, RangeFull, -1..=1), 81),
            (
                CellNeighbors::new(RangeFull, -1..=1, NeighborsRange::default()),
                90,
            ),
            (CellNeighbors::new(-4..=2, 0..3, -2..=2), 105),
            (CellNeighbors::new(4..=6, -3..-1, 0..=7), 48),
            (CellNeighbors::new(0..5, -4..=0, -1..2), 75),
            (CellNeighbors::new(-3..=5, 0..8, -2..10), 864),
            (CellNeighbors::new(0..11, -15..=-4, 7..23), 1170),
            (CellNeighbors::new(2..=15, -1345..2947, 5..=18), 1170),
            (CellNeighbors::new(RangeFull, RangeFull, RangeFull), 1170),
        ] {
            for (x, y, z) in iproduct!(-1..11, -1..10, -1..14) {
                let (dx, dy, dz): (f32, f32, f32) = (rng.gen(), rng.gen(), rng.gen());
                assert_eq!(
                    cellgrid
                        .neighbors_iter(
                            Vector3D::new(x as f32 + dx, y as f32 + dy, z as f32 + dz),
                            neighbors.clone()
                        )
                        .count(),
                    expected
                );
            }
        }
    }

    fn compare_nearby_atoms(lhs: &Vec<Vec<usize>>, rhs: &Vec<Vec<usize>>) {
        for (a1, a2) in lhs.iter().zip(rhs.iter()) {
            for nearby in a1.iter() {
                assert!(a2.contains(nearby), "Mismatch between nearby atoms.");
            }
        }
    }

    fn add_nearby_atoms_naive(
        system: &System,
        geometry: impl Shape,
        nearby_atoms: &mut Vec<Vec<usize>>,
    ) {
        nearby_atoms.push(
            system
                .group_iter("Group")
                .unwrap()
                .filter_geometry(geometry)
                .map(|atom| atom.get_index())
                .collect::<Vec<usize>>(),
        );
    }

    fn collect_nearby_sphere_naive(system: &System, radius: f32) -> Vec<Vec<usize>> {
        let mut nearby_atoms = Vec::new();
        for atom in system.group_iter("Group").unwrap() {
            let sphere = Sphere::new(atom.get_position().unwrap().clone(), radius);
            add_nearby_atoms_naive(&system, sphere, &mut nearby_atoms);
        }

        nearby_atoms
    }

    fn collect_nearby_cylinder_naive(
        system: &System,
        radius: f32,
        height: f32,
        dimension: Dimension,
    ) -> Vec<Vec<usize>> {
        let mut nearby_atoms = Vec::new();
        for atom in system.group_iter("Group").unwrap() {
            let cylinder = Cylinder::new(
                atom.get_position().unwrap().clone(),
                radius,
                height,
                dimension,
            );
            add_nearby_atoms_naive(&system, cylinder, &mut nearby_atoms);
        }

        nearby_atoms
    }

    fn collect_nearby_rectangular_naive(
        system: &System,
        x: f32,
        y: f32,
        z: f32,
    ) -> Vec<Vec<usize>> {
        let mut nearby_atoms = Vec::new();
        for atom in system.group_iter("Group").unwrap() {
            let rectangular = Rectangular::new(atom.get_position().unwrap().clone(), x, y, z);
            add_nearby_atoms_naive(&system, rectangular, &mut nearby_atoms);
        }

        nearby_atoms
    }

    fn add_nearby_atoms_cellgrid(
        cellgrid: &CellGrid,
        position: Vector3D,
        geometry: impl Shape,
        neighbors: CellNeighbors,
        nearby_atoms: &mut Vec<Vec<usize>>,
    ) {
        nearby_atoms.push(
            cellgrid
                .neighbors_iter(position, neighbors)
                .filter_geometry(geometry)
                .map(|atom| atom.get_index())
                .collect::<Vec<usize>>(),
        );
    }

    // this tries several different cell sizes and checks that they all return the same results
    fn collect_nearby_sphere_cellgrid(system: &System, radius: f32) -> Vec<Vec<usize>> {
        let mut last_nearby_atoms = None;
        for (modifier, range) in [(0.2, 5), (0.3, 4), (0.5, 2), (1.0, 1), (1.0, 2), (2.0, 1)] {
            let mut nearby_atoms = Vec::new();
            let cellgrid = CellGrid::new(system, "Group", radius * modifier).unwrap();

            for atom in system.group_iter("Group").unwrap() {
                let sphere = Sphere::new(atom.get_position().unwrap().clone(), radius);
                add_nearby_atoms_cellgrid(
                    &cellgrid,
                    atom.get_position().unwrap().clone(),
                    sphere,
                    CellNeighbors::new(-range..=range, -range..=range, -range..=range),
                    &mut nearby_atoms,
                );
            }

            match last_nearby_atoms {
                Some(x) => compare_nearby_atoms(&nearby_atoms, &x),
                None => (),
            }

            last_nearby_atoms = Some(nearby_atoms);
        }

        last_nearby_atoms.unwrap()
    }

    fn get_cylinder_ranges(
        radius: f32,
        height: f32,
        dimension: Dimension,
        cell_size: f32,
    ) -> CellNeighbors {
        let radius_range = (radius / cell_size).round() as isize;
        let height_range = (height / cell_size).round() as isize;

        match dimension {
            Dimension::X => CellNeighbors::new(
                0..=height_range,
                -radius_range..=radius_range,
                -radius_range..=radius_range,
            ),
            Dimension::Y => CellNeighbors::new(
                -radius_range..=radius_range,
                0..=height_range,
                -radius_range..=radius_range,
            ),
            Dimension::Z => CellNeighbors::new(
                -radius_range..=radius_range,
                -radius_range..=radius_range,
                0..=height_range,
            ),
            _ => panic!("Invalid dimension."),
        }
    }

    fn collect_nearby_cylinder_cellgrid(
        system: &System,
        radius: f32,
        height: f32,
        dimension: Dimension,
        cell_size: f32,
    ) -> Vec<Vec<usize>> {
        let mut nearby_atoms = Vec::new();
        let cellgrid = CellGrid::new(system, "Group", cell_size).unwrap();
        let ranges = get_cylinder_ranges(radius, height, dimension, cell_size);

        for atom in system.group_iter("Group").unwrap() {
            let cylinder = Cylinder::new(
                atom.get_position().unwrap().clone(),
                radius,
                height,
                dimension,
            );
            add_nearby_atoms_cellgrid(
                &cellgrid,
                atom.get_position().unwrap().clone(),
                cylinder,
                ranges.clone(),
                &mut nearby_atoms,
            );
        }

        nearby_atoms
    }

    fn get_rectangular_ranges(x: f32, y: f32, z: f32, cell_size: f32) -> CellNeighbors {
        let xcells = (x / cell_size).round() as isize;
        let ycells = (y / cell_size).round() as isize;
        let zcells = (z / cell_size).round() as isize;

        CellNeighbors::new(0..=xcells, 0..=ycells, 0..=zcells)
    }

    fn collect_nearby_rectangular_cellgrid(
        system: &System,
        x: f32,
        y: f32,
        z: f32,
        cell_size: f32,
    ) -> Vec<Vec<usize>> {
        let mut nearby_atoms = Vec::new();
        let cellgrid = CellGrid::new(system, "Group", cell_size).unwrap();
        let ranges = get_rectangular_ranges(x, y, z, cell_size);

        for atom in system.group_iter("Group").unwrap() {
            let rectangular = Rectangular::new(atom.get_position().unwrap().clone(), x, y, z);
            add_nearby_atoms_cellgrid(
                &cellgrid,
                atom.get_position().unwrap().clone(),
                rectangular,
                ranges.clone(),
                &mut nearby_atoms,
            );
        }

        nearby_atoms
    }

    #[test]
    fn test_real_geometry_sphere() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let cutoffs = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 15.0];
        for cutoff in cutoffs {
            let selected_atoms = system
                .get_atoms()
                .choose_multiple(&mut rng, 500)
                .map(|atom| atom.get_index())
                .collect::<Vec<usize>>();

            // ignore the error
            let _ = system.group_create_from_indices("Group", selected_atoms);
            assert_eq!(system.group_get_n_atoms("Group").unwrap(), 500);

            let nearby_naive = collect_nearby_sphere_naive(&system, cutoff);
            let nearby_cellgrid = collect_nearby_sphere_cellgrid(&system, cutoff);

            compare_nearby_atoms(&nearby_naive, &nearby_cellgrid);
            /*for (a1, a2) in nearby_naive.iter().zip(nearby_cellgrid.iter()) {
                println!("{:?} {:?}", a1, a2);
            }*/
        }
    }

    #[test]
    fn test_real_geometry_cylinder() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let radii = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 2.0];
        let heights = [1.0, 0.75, 2.0, 0.5, 3.0, 0.5, 16.0];
        let dimensions = [
            Dimension::Z,
            Dimension::X,
            Dimension::X,
            Dimension::Y,
            Dimension::Z,
            Dimension::Y,
        ];
        let cell_sizes = [0.5, 0.75, 1.0, 0.5, 1.0, 0.5, 2.0];

        for (radius, height, dimension, cell_size) in izip!(radii, heights, dimensions, cell_sizes)
        {
            let selected_atoms = system
                .get_atoms()
                .choose_multiple(&mut rng, 1000)
                .map(|atom| atom.get_index())
                .collect::<Vec<usize>>();

            // ignore the error
            let _ = system.group_create_from_indices("Group", selected_atoms);
            assert_eq!(system.group_get_n_atoms("Group").unwrap(), 1000);

            let nearby_naive = collect_nearby_cylinder_naive(&system, radius, height, dimension);
            let nearby_cellgrid =
                collect_nearby_cylinder_cellgrid(&system, radius, height, dimension, cell_size);

            compare_nearby_atoms(&nearby_naive, &nearby_cellgrid);
        }
    }

    #[test]
    fn test_real_geometry_rectangular() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let xs = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 15.0, 0.5, 1.0, 6.0];
        let ys = [1.0, 0.75, 0.5, 1.5, 1.0, 0.5, 1.0, 18.0, 2.0, 8.0];
        let zs = [2.0, 0.75, 1.0, 1.5, 3.0, 2.0, 3.0, 2.0, 13.0, 4.0];
        let cell_sizes = [0.5, 0.75, 0.5, 1.5, 1.0, 0.5, 1.0, 0.5, 1.0, 2.0];

        for (x, y, z, cell_size) in izip!(xs, ys, zs, cell_sizes) {
            let selected_atoms = system
                .get_atoms()
                .choose_multiple(&mut rng, 1000)
                .map(|atom| atom.get_index())
                .collect::<Vec<usize>>();

            // ignore the error
            let _ = system.group_create_from_indices("Group", selected_atoms);
            assert_eq!(system.group_get_n_atoms("Group").unwrap(), 1000);

            let nearby_naive = collect_nearby_rectangular_naive(&system, x, y, z);
            let nearby_cellgrid = collect_nearby_rectangular_cellgrid(&system, x, y, z, cell_size);

            compare_nearby_atoms(&nearby_naive, &nearby_cellgrid);
            /*for (a1, a2) in nearby_naive.iter().zip(nearby_cellgrid.iter()) {
                println!("{:?} {:?}", a1, a2);
            }*/
        }
    }
}
