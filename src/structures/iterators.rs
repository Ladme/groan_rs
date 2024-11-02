// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of iterators over atoms and filter functions.

use std::marker::PhantomData;

use crate::{
    auxiliary::PI_X2,
    errors::{AtomError, MassError, PositionError},
    structures::{
        atom::Atom,
        container::{AtomContainer, AtomContainerIterator, OwnedAtomContainerIterator},
        shape::Shape,
        simbox::SimBox,
    },
};

use super::{simbox::simbox_check, vector3d::Vector3D};

/**************************/
/*  IMMUTABLE ITERATORS   */
/**************************/

/// Immutable iterator over atoms. Constructed using `System::atoms_iter()` or `System::group_iter()`.
/// Is guaranteed to iterate over atoms in the same order in which they are defined in the `System` structure.
#[derive(Debug, Clone)]
pub struct AtomIterator<'a> {
    atoms: &'a [Atom],
    container_iterator: AtomContainerIterator<'a>,
    simbox: Option<&'a SimBox>,
}

/// Iteration over atoms.
impl<'a> Iterator for AtomIterator<'a> {
    type Item = &'a Atom;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(index) = self.container_iterator.next() {
            unsafe { Some(self.atoms.get_unchecked(index)) }
        } else {
            None
        }
    }
}

impl<'a> AtomIterator<'a> {
    /// Create a new `AtomIterator`.
    ///
    /// ## Parameters
    /// - `atoms`: reference to the atoms of the `System`
    /// - `atom_container`: `AtomContainer` specifying the atoms to iterate through
    /// - `simbox`: current dimensions of the simulation box
    pub fn new(
        atoms: &'a [Atom],
        atom_container: &'a AtomContainer,
        simbox: Option<&'a SimBox>,
    ) -> Self {
        AtomIterator {
            atoms,
            container_iterator: atom_container.iter(),
            simbox,
        }
    }
}

impl<'a> AtomIteratorWithBox<'a> for AtomIterator<'a> {
    #[inline(always)]
    fn get_simbox(&self) -> Option<&SimBox> {
        self.simbox
    }
}

impl<'a> OrderedAtomIterator<'a> for AtomIterator<'a> {}

/// Immutable iterator over atoms. Same as `AtomIterator` but can be constructed for the system
/// without having to create a group. Constructed using `System::selection_iter()`.
/// Note that the `OwnedAtomIterator` still does NOT own the atoms of the system.
#[derive(Debug, Clone)]
pub struct OwnedAtomIterator<'a> {
    atoms: &'a [Atom],
    container_iterator: OwnedAtomContainerIterator,
    simbox: Option<&'a SimBox>,
}

/// Iteration over atoms.
impl<'a> Iterator for OwnedAtomIterator<'a> {
    type Item = &'a Atom;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(index) = self.container_iterator.next() {
            unsafe { Some(self.atoms.get_unchecked(index)) }
        } else {
            None
        }
    }
}

impl<'a> OwnedAtomIterator<'a> {
    /// Create a new `OwnedAtomIterator`.
    ///
    /// ## Parameters
    /// - `atoms`: reference to the atoms of the `System`
    /// - `atom_container`: `AtomContainer` specifying the atoms to iterate through
    /// - `simbox`: current dimensions of the simulation box
    pub fn new(
        atoms: &'a [Atom],
        atom_container: AtomContainer,
        simbox: Option<&'a SimBox>,
    ) -> Self {
        OwnedAtomIterator {
            atoms,
            container_iterator: atom_container.into_iter(),
            simbox,
        }
    }
}

impl<'a> AtomIteratorWithBox<'a> for OwnedAtomIterator<'a> {
    #[inline(always)]
    fn get_simbox(&self) -> Option<&SimBox> {
        self.simbox
    }
}

impl<'a> OrderedAtomIterator<'a> for OwnedAtomIterator<'a> {}

/// Immutable iterator over atoms with applied geometry filter.
/// Constructed by calling `filter_geometry` method on `AtomIterator` or `FilterAtomIterator`.
///
/// ## Notes
/// - Atoms with no positions set are never inside any geometric shape.
#[derive(Debug, Clone)]
pub struct FilterAtomIterator<'a, I, S>
where
    I: Iterator<Item = &'a Atom>,
    S: Shape,
{
    iterator: I,
    geometry: S,
    simbox: SimBox,
}

impl<'a, I, S> AtomIteratorWithBox<'a> for FilterAtomIterator<'a, I, S>
where
    I: Iterator<Item = &'a Atom>,
    S: Shape,
{
    #[inline(always)]
    fn get_simbox(&self) -> Option<&SimBox> {
        Some(&self.simbox)
    }
}

impl<'a, I, S> Iterator for FilterAtomIterator<'a, I, S>
where
    I: Iterator<Item = &'a Atom>,
    S: Shape,
{
    type Item = &'a Atom;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.iterator.find(|atom| {
            atom.has_position()
                && self.geometry.inside(
                    atom.get_position().expect(
                        "FATAL GROAN ERROR | FilterAtomIterator::next | Atom should have position.",
                    ),
                    &self.simbox,
                )
        })
    }
}

impl<'a, I, S> OrderedAtomIterator<'a> for FilterAtomIterator<'a, I, S>
where
    I: Iterator<Item = &'a Atom>,
    S: Shape,
{
}

/// Immutable iterator over atoms of a molecule.
/// Iterates over atoms respecting the topology of the molecule.
#[derive(Debug, Clone)]
pub struct MoleculeIterator<'a> {
    atoms: &'a [Atom],
    container: Vec<usize>,
    current_index: usize,
    simbox: Option<&'a SimBox>,
}

impl<'a> MoleculeIterator<'a> {
    /// Create a new `MoleculeIterator`
    ///
    /// ## Parameters
    /// - `atoms`: reference to the atoms of the `System`
    /// - `container`: vector specifying indices of the atoms to iterate through
    /// - `simbox`: current dimensions of the simulation box
    pub fn new(atoms: &'a [Atom], container: Vec<usize>, simbox: Option<&'a SimBox>) -> Self {
        MoleculeIterator {
            atoms,
            container,
            current_index: 0,
            simbox,
        }
    }
}

impl<'a> AtomIteratorWithBox<'a> for MoleculeIterator<'a> {
    #[inline(always)]
    fn get_simbox(&self) -> Option<&SimBox> {
        self.simbox
    }
}

/// Iteration over atoms of the molecule.
impl<'a> Iterator for MoleculeIterator<'a> {
    type Item = &'a Atom;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(index) = self.container.get(self.current_index) {
            self.current_index += 1;
            unsafe { Some(self.atoms.get_unchecked(*index)) }
        } else {
            None
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum IteratorOrigin {
    First,
    Second,
}

/// Structure storing an atom and what iterator it was yielded by.
#[derive(Debug, Clone)]
struct AtomOrigin<'a> {
    atom: &'a Atom,
    origin: IteratorOrigin,
}

impl<'a> AtomOrigin<'a> {
    fn try_new(atom: Option<&Atom>, origin: IteratorOrigin) -> Option<AtomOrigin> {
        match atom {
            Some(x) => Some(AtomOrigin {
                atom: x,
                origin: origin,
            }),
            None => None,
        }
    }
}

/// Sort two atoms based on their index.
fn sort_atoms<'a>(a: AtomOrigin<'a>, b: AtomOrigin<'a>) -> (AtomOrigin<'a>, AtomOrigin<'a>, bool) {
    if a.atom.get_index() < b.atom.get_index() {
        (a, b, false)
    } else if a.atom.get_index() > b.atom.get_index() {
        (b, a, false)
    } else {
        (a, b, true)
    }
}

/// Immutable iterator over atoms of two iterators.
/// The atoms are guaranteed to be yielded in the order in which they are defined in the System.
/// Every atom will be yielded at most once.
#[derive(Debug)]
pub struct UnionAtomIterator<'a, I1, I2>
where
    I1: OrderedAtomIterator<'a> + AtomIteratorWithBox<'a>,
    I2: OrderedAtomIterator<'a> + AtomIteratorWithBox<'a>,
{
    iterator_1: I1,
    iterator_2: I2,
    buffer: Option<AtomOrigin<'a>>,
}

impl<'a, I1, I2> Iterator for UnionAtomIterator<'a, I1, I2>
where
    I1: OrderedAtomIterator<'a> + AtomIteratorWithBox<'a>,
    I2: OrderedAtomIterator<'a> + AtomIteratorWithBox<'a>,
{
    type Item = &'a Atom;

    fn next(&mut self) -> Option<Self::Item> {
        match self.buffer.clone() {
            // if we have any value stored in the buffer, we propagate the other iterator and get a second value
            Some(a) => {
                let b = match a.origin {
                    IteratorOrigin::First => {
                        AtomOrigin::try_new(self.iterator_2.next(), IteratorOrigin::Second)
                    }
                    IteratorOrigin::Second => {
                        AtomOrigin::try_new(self.iterator_1.next(), IteratorOrigin::First)
                    }
                };

                match b {
                    // compare the new yielded value with the buffer, yield the atom with lower index and store the other in the buffer
                    Some(b) => {
                        let (min, max, same) = sort_atoms(a, b);

                        // store the other atom into the buffer
                        if !same {
                            self.buffer = Some(max);
                        // clear the buffer if both atoms have the same index (can't be yielded twice)
                        } else {
                            self.buffer = None;
                        }

                        Some(min.atom)
                    }
                    // we return the value from the buffer and clear the buffer the buffer
                    None => {
                        self.buffer = None;
                        Some(a.atom)
                    }
                }
            }

            // if there is no value in the buffer, propagate both iterators and select the lower atom from them
            None => {
                let first = AtomOrigin::try_new(self.iterator_1.next(), IteratorOrigin::First);
                let second = AtomOrigin::try_new(self.iterator_2.next(), IteratorOrigin::Second);

                match (first, second) {
                    (None, None) => None,
                    (Some(a), None) => Some(a.atom), // buffer is already empty
                    (None, Some(b)) => Some(b.atom),
                    (Some(a), Some(b)) => {
                        let (min, max, same) = sort_atoms(a, b);

                        // store the other atom into the buffer
                        if !same {
                            self.buffer = Some(max);
                        }

                        Some(min.atom)
                    }
                }
            }
        }
    }
}

impl<'a, I1, I2> AtomIteratorWithBox<'a> for UnionAtomIterator<'a, I1, I2>
where
    I1: OrderedAtomIterator<'a> + AtomIteratorWithBox<'a>,
    I2: OrderedAtomIterator<'a> + AtomIteratorWithBox<'a>,
{
    #[inline(always)]
    fn get_simbox(&self) -> Option<&SimBox> {
        self.iterator_1.get_simbox()
    }
}

impl<'a, I1, I2> OrderedAtomIterator<'a> for UnionAtomIterator<'a, I1, I2>
where
    I1: OrderedAtomIterator<'a> + AtomIteratorWithBox<'a>,
    I2: OrderedAtomIterator<'a> + AtomIteratorWithBox<'a>,
{
}

/// Immutable iterator over atoms shared by two iterators.
/// The atoms are guaranteed to be yielded in the order in which they are defined in the System.
/// Every atom will be yielded at most once.
pub struct IntersectionAtomIterator {}

/**************************/
/*    MUTABLE ITERATORS   */
/**************************/

/// Mutable iterator over atoms. Constructed using `System::atoms_iter_mut()` or `System::group_iter_mut()`.
/// Is guaranteed to iterate over atoms in the same order in which they are defined in the `System` structure.
#[derive(Debug, Clone)]
pub struct MutAtomIterator<'a> {
    atoms: *mut [Atom],
    container_iterator: AtomContainerIterator<'a>,
    simbox: Option<&'a SimBox>,
}

impl<'a> MutAtomIterator<'a> {
    /// Create a new `MutAtomIterator`.
    ///
    /// ## Parameters
    /// - `atoms`: mutable reference to the atoms of the `System`
    /// - `atom_container`: `AtomContainer` specifying the atoms to iterate through
    /// - `simbox`: current dimensions of the simulation box
    pub fn new(
        atoms: &'a mut [Atom],
        atom_container: &'a AtomContainer,
        simbox: Option<&'a SimBox>,
    ) -> Self {
        MutAtomIterator {
            atoms: atoms as *mut [Atom],
            container_iterator: atom_container.iter(),
            simbox,
        }
    }
}

impl<'a> MutAtomIteratorWithBox<'a> for MutAtomIterator<'a> {
    #[inline(always)]
    fn get_simbox(&self) -> Option<&SimBox> {
        self.simbox
    }
}

impl<'a> Iterator for MutAtomIterator<'a> {
    type Item = &'a mut Atom;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(index) = self.container_iterator.next() {
            unsafe { Some((*self.atoms).get_unchecked_mut(index)) }
        } else {
            None
        }
    }
}

/// Mutable iterator over atoms. Same as `MutAtomIterator` but can be constructed for the system
/// without having to create a group. Constructed using `System::selection_iter_mut()`.
/// Note that the `OwnedMutAtomIterator` still does NOT own the atoms of the system.
#[derive(Debug, Clone)]
pub struct OwnedMutAtomIterator<'a> {
    atoms: *mut [Atom],
    container_iterator: OwnedAtomContainerIterator,
    simbox: Option<&'a SimBox>,
}

impl<'a> OwnedMutAtomIterator<'a> {
    /// Create a new `MutOwnedAtomIterator`.
    ///
    /// ## Parameters
    /// - `atoms`: mutable reference to the atoms of the `System`
    /// - `atom_container`: `AtomContainer` specifying the atoms to iterate through
    /// - `simbox`: current dimensions of the simulation box
    pub fn new(
        atoms: &'a mut [Atom],
        atom_container: AtomContainer,
        simbox: Option<&'a SimBox>,
    ) -> Self {
        OwnedMutAtomIterator {
            atoms: atoms as *mut [Atom],
            container_iterator: atom_container.into_iter(),
            simbox,
        }
    }
}

impl<'a> MutAtomIteratorWithBox<'a> for OwnedMutAtomIterator<'a> {
    #[inline(always)]
    fn get_simbox(&self) -> Option<&SimBox> {
        self.simbox
    }
}

impl<'a> Iterator for OwnedMutAtomIterator<'a> {
    type Item = &'a mut Atom;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(index) = self.container_iterator.next() {
            unsafe { Some((*self.atoms).get_unchecked_mut(index)) }
        } else {
            None
        }
    }
}

/// Mutable iterator over atoms with applied geometry filter.
/// Constructed by calling `filter_geometry` method on `MutAtomIterator` or `MutFilterAtomIterator`.
///
/// ## Notes
/// - Atoms with no positions set are never inside any geometric shape.
#[derive(Debug, Clone)]
pub struct MutFilterAtomIterator<'a, I, S>
where
    I: Iterator<Item = &'a mut Atom>,
    S: Shape,
{
    iterator: I,
    geometry: S,
    simbox: SimBox,
}

impl<'a, I, S> MutAtomIteratorWithBox<'a> for MutFilterAtomIterator<'a, I, S>
where
    I: Iterator<Item = &'a mut Atom>,
    S: Shape,
{
    #[inline(always)]
    fn get_simbox(&self) -> Option<&SimBox> {
        Some(&self.simbox)
    }
}

impl<'a, I, S> Iterator for MutFilterAtomIterator<'a, I, S>
where
    I: Iterator<Item = &'a mut Atom>,
    S: Shape,
{
    type Item = &'a mut Atom;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.iterator
            .find(|atom| atom.has_position() &&
                self.geometry.inside(atom
                    .get_position()
                    .expect("FATAL GROAN ERROR | MutFilterAtomIterator::next | Atom should have position."), 
                &self.simbox))
    }
}

/// Mutable iterator over atoms of a molecule.
/// Iterates over atoms respecting the topology of the molecule.
#[derive(Debug, Clone)]
pub struct MutMoleculeIterator<'a> {
    atoms: *mut [Atom],
    container: Vec<usize>,
    current_index: usize,
    simbox: Option<&'a SimBox>,
}

impl<'a> MutMoleculeIterator<'a> {
    /// Create a new `MoleculeIterator`
    ///
    /// ## Parameters
    /// - `atoms`: mutable reference to the atoms of the `System`
    /// - `container`: vector specifying indices of the atoms to iterate through
    /// - `simbox`: current dimensions of the simulation box
    pub fn new(atoms: &'a mut [Atom], container: Vec<usize>, simbox: Option<&'a SimBox>) -> Self {
        MutMoleculeIterator {
            atoms,
            container,
            current_index: 0,
            simbox,
        }
    }
}

impl<'a> MutAtomIteratorWithBox<'a> for MutMoleculeIterator<'a> {
    #[inline(always)]
    fn get_simbox(&self) -> Option<&SimBox> {
        self.simbox
    }
}

/// Iteration over atoms of the molecule.
impl<'a> Iterator for MutMoleculeIterator<'a> {
    type Item = &'a mut Atom;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(index) = self.container.get(self.current_index) {
            self.current_index += 1;
            unsafe { Some((*self.atoms).get_unchecked_mut(*index)) }
        } else {
            None
        }
    }
}

/**************************/
/*     ITERATOR TRAITS    */
/**************************/

/// Trait implemented by all immutable iterators over atoms that contain information about the simulation box.
pub trait AtomIteratorWithBox<'a>: Iterator<Item = &'a Atom> + Sized {
    /// Get reference to the simulation box inside the iterator.
    /// Returns an option.
    fn get_simbox(&self) -> Option<&SimBox>;

    /// Get reference to the simulation box inside the iterator.
    ///
    /// ## Panics
    /// Panics if the simulation box is undefined.
    #[inline(always)]
    fn get_simbox_unwrap(&self) -> &SimBox {
        match self.get_simbox() {
            Some(x) => x,
            None => panic!("FATAL GROAN ERROR | MasterAtomIterator::get_simbox_unwrap | No simulation box associated with the atom iterator.")
        }
    }

    /// Filter atoms located inside the specified geometric shape.
    ///
    /// ## Panics
    /// Panics if the iterator has no associated simulation box.
    ///
    /// ## Example
    /// Iterating over all atoms of the system
    /// that are located in a sphere around a specific point.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let system = System::from_file("system.gro").unwrap();
    ///
    /// // construct a sphere located at x = 1, y = 2, z = 3 with a radius of 2.5 nm
    /// let sphere = Sphere::new([1.0, 2.0, 3.0].into(), 2.5);
    ///
    /// for atom in system.atoms_iter().filter_geometry(sphere) {
    ///     println!("{:?}", atom);
    /// }
    /// ```
    fn filter_geometry(self, geometry: impl Shape) -> FilterAtomIterator<'a, Self, impl Shape> {
        let simbox = self.get_simbox_unwrap().clone();

        FilterAtomIterator {
            iterator: self,
            geometry,
            simbox,
        }
    }

    /// Calculate center of geometry of a group of atoms selected by an iterator.
    /// Takes periodic boundary conditions into consideration.
    /// Useful for the calculation of local center of geometry.
    ///
    /// ## Returns
    /// - `Vector3D` corresponding to the geometric center of the selected atoms.
    /// - `AtomError::InvalidSimBox` if the iterator has no simulation box
    ///   or the simulation box is not orthogonal.
    /// - `AtomError::InvalidPosition` if any of the atoms of the iterator has no position.
    ///
    /// ## Notes
    /// - This calculation approach is adapted from Linge Bai & David Breen (2008).
    /// - It is able to calculate correct center of geometry for any distribution of atoms
    ///   that is not completely homogeneous.
    /// - In case the iterator is empty, the center of geometry is NaN.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// let sphere = Sphere::new(Vector3D::new(1.0, 2.0, 3.0), 2.5);
    ///
    /// // select atoms of the group "Group" located inside a sphere with
    /// // a radius of 2.5 nm centered at [1.0, 2.0, 3.0]
    /// let atoms = system
    ///     .group_iter("Group")
    ///     .unwrap()
    ///     .filter_geometry(sphere);
    ///
    /// // calculate center of geometry of "atoms"
    /// let center = match atoms.get_center() {
    ///     Ok(x) => x,
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;    
    ///     }
    /// };
    /// ```
    fn get_center(self) -> Result<Vector3D, AtomError> {
        let simbox = simbox_check(self.get_simbox()).map_err(AtomError::InvalidSimBox)?;
        let scaling = Vector3D::new(PI_X2 / simbox.x, PI_X2 / simbox.y, PI_X2 / simbox.z);
        let simbox = simbox as *const SimBox;

        let mut sum_xi = Vector3D::default();
        let mut sum_zeta = Vector3D::default();

        let mut empty = true;
        for atom in self {
            match atom.get_position() {
                Some(x) => unsafe {
                    crate::auxiliary::center_atom_contribution(
                        x.clone(),
                        &scaling,
                        &*simbox,
                        1.0, // use 1 as mass since we are calculating center of geometry
                        &mut sum_xi,
                        &mut sum_zeta,
                    )
                },
                None => {
                    return Err(AtomError::InvalidPosition(PositionError::NoPosition(
                        atom.get_index(),
                    )));
                }
            }

            empty = false;
        }

        if empty {
            return Ok(Vector3D::new(f32::NAN, f32::NAN, f32::NAN));
        }

        // convert to real coordinates
        Ok(crate::auxiliary::from_circle_to_line(
            sum_zeta, sum_xi, &scaling,
        ))
    }

    /// Calculate center of mass of a group of atoms selected by an iterator.
    /// Takes periodic boundary conditions into consideration.
    /// Useful for the calculation of local center of geometry.
    ///
    /// ## Returns
    /// - `Vector3D` corresponding to the center of mass of the selected atoms.
    /// - `AtomError::InvalidSimBox` if the iterator has no simulation box
    ///   or the simulation box is not orthogonal.
    /// - `AtomError::InvalidPosition` if any of the atoms of the iterator has no position.
    /// - `AtomError::InvalidMass` if any of the atoms of the iterator has no mass.
    ///
    /// ## Notes
    /// - This calculation approach is adapted from Linge Bai & David Breen (2008).
    /// - It is able to calculate correct center of mass for any distribution of atoms
    ///   that is not completely homogeneous.
    /// - In case the iterator is empty, the center of mass is NaN.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.read_ndx("index.ndx").unwrap();
    ///
    /// let sphere = Sphere::new(Vector3D::new(1.0, 2.0, 3.0), 2.5);
    ///
    /// // select atoms of the group "Group" located inside a sphere with
    /// // a radius of 2.5 nm centered at [1.0, 2.0, 3.0]
    /// let atoms = system
    ///     .group_iter("Group")
    ///     .unwrap()
    ///     .filter_geometry(sphere);
    ///
    /// // calculate center of mass of "atoms"
    /// let center = match atoms.get_com() {
    ///     Ok(x) => x,
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;    
    ///     }
    /// };
    /// ```
    fn get_com(self) -> Result<Vector3D, AtomError> {
        let simbox = simbox_check(self.get_simbox()).map_err(AtomError::InvalidSimBox)?;
        let scaling = Vector3D::new(PI_X2 / simbox.x, PI_X2 / simbox.y, PI_X2 / simbox.z);
        let simbox = simbox as *const SimBox;

        let mut sum_xi = Vector3D::default();
        let mut sum_zeta = Vector3D::default();

        let mut empty = true;
        for atom in self {
            let mass = atom
                .get_mass()
                .ok_or(AtomError::InvalidMass(MassError::NoMass(atom.get_index())))?;

            match atom.get_position() {
                Some(x) => unsafe {
                    crate::auxiliary::center_atom_contribution(
                        x.clone(),
                        &scaling,
                        &*simbox,
                        mass,
                        &mut sum_xi,
                        &mut sum_zeta,
                    )
                },
                None => {
                    return Err(AtomError::InvalidPosition(PositionError::NoPosition(
                        atom.get_index(),
                    )));
                }
            }

            empty = false;
        }

        if empty {
            return Ok(Vector3D::new(f32::NAN, f32::NAN, f32::NAN));
        }

        // convert to real coordinates
        Ok(crate::auxiliary::from_circle_to_line(
            sum_zeta, sum_xi, &scaling,
        ))
    }
}

/// Trait implemented by all mutable iterators over atoms that contain information about the simulation box.
pub trait MutAtomIteratorWithBox<'a>: Iterator<Item = &'a mut Atom> + Sized {
    /// Get reference to the simulation box inside the iterator.
    /// Returns an option.
    fn get_simbox(&self) -> Option<&SimBox>;

    /// Get reference to the simulation box inside the iterator.
    ///
    /// ## Panics
    /// Panics if the simulation box is undefined.
    fn get_simbox_unwrap(&self) -> &SimBox {
        match self.get_simbox() {
            Some(x) => x,
            None => panic!("FATAL GROAN ERROR | MasterMutAtomIterator::get_simbox_unwrap | No simulation box associated with the atom iterator.")
        }
    }

    /// Filter atoms located inside the specified geometric shape.
    ///
    /// ## Panics
    /// Panics if the iterator has no associated simulation box.
    ///
    /// ## Example
    /// Iterating over all atoms of the system
    /// that are located in a sphere around a specific point
    /// and change their names.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// // construct a sphere located at x = 1, y = 2, z = 3 with a radius of 2.5 nm
    /// let sphere = Sphere::new([1.0, 2.0, 3.0].into(), 2.5);
    ///
    /// for atom in system.atoms_iter_mut().filter_geometry(sphere) {
    ///     atom.set_atom_name("XYZ");
    /// }
    /// ```
    fn filter_geometry(self, geometry: impl Shape) -> MutFilterAtomIterator<'a, Self, impl Shape> {
        let simbox = self.get_simbox_unwrap().clone();

        MutFilterAtomIterator {
            iterator: self,
            geometry,
            simbox,
        }
    }

    /// Translate atoms of the iterator by target vector.
    ///
    /// ## Returns
    /// - `Ok` if everything was successful.
    /// - `AtomError::InvalidSimBox` if the iterator has no simulation box
    ///   or the simulation box is not orthogonal.
    /// - `AtomError::InvalidPosition` if any of the atoms of the iterator
    ///   has an undefined position.
    ///
    /// ## Example
    /// Translate all atoms with the name CA by 1.3 nm in the x-dimension.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system
    ///     .selection_iter_mut("name CA")
    ///     .unwrap()
    ///     .translate(&Vector3D::new(1.3, 0.0, 0.0))
    ///     .unwrap();
    /// ```
    fn translate(mut self, vector: &Vector3D) -> Result<(), AtomError> {
        let simbox =
            simbox_check(self.get_simbox()).map_err(AtomError::InvalidSimBox)? as *const SimBox;

        unsafe { self.try_for_each(|atom| atom.translate(vector, &*simbox)) }
    }

    /// Wrap atoms of the iterator into the simulation box.
    ///
    /// ## Returns
    /// `Ok` if everything was successful.
    /// - `AtomError::InvalidSimBox` if the iterator has no simulation box
    ///   or the simulation box is not orthogonal.
    /// - `AtomError::InvalidPosition` if any of the atoms of the iterator
    ///   has an undefined position.
    ///
    /// ## Example
    /// Wrap all atoms with the name CA into the simulation box.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system
    ///     .selection_iter_mut("name CA")
    ///     .unwrap()
    ///     .wrap()
    ///     .unwrap();
    /// ```
    fn wrap(mut self) -> Result<(), AtomError> {
        let simbox =
            simbox_check(self.get_simbox()).map_err(AtomError::InvalidSimBox)? as *const SimBox;

        unsafe { self.try_for_each(|atom| atom.wrap(&*simbox)) }
    }
}

/// Trait implemented by immutable iterators which yield atoms in the order
/// in which they are defined in the parent System.
pub trait OrderedAtomIterator<'a>: Iterator<Item = &'a Atom> + AtomIteratorWithBox<'a> {
    fn union<T>(self, other: T) -> UnionAtomIterator<'a, Self, T>
    where
        Self: Sized,
        T: OrderedAtomIterator<'a> + Sized,
    {
        UnionAtomIterator {
            iterator_1: self,
            iterator_2: other,
            buffer: None,
        }
    }

    /*fn intersection(self, other: impl OrderedAtomIterator<'a>) -> IntersectionAtomIterator<'a>
    where
        Self: Sized,
    {
        IntersectionAtomIterator {
            iterator_1: Box::new(self),
            iterator_2: Box::new(other),
        }
    }*/
}

/**************************/
/*     PAIR ITERATORS     */
/**************************/

/// Iterator over pairs of immutable atoms in specified order.
#[derive(Debug, Clone)]
pub struct AtomPairIterator<'a> {
    atoms: &'a [Atom],
    container: Vec<(usize, usize)>,
    current_index: usize,
}

impl<'a> AtomPairIterator<'a> {
    pub fn new(atoms: &'a [Atom], container: Vec<(usize, usize)>) -> Self {
        AtomPairIterator {
            atoms,
            container,
            current_index: 0,
        }
    }
}

impl<'a> Iterator for AtomPairIterator<'a> {
    type Item = (&'a Atom, &'a Atom);

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some((i, j)) = self.container.get(self.current_index) {
            self.current_index += 1;
            unsafe { Some((self.atoms.get_unchecked(*i), self.atoms.get_unchecked(*j))) }
        } else {
            None
        }
    }
}

/// Iterator over pairs of mutable atoms in specified order.
#[derive(Debug, Clone)]
pub struct MutAtomPairIterator<'a> {
    atoms: *mut [Atom],
    container: Vec<(usize, usize)>,
    current_index: usize,
    phantom: PhantomData<&'a Atom>, // remove, if we ever use simbox here
}

impl<'a> MutAtomPairIterator<'a> {
    pub fn new(atoms: &mut [Atom], container: Vec<(usize, usize)>) -> Self {
        MutAtomPairIterator {
            atoms: atoms as *mut [Atom],
            container,
            current_index: 0,
            phantom: PhantomData,
        }
    }
}

impl<'a> Iterator for MutAtomPairIterator<'a> {
    type Item = (&'a mut Atom, &'a mut Atom);

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        if let Some((i, j)) = self.container.get(self.current_index) {
            self.current_index += 1;
            unsafe {
                Some((
                    (*self.atoms).get_unchecked_mut(*i),
                    (*self.atoms).get_unchecked_mut(*j),
                ))
            }
        } else {
            None
        }
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use float_cmp::assert_approx_eq;

    use crate::{
        prelude::{Cylinder, Dimension},
        structures::{element::Elements, shape::Sphere},
        system::System,
    };

    use super::*;

    #[test]
    fn iterator_get_center() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let sphere_pos = system.group_get_center("Protein").unwrap();
        let sphere = Sphere::new(sphere_pos.clone(), 2.0);

        let center = system
            .group_iter("Membrane")
            .unwrap()
            .filter_geometry(sphere)
            .get_center()
            .unwrap();

        assert_approx_eq!(f32, center.x, 9.8453);
        assert_approx_eq!(f32, center.y, 2.4803874);
        assert_approx_eq!(f32, center.z, 5.434977);
    }

    #[test]
    fn iterator_get_center_empty() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.group_create("EmptyGroup", "not all").unwrap();

        let center = system
            .group_iter("EmptyGroup")
            .unwrap()
            .get_center()
            .unwrap();

        assert!(center.x.is_nan());
        assert!(center.y.is_nan());
        assert!(center.z.is_nan());
    }

    #[test]
    fn iterator_get_com() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();

        system.group_create("Peptide", "@protein").unwrap();
        system.group_create("Membrane", "@membrane").unwrap();

        system.guess_elements(Elements::default()).unwrap();

        let sphere_pos = system.group_get_center("Peptide").unwrap();
        let sphere = Sphere::new(sphere_pos.clone(), 1.0);

        let com = system
            .group_iter("Membrane")
            .unwrap()
            .filter_geometry(sphere)
            .get_com()
            .unwrap();

        assert_approx_eq!(f32, com.x, 4.0072813);
        assert_approx_eq!(f32, com.y, 3.7480402);
        assert_approx_eq!(f32, com.z, 3.3228612);
    }

    #[test]
    fn iterator_get_com_empty() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.group_create("EmptyGroup", "not all").unwrap();

        system.guess_elements(Elements::default()).unwrap();

        let center = system.group_iter("EmptyGroup").unwrap().get_com().unwrap();

        assert!(center.x.is_nan());
        assert!(center.y.is_nan());
        assert!(center.z.is_nan());
    }

    #[test]
    fn iterator_translate() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system
            .selection_iter_mut("resname ALA")
            .unwrap()
            .translate(&Vector3D::new(3.5, -1.1, 5.4))
            .unwrap();

        let first = system.get_atom(31).unwrap().get_position().unwrap();
        let last = system.get_atom(52).unwrap().get_position().unwrap();

        assert_approx_eq!(f32, first.x, 0.23069);
        assert_approx_eq!(f32, first.y, 1.567);
        assert_approx_eq!(f32, first.z, 10.745);

        assert_approx_eq!(f32, last.x, 0.28168964);
        assert_approx_eq!(f32, last.y, 1.231);
        assert_approx_eq!(f32, last.z, 9.237);
    }

    #[test]
    fn iterator_wrap() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.atoms_iter_mut().for_each(|atom| {
            atom.translate_nopbc(&Vector3D::new(1000.0, 1000.0, 1000.0))
                .unwrap()
        });

        system.group_create("Alanines", "resname ALA").unwrap();
        system.group_iter_mut("Alanines").unwrap().wrap().unwrap();

        let simbox = system.get_box().unwrap();
        for (a, atom) in system.atoms_iter().enumerate() {
            let pos = atom.get_position().unwrap();

            if system.group_isin("Alanines", a).unwrap() {
                assert!(pos.x <= simbox.x);
                assert!(pos.y <= simbox.y);
                assert!(pos.z <= simbox.z);
            } else {
                assert!(pos.x >= 1000.0);
                assert!(pos.y >= 1000.0);
                assert!(pos.z >= 1000.0);
            }
        }
    }

    #[test]
    fn iterator_union() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.group_create("None", "not all").unwrap();
        let iterator1 = system.group_iter("None").unwrap();
        let iterator2 = system.group_iter("None").unwrap();
        assert!(iterator1.union(iterator2).next().is_none());

        let iterator1 = system.group_iter("None").unwrap();
        let iterator2 = system.selection_iter("serial 1 2 3 7 8 13").unwrap();
        let expected = [1, 2, 3, 7, 8, 13];
        for (atom, expected) in iterator1.union(iterator2).zip(expected.into_iter()) {
            assert_eq!(atom.get_atom_number(), expected);
        }

        system.group_create("Group", "serial 1 2 3 7 8 13").unwrap();
        let iterator1 = system.group_iter("Group").unwrap();
        let iterator2 = system.selection_iter("serial 1 2 3 7 8 13").unwrap();
        let expected = [1, 2, 3, 7, 8, 13];
        for (atom, expected) in iterator1.union(iterator2).zip(expected.into_iter()) {
            assert_eq!(atom.get_atom_number(), expected);
        }

        let iterator1 = system.selection_iter("serial 1 2 3 7 8 13").unwrap();
        let iterator2 = system
            .selection_iter("serial 10 11 12 13 14 5 6 7 8")
            .unwrap();
        let expected = [1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13, 14];
        for (atom, expected) in iterator1.union(iterator2).zip(expected.into_iter()) {
            assert_eq!(atom.get_atom_number(), expected);
        }
    }

    #[test]
    fn iterator_union_filter_geometry() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let iterator1 = system.selection_iter("@membrane").unwrap();
        let iterator2 = system.selection_iter("@water").unwrap();

        let cylinder = Cylinder::new(
            system.group_get_center("Protein").unwrap(),
            2.5,
            4.0,
            Dimension::Z,
        );

        for (a1, a2) in iterator1
            .union(iterator2)
            .filter_geometry(cylinder.clone())
            .zip(
                system
                    .selection_iter("@membrane or @water")
                    .unwrap()
                    .filter_geometry(cylinder),
            )
        {
            assert_eq!(a1.get_index(), a2.get_index());
        }
    }

    #[test]
    fn iterator_filter_geometry_union() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        let zcylinder = Cylinder::new(
            system.group_get_center("Protein").unwrap(),
            2.5,
            4.0,
            Dimension::Z,
        );

        let xcylinder = Cylinder::new(
            system.group_get_center("Protein").unwrap(),
            3.5,
            2.0,
            Dimension::X,
        );

        system
            .group_create_from_geometry("Zcylinder", "@membrane", zcylinder.clone())
            .unwrap();

        system
            .group_create_from_geometry("Xcylinder", "@membrane", xcylinder.clone())
            .unwrap();

        system
            .group_union("Xcylinder", "Zcylinder", "Geometry")
            .unwrap();

        let iterator1 = system
            .selection_iter("@membrane")
            .unwrap()
            .filter_geometry(zcylinder);
        let iterator2 = system
            .selection_iter("@membrane")
            .unwrap()
            .filter_geometry(xcylinder);

        for (a1, a2) in iterator1
            .union(iterator2)
            .zip(system.group_iter("Geometry").unwrap())
        {
            assert_eq!(a1.get_index(), a2.get_index());
        }
    }

    #[test]
    fn iterator_union_union() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system
            .group_create("ProtMemWat", "@protein or @membrane or @water")
            .unwrap();

        let iterator1 = system.selection_iter("@protein").unwrap();
        let iterator2 = system.selection_iter("@water").unwrap();
        let iterator3 = system.selection_iter("@membrane").unwrap();

        for (a1, a2) in iterator1
            .clone()
            .union(iterator2.clone())
            .union(iterator3.clone())
            .zip(system.group_iter("ProtMemWat").unwrap())
        {
            assert_eq!(a1.get_index(), a2.get_index());
        }

        for (a1, a2) in iterator3
            .clone()
            .union(iterator1.clone())
            .union(iterator2.clone())
            .zip(system.group_iter("ProtMemWat").unwrap())
        {
            assert_eq!(a1.get_index(), a2.get_index());
        }
    }
}
