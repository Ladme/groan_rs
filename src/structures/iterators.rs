// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of iterators over atoms and filter functions.

use crate::{
    errors::AtomError,
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
pub struct AtomIterator<'a> {
    atoms: &'a [Atom],
    container_iterator: AtomContainerIterator<'a>,
    simbox: Option<&'a SimBox>,
}

/// Iteration over atoms.
impl<'a> Iterator for AtomIterator<'a> {
    type Item = &'a Atom;

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
    /// - `simbox`: current dimensions of the simulation box (required only for geometric filtering)
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

impl<'a> MasterAtomIterator<'a> for AtomIterator<'a> {
    fn get_simbox(&self) -> Option<&SimBox> {
        self.simbox
    }
}

/// Immutable iterator over atoms. Same as `AtomIterator` but can be constructed for the system
/// without having to create a group. Constructed using `System::selection_iter()`.
/// Note that the `OwnedAtomIterator` still does NOT own the atoms of the system.
pub struct OwnedAtomIterator<'a> {
    atoms: &'a [Atom],
    container_iterator: OwnedAtomContainerIterator,
    simbox: Option<&'a SimBox>,
}

/// Iteration over atoms.
impl<'a> Iterator for OwnedAtomIterator<'a> {
    type Item = &'a Atom;

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
    /// - `simbox`: current dimensions of the simulation box (required only for geometric filtering)
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

impl<'a> MasterAtomIterator<'a> for OwnedAtomIterator<'a> {
    fn get_simbox(&self) -> Option<&SimBox> {
        self.simbox
    }
}

/// Immutable iterator over atoms with applied geometry filter.
/// Constructed by calling `filter_geometry` method on `AtomIterator` or `FilterAtomIterator`.
///
/// ## Notes
/// - Atoms with no positions set are never inside any geometric shape.
pub struct FilterAtomIterator<'a, I, S>
where
    I: Iterator<Item = &'a Atom>,
    S: Shape,
{
    iterator: I,
    geometry: S,
    simbox: SimBox,
}

impl<'a, I, S> MasterAtomIterator<'a> for FilterAtomIterator<'a, I, S>
where
    I: Iterator<Item = &'a Atom>,
    S: Shape,
{
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

    fn next(&mut self) -> Option<Self::Item> {
        self.iterator
            .find(|atom| atom.has_position() &&
                self.geometry.inside(atom
                    .get_position()
                    .expect("FATAL GROAN ERROR | MutFilterAtomIterator::next | Atom should have position."), 
                &self.simbox))
    }
}

/// Immutable iterator over atoms of a molecule.
/// Iterates over atoms respecting the topology of the molecule.
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
    /// - `simbox`: current dimensions of the simulation box (required only for geometric filtering)
    pub fn new(atoms: &'a [Atom], container: Vec<usize>, simbox: Option<&'a SimBox>) -> Self {
        MoleculeIterator {
            atoms,
            container,
            current_index: 0,
            simbox,
        }
    }
}

impl<'a> MasterAtomIterator<'a> for MoleculeIterator<'a> {
    fn get_simbox(&self) -> Option<&SimBox> {
        self.simbox
    }
}

/// Iteration over atoms of the molecule.
impl<'a> Iterator for MoleculeIterator<'a> {
    type Item = &'a Atom;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(index) = self.container.get(self.current_index) {
            self.current_index += 1;
            unsafe { Some(self.atoms.get_unchecked(*index)) }
        } else {
            None
        }
    }
}

/**************************/
/*    MUTABLE ITERATORS   */
/**************************/

/// Mutable iterator over atoms. Constructed using `System::atoms_iter_mut()` or `System::group_iter_mut()`.
/// Is guaranteed to iterate over atoms in the same order in which they are defined in the `System` structure.
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
    /// - `simbox`: current dimensions of the simulation box (required only for geometric filtering)
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

impl<'a> MasterMutAtomIterator<'a> for MutAtomIterator<'a> {
    fn get_simbox(&self) -> Option<&SimBox> {
        self.simbox
    }
}

impl<'a> Iterator for MutAtomIterator<'a> {
    type Item = &'a mut Atom;

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
    /// - `simbox`: current dimensions of the simulation box (required only for geometric filtering)
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

impl<'a> MasterMutAtomIterator<'a> for OwnedMutAtomIterator<'a> {
    fn get_simbox(&self) -> Option<&SimBox> {
        self.simbox
    }
}

impl<'a> Iterator for OwnedMutAtomIterator<'a> {
    type Item = &'a mut Atom;

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
pub struct MutFilterAtomIterator<'a, I, S>
where
    I: Iterator<Item = &'a mut Atom>,
    S: Shape,
{
    iterator: I,
    geometry: S,
    simbox: SimBox,
}

impl<'a, I, S> MasterMutAtomIterator<'a> for MutFilterAtomIterator<'a, I, S>
where
    I: Iterator<Item = &'a mut Atom>,
    S: Shape,
{
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

impl<'a> MasterMutAtomIterator<'a> for MutMoleculeIterator<'a> {
    fn get_simbox(&self) -> Option<&SimBox> {
        self.simbox
    }
}

/// Iteration over atoms of the molecule.
impl<'a> Iterator for MutMoleculeIterator<'a> {
    type Item = &'a mut Atom;

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
/* MASTER ITERATOR TRAITS */
/**************************/

/// Trait implemented by all immutable AtomIterators.
pub trait MasterAtomIterator<'a>: Iterator<Item = &'a Atom> + Sized {
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
}

pub trait MasterMutAtomIterator<'a>: Iterator<Item = &'a mut Atom> + Sized {
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
    /// or the simulation box is not orthogonal.
    /// - `AtomError::InvalidPosition` if any of the atoms of the iterator
    /// has an undefined position.
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
    /// or the simulation box is not orthogonal.
    /// - `AtomError::InvalidPosition` if any of the atoms of the iterator
    /// has an undefined position.
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
