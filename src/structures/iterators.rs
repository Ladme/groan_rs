// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of iterators over atoms and filter functions.

use crate::structures::{
    atom::Atom,
    container::{AtomContainer, AtomContainerIterator},
    shape::Shape,
    simbox::SimBox,
};

/**************************/
/*  IMMUTABLE ITERATORS   */
/**************************/

/// Immutable iterator over atoms. Constructed using `System::atoms_iter()` or `System::group_iter()`.
/// Is guaranteed to iterate over atoms in the same order in which they are defined in the `System` structure.
pub struct AtomIterator<'a> {
    atoms: &'a [Atom],
    container_iterator: AtomContainerIterator<'a>,
    simbox: &'a SimBox,
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
    /// - `simbox`: current dimensions of the simulation box
    pub fn new(atoms: &'a [Atom], atom_container: &'a AtomContainer, simbox: &'a SimBox) -> Self {
        AtomIterator {
            atoms,
            container_iterator: atom_container.iter(),
            simbox,
        }
    }

    /// Filter atoms located inside the specified geometric shape.
    ///
    /// ## Example
    /// Iterating over all atoms of the system
    /// that are located in a sphere around a specific point.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let system = System::from_file("system.gro").unwrap();
    ///
    /// // construct a sphere located at x = 1, y = 2, z = 3 with a radius of 2.5 nm
    /// let sphere = Sphere::new([1.0, 2.0, 3.0].into(), 2.5);
    ///
    /// for atom in system.atoms_iter().filter_geometry(sphere) {
    ///     println!("{:?}", atom);
    /// }
    /// ```
    pub fn filter_geometry(
        self,
        geometry: impl Shape,
    ) -> FilterAtomIterator<'a, AtomIterator<'a>, impl Shape> {
        let simbox = self.simbox;
        FilterAtomIterator {
            iterator: self,
            geometry,
            simbox,
        }
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
    simbox: &'a SimBox,
}

impl<'a, I, S> FilterAtomIterator<'a, I, S>
where
    I: Iterator<Item = &'a Atom>,
    S: Shape,
{
    pub fn filter_geometry(
        self,
        geometry: impl Shape,
    ) -> FilterAtomIterator<'a, FilterAtomIterator<'a, I, S>, impl Shape> {
        let simbox = self.simbox;
        FilterAtomIterator {
            iterator: self,
            geometry,
            simbox,
        }
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
                self.simbox))
    }
}

/// Immutable iterator over atoms of a molecule.
/// Iterates over atoms respecting the topology of the molecule.
pub struct MoleculeIterator<'a> {
    atoms: &'a [Atom],
    container: Vec<usize>,
    current_index: usize,
    simbox: &'a SimBox,
}

impl<'a> MoleculeIterator<'a> {
    /// Create a new `MoleculeIterator`
    ///
    /// ## Parameters
    /// - `atoms`: reference to the atoms of the `System`
    /// - `container`: vector specifying indices of the atoms to iterate through
    /// - `simbox`: current dimensions of the simulation box
    pub fn new(atoms: &'a [Atom], container: Vec<usize>, simbox: &'a SimBox) -> Self {
        MoleculeIterator {
            atoms,
            container,
            current_index: 0,
            simbox,
        }
    }

    /// Filter atoms located inside the specified geometric shape.
    pub fn filter_geometry(
        self,
        geometry: impl Shape,
    ) -> FilterAtomIterator<'a, MoleculeIterator<'a>, impl Shape> {
        let simbox = self.simbox;
        FilterAtomIterator {
            iterator: self,
            geometry,
            simbox,
        }
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
    simbox: &'a SimBox,
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
        simbox: &'a SimBox,
    ) -> Self {
        MutAtomIterator {
            atoms: atoms as *mut [Atom],
            container_iterator: atom_container.iter(),
            simbox,
        }
    }

    /// Filter atoms located inside the specified geometric shape.
    ///
    /// ## Example
    /// Iterating over all atoms of the system
    /// that are located in a sphere around a specific point
    /// and change their names.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// // construct a sphere located at x = 1, y = 2, z = 3 with a radius of 2.5 nm
    /// let sphere = Sphere::new([1.0, 2.0, 3.0].into(), 2.5);
    ///
    /// for atom in system.atoms_iter_mut().filter_geometry(sphere) {
    ///     atom.set_atom_name("XYZ");
    /// }
    /// ```
    pub fn filter_geometry(
        self,
        geometry: impl Shape,
    ) -> MutFilterAtomIterator<'a, MutAtomIterator<'a>, impl Shape> {
        let simbox = self.simbox;
        MutFilterAtomIterator {
            iterator: self,
            geometry,
            simbox,
        }
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

/// Mutable iterator over atoms with applied geometry filter.
/// Constructed by calling `filter_geometry` method on `MutAtomIterator` or `MutFilterAtomIterator`.
/// ## Notes
/// - Atoms with no positions set are never inside any geometric shape.
pub struct MutFilterAtomIterator<'a, I, S>
where
    I: Iterator<Item = &'a mut Atom>,
    S: Shape,
{
    iterator: I,
    geometry: S,
    simbox: &'a SimBox,
}

impl<'a, I, S> MutFilterAtomIterator<'a, I, S>
where
    I: Iterator<Item = &'a mut Atom>,
    S: Shape,
{
    pub fn filter_geometry(
        self,
        geometry: impl Shape,
    ) -> MutFilterAtomIterator<'a, MutFilterAtomIterator<'a, I, S>, impl Shape> {
        let simbox = self.simbox;
        MutFilterAtomIterator {
            iterator: self,
            geometry,
            simbox,
        }
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
                self.simbox))
    }
}

/// Mutable iterator over atoms of a molecule.
/// Iterates over atoms respecting the topology of the molecule.
pub struct MutMoleculeIterator<'a> {
    atoms: *mut [Atom],
    container: Vec<usize>,
    current_index: usize,
    simbox: &'a SimBox,
}

impl<'a> MutMoleculeIterator<'a> {
    /// Create a new `MoleculeIterator`
    ///
    /// ## Parameters
    /// - `atoms`: mutable reference to the atoms of the `System`
    /// - `container`: vector specifying indices of the atoms to iterate through
    /// - `simbox`: current dimensions of the simulation box
    pub fn new(atoms: &'a mut [Atom], container: Vec<usize>, simbox: &'a SimBox) -> Self {
        MutMoleculeIterator {
            atoms,
            container,
            current_index: 0,
            simbox,
        }
    }

    /// Filter atoms located inside the specified geometric shape.
    pub fn filter_geometry(
        self,
        geometry: impl Shape,
    ) -> MutFilterAtomIterator<'a, MutMoleculeIterator<'a>, impl Shape> {
        let simbox = self.simbox;
        MutFilterAtomIterator {
            iterator: self,
            geometry,
            simbox,
        }
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
