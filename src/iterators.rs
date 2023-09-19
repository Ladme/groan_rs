// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of iterators over atoms and filter functions.

use crate::errors::GroupError;
use crate::structures::{atom::Atom, group::Group, shape::Shape, simbox::SimBox};
use crate::system::general::System;

/**************************/
/*  IMMUTABLE ITERATORS   */
/**************************/

/// Immutable iterator over atoms. Constructed using `System::atom_iter()` or `System::group_iter()`.
pub struct AtomIterator<'a> {
    atoms: &'a [Atom],
    atom_ranges: &'a [(usize, usize)],
    current_range_index: usize,
    current_index: usize,
    simbox: &'a SimBox,
}

impl<'a> AtomIterator<'a> {
    pub fn new(atoms: &'a [Atom], atom_ranges: &'a [(usize, usize)], simbox: &'a SimBox) -> Self {
        AtomIterator {
            atoms,
            atom_ranges,
            current_range_index: 0,
            current_index: 0,
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

impl<'a> Iterator for AtomIterator<'a> {
    type Item = &'a Atom;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some((start, end)) = self.atom_ranges.get(self.current_range_index) {
            if self.current_index < *start {
                self.current_index = *start;
            }

            if self.current_index <= *end {
                let atom = unsafe { self.atoms.get_unchecked(self.current_index) };
                self.current_index += 1;
                return Some(atom);
            }

            self.current_range_index += 1;
        }

        None
    }
}

/// Immutable iterator over atoms with applied geometry filter.
/// Constructed by calling `filter_geometry` method on `AtomIterator` or `FilterAtomIterator`.
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
            .find(|atom| self.geometry.inside(atom.get_position(), self.simbox))
    }
}

/**************************/
/*    MUTABLE ITERATORS   */
/**************************/

/// Mutable iterator over atoms. Constructed using `System::atom_iter_mut()` or `System::group_iter_mut()`.
pub struct MutAtomIterator<'a> {
    atoms: *mut [Atom],
    atom_ranges: &'a [(usize, usize)],
    current_range_index: usize,
    current_index: usize,
    simbox: &'a SimBox,
}

impl<'a> MutAtomIterator<'a> {
    pub fn new(
        atoms: &'a mut [Atom],
        atom_ranges: &'a [(usize, usize)],
        simbox: &'a SimBox,
    ) -> Self {
        MutAtomIterator {
            atoms: atoms as *mut [Atom],
            atom_ranges,
            current_range_index: 0,
            current_index: 0,
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
        while let Some((start, end)) = self.atom_ranges.get(self.current_range_index) {
            if self.current_index < *start {
                self.current_index = *start;
            }

            if self.current_index <= *end {
                let atom = unsafe { (*self.atoms).get_unchecked_mut(self.current_index) };
                self.current_index += 1;
                return Some(atom);
            }

            self.current_range_index += 1;
        }

        None
    }
}

/// Mutable iterator over atoms with applied geometry filter.
/// Constructed by calling `filter_geometry` method on `MutAtomIterator` or `MutFilterAtomIterator`.
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
            .find(|atom| self.geometry.inside(atom.get_position(), self.simbox))
    }
}

/**************************/
/*    ITERATING SYSTEM    */
/**************************/

impl System {
    /// Create an iterator over a group of atoms. The atoms are immutable.
    ///
    /// ## Returns
    /// `AtomIterator` or `GroupError::NotFound` in case the group does not exist.
    ///
    /// ## Example
    /// Printing the atoms of group "Protein".
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("Protein", "resid 1 to 29").unwrap();
    ///
    /// match system.group_iter("Protein") {
    ///     Ok(iterator) => {
    ///         for atom in iterator {
    ///             println!("{:?}", atom);
    ///         }
    ///     },
    ///     Err(e) => eprintln!("{}", e),
    /// };
    /// ```
    pub fn group_iter(&self, name: &str) -> Result<AtomIterator, GroupError> {
        let group = self
            .get_groups_as_ref()
            .get(name)
            .ok_or(GroupError::NotFound(name.to_string()))?;

        Ok(AtomIterator::new(
            self.get_atoms_as_ref(),
            group.get_atom_ranges(),
            self.get_box_as_ref(),
        ))
    }

    /// Create an iterator over a group of atoms. The atoms are mutable.
    ///
    /// ## Returns
    /// `MutAtomIterator` or `GroupError::NotFound` in case the group does not exist.
    ///
    /// ## Example
    /// Translating the atoms of the group "Protein".
    /// Note that using `system.group_translate()` may be faster.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("Protein", "resid 1 to 29").unwrap();
    /// let simulation_box = system.get_box_copy();
    ///
    /// match system.group_iter_mut("Protein") {
    ///     Ok(iterator) => {
    ///         for atom in iterator {
    ///             atom.translate(&[1.0, 2.0, -1.0].into(), &simulation_box);
    ///         }
    ///     },
    ///     Err(e) => eprintln!("{}", e),
    /// };
    /// ```
    pub fn group_iter_mut(&mut self, name: &str) -> Result<MutAtomIterator, GroupError> {
        unsafe {
            let simbox = self.get_box_as_ref() as *const SimBox;

            let group =
                self.get_groups_as_ref_mut()
                    .get_mut(name)
                    .ok_or(GroupError::NotFound(name.to_string()))? as *mut Group;

            Ok(MutAtomIterator::new(
                self.get_atoms_as_ref_mut(),
                (*group).get_atom_ranges(),
                simbox.as_ref().unwrap(),
            ))
        }
    }

    /// Create an iterator over all atoms in the system. The atoms are immutable.
    ///
    /// ## Example
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// for atom in system.atoms_iter() {
    ///     println!("{:?}", atom);
    /// }
    /// ```
    ///
    /// ## Note on performance
    /// It might be slightly faster to iterate using `system.get_atoms_as_ref().iter()` if you do
    /// not care about the additional methods `AtomIterator` implements.
    pub fn atoms_iter(&self) -> AtomIterator {
        self.group_iter("all")
            .expect("Groan error. Default group `all` does not exist but it should.")
    }

    /// Create an iterator over all atoms in the system. The atoms are mutable.
    ///
    /// ## Example
    /// Translating all the atoms in the system by a specified vector.
    /// Note that using `system.atoms_translate()` may be faster.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// let mut system = System::from_file("system.gro").unwrap();
    /// let simulation_box = system.get_box_copy();
    ///
    /// for atom in system.atoms_iter_mut() {
    ///     atom.translate(&[1.0, -1.0, 2.5].into(), &simulation_box);
    /// }
    /// ```
    pub fn atoms_iter_mut(&mut self) -> MutAtomIterator {
        self.group_iter_mut("all")
            .expect("Groan error. Default group `all` does not exist but it should.")
    }
}

/**************************/
/*       UNIT TESTS       */
/**************************/

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structures::{dimension::Dimension, shape::*, vector3d::Vector3D};

    use float_cmp::assert_approx_eq;

    #[test]
    fn group_iter_protein() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        for (group_atom, system_atom) in system
            .group_iter("Protein")
            .unwrap()
            .zip(system.get_atoms_as_ref().iter().take(61))
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }
    }

    #[test]
    fn group_iter_membrane() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        for (group_atom, system_atom) in system
            .group_iter("Membrane")
            .unwrap()
            .zip(system.get_atoms_as_ref().iter().skip(61))
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }
    }

    #[test]
    fn group_iter_ion() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        for (group_atom, system_atom) in system
            .group_iter("ION")
            .unwrap()
            .zip(system.get_atoms_as_ref().iter().skip(16604))
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());
        }
    }

    #[test]
    fn atoms_iter() {
        let system = System::from_file("test_files/example.gro").unwrap();

        let atoms = system.get_atoms_as_ref();

        for (extracted_atom, system_atom) in atoms.iter().zip(system.atoms_iter()) {
            assert_eq!(
                system_atom.get_atom_number(),
                extracted_atom.get_atom_number()
            );
        }
    }

    #[test]
    fn group_iter_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let extracted = system.get_atoms_copy();

        for (group_atom, system_atom) in system
            .group_iter_mut("Protein")
            .unwrap()
            .zip(extracted.iter().take(61))
        {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());

            group_atom.translate_nopbc(&Vector3D::from([0.5, -1.1, 2.4]));
            assert_approx_eq!(
                f32,
                group_atom.get_position().x,
                system_atom.get_position().x + 0.5
            );
            assert_approx_eq!(
                f32,
                group_atom.get_position().y,
                system_atom.get_position().y - 1.1
            );
            assert_approx_eq!(
                f32,
                group_atom.get_position().z,
                system_atom.get_position().z + 2.4
            );
        }
    }

    #[test]
    fn atoms_iter_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let extracted = system.get_atoms_copy();

        for (group_atom, system_atom) in system.atoms_iter_mut().zip(extracted.iter()) {
            assert_eq!(system_atom.get_atom_number(), group_atom.get_atom_number());

            group_atom.translate_nopbc(&Vector3D::from([0.5, -1.1, 2.4]));
            assert_approx_eq!(
                f32,
                group_atom.get_position().x,
                system_atom.get_position().x + 0.5
            );
            assert_approx_eq!(
                f32,
                group_atom.get_position().y,
                system_atom.get_position().y - 1.1
            );
            assert_approx_eq!(
                f32,
                group_atom.get_position().z,
                system_atom.get_position().z + 2.4
            );
        }
    }

    #[test]
    fn filter_sphere() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let sphere_pos = system.group_get_center("Protein").unwrap();
        let sphere = Sphere::new(sphere_pos.clone(), 5.0);

        for atom in system.atoms_iter().filter_geometry(sphere.clone()) {
            assert!(
                atom.get_position()
                    .distance(&sphere_pos, Dimension::XYZ, system.get_box_as_ref())
                    < 5.0
            );
        }

        let sphere_pos2 = Vector3D::from([9.0, 3.0, 1.0]);
        let sphere2 = Sphere::new(sphere_pos2.clone(), 4.0);

        for atom in system
            .atoms_iter()
            .filter_geometry(sphere)
            .filter_geometry(sphere2)
        {
            assert!(
                atom.get_position()
                    .distance(&sphere_pos, Dimension::XYZ, system.get_box_as_ref())
                    < 5.0
            );
            assert!(
                atom.get_position()
                    .distance(&sphere_pos2, Dimension::XYZ, system.get_box_as_ref())
                    < 4.0
            );
        }
    }

    #[test]
    fn filter_sphere_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let sphere = Sphere::new(system.group_get_center("Protein").unwrap(), 1.5);
        let sphere2 = Sphere::new(system.group_get_center("Protein").unwrap(), 0.4);

        let sbox = system.get_box_copy();

        let orig_atoms = system.get_atoms_copy();

        for atom in system
            .atoms_iter_mut()
            .filter_geometry(sphere.clone())
            .filter_geometry(sphere2.clone())
        {
            atom.set_atom_name("XYZ");
        }

        for (i, atom) in system.atoms_iter().enumerate() {
            if sphere.inside(atom.get_position(), &sbox)
                && sphere2.inside(atom.get_position(), &sbox)
            {
                assert_eq!(atom.get_atom_name(), "XYZ");
            } else {
                assert_eq!(atom.get_atom_name(), orig_atoms[i].get_atom_name());
            }
        }
    }

    #[test]
    fn filter_sphere_count_water() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let sphere = Sphere::new(system.group_get_center("Protein").unwrap(), 4.5);
        let count = system
            .group_iter("W")
            .unwrap()
            .filter_geometry(sphere)
            .count();

        assert_eq!(count, 1304);
    }

    #[test]
    fn filter_sphere_count_phosphates() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.group_create("Phosphates", "name PO4").unwrap();

        let sphere = Sphere::new(system.group_get_center("Protein").unwrap(), 2.1);
        let count = system
            .group_iter("Phosphates")
            .unwrap()
            .filter_geometry(sphere)
            .count();

        assert_eq!(count, 6);
    }

    #[test]
    fn filter_xcylinder_count_water() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let cylinder = Cylinder::new(
            system.group_get_center("Protein").unwrap(),
            2.0,
            3.0,
            Dimension::X,
        );
        let count = system
            .group_iter("W")
            .unwrap()
            .filter_geometry(cylinder)
            .count();

        assert_eq!(count, 29);
    }

    #[test]
    fn filter_ycylinder_count_water() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let cylinder = Cylinder::new(
            system.group_get_center("Protein").unwrap(),
            2.0,
            3.0,
            Dimension::Y,
        );
        let count = system
            .group_iter("W")
            .unwrap()
            .filter_geometry(cylinder)
            .count();

        assert_eq!(count, 16);
    }

    #[test]
    fn filter_zcylinder_count_water() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let cylinder = Cylinder::new(
            system.group_get_center("Protein").unwrap(),
            2.0,
            3.0,
            Dimension::Z,
        );
        let count = system
            .group_iter("W")
            .unwrap()
            .filter_geometry(cylinder)
            .count();

        assert_eq!(count, 79);
    }

    #[test]
    fn filter_rectangular_count_water() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let rect = Rectangular::new(system.group_get_center("Protein").unwrap(), 2.0, 3.0, 4.0);
        let count = system
            .group_iter("W")
            .unwrap()
            .filter_geometry(rect)
            .count();

        assert_eq!(count, 92);
    }

    #[test]
    fn filter_rectangular_full_count_water() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let rect = Rectangular::new([0.0, 0.0, 0.0].into(), 100.0, 100.0, 100.0);
        let count = system
            .group_iter("W")
            .unwrap()
            .filter_geometry(rect)
            .count();

        assert_eq!(count, system.group_get_n_atoms("W").unwrap());
    }
}
