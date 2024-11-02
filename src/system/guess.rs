// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of System methods for guessing properties of atoms.

use std::fmt;

use colored::Colorize;
use indexmap::IndexMap;

use crate::errors::ElementError;
use crate::structures::dimension::Dimension;
use crate::structures::element::Elements;
use crate::system::System;

/// Default factor used to quess bonds between atoms.
const DEFAULT_RADIUS_FACTOR: f32 = 0.55;

/// ## Methods for guessing properties of the System.
impl System {
    /// Guess elements of the atoms in the system based on the provided `Elements` structure.
    ///
    /// ## Returns
    /// - Returns `Ok` if the assignment was successful, all atoms were assigned an element and no
    ///   atom matched multiple elements.
    /// - Returns `ElementError::ElementGuessWarning` if there is at least one atom which was not assigned
    ///   an element and/or at least one atom which matches multiple elements. This does not
    ///   indicate an error as other atoms are assigned an element.
    /// - Returns a different `ElementError` if the assignment was not successful. In case
    ///   an error occurs, the system is not modified. Errors take precedence over warnings.
    ///
    /// ## Examples
    /// Guess elements for the atoms of the system using **default** element definitions.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.guess_elements(Elements::default());
    /// ```
    /// Guess elements for the atoms of the system using **custom** element definitions.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.guess_elements(Elements::from_file("my_elements.yaml").unwrap());
    /// ```
    ///
    /// The `my_elements.yaml` is a YAML file which should look something like this:
    /// ```yaml
    /// ---
    /// # name of the element
    /// element1:
    ///    # symbol of the element
    ///    symbol: E1
    ///    # groan selection language query identifying the atoms of the element
    ///    query: name r'^[Ee]1.*'
    ///    # mass of the atom of the element
    ///    mass: 16.743
    ///    # van der Waals radius of the atom of the element
    ///    vdw: 0.11
    ///    # expected maximal number of bonds
    ///    expected_max_bonds: 3
    /// element2:
    ///    query: name r'^[Ee]2.*'
    ///    mass: 13.001
    ///    # all properties are optional
    /// ...
    /// ```
    ///
    /// Guess elements for the atoms of the system using default element definitions
    /// with some **minor modifications**.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// let mut elements = Elements::default();
    /// elements.update(Elements::from_file("my_modifications.yaml").unwrap());
    /// system.guess_elements(elements);
    /// ```
    /// The `my_modifications.yaml` is a YAML file containing only the changes
    /// to be made in the default `Elements` structure. For instance, if you want
    /// to change the selection query identifying the carbon atoms, you can use:
    /// ```yaml
    /// ---
    /// carbon:
    ///    query: name r'^[Cc].*$'
    /// ...
    /// ```
    ///
    /// ## Notes
    /// - This function will assign an element to each atom if it finds a matching query in the `Elements`
    ///   structure. If a match is found, the function will set the `element_name` field.
    ///   If the matched element is assigned a symbol, the `element_symbol` field of the atom is also set.
    ///   In case the `element_name` and/or `element_symbol` fields have been previously set for the atom,
    ///   they are overwritten.
    /// - The `mass`, `vdw`, and `expected_max_bonds` are set based on the properties of the
    ///   element **only** if they were not previously set. I.e. if the `mass` of an atom is
    ///   `None`, the `mass` will be assigned based on the mass of the element. If the
    ///   `mass` is already assigned and thus **not** `None`, it is not changed. This goes for all
    ///   three of the aforementioned properties.
    /// - In case you want to set new masses (or other properties) for your atoms based on the newly
    ///   assigned elements, run `reset_mass` (or similar function for other properties) on each atom of
    ///   the system before running `System::guess_elements` or run `System::guess_properties` after
    ///   running `System::guess_elements`.
    /// - This function checks for all potential matches between an atom and the provided elements.
    ///   In case there are multiple matches, the function assigns the first element from the list
    ///   and returns a warning.
    /// - Unless an error is returned, the time complexity of this function is _always_ O(m * n + n),
    ///   where `m` is the number of atoms in the system and `n` is the number of available elements.
    /// - If an error (not a warning!) is returned, the `System` structure is not modified.
    pub fn guess_elements(&mut self, elements: Elements) -> Result<(), ElementError> {
        self.guess_elements_partial(elements, true)
    }

    /// Guess elements of the atoms in the system which do not have an assigned element
    /// based on the provided `Elements` structure.
    ///
    /// Behaves similarly to [`System::guess_elements`] but only assigns elements to atoms
    /// that do not have assigned elements already.
    /// No properties are guessed for atoms which already have elements assigned.
    pub fn guess_elements_unknown(&mut self, elements: Elements) -> Result<(), ElementError> {
        self.guess_elements_partial(elements, false)
    }

    /// Guess elements either for all atoms (for_all == true) or only for atoms which are not assigned
    /// an element (for_all == false).
    fn guess_elements_partial(
        &mut self,
        elements: Elements,
        for_all: bool,
    ) -> Result<(), ElementError> {
        // check that all select trees are valid
        // this has to be done before modifying the system because
        // an error should not cause the system to be in an undefined state
        self.validate_queries(&elements)?;

        // note that we do not explicitly expand the regular expression group names in the select trees
        // as these can be expected to be rare in element definitions
        // so they will be expanded individually for each atom

        let mut no_elements = Vec::new();
        let mut multiple_elements: IndexMap<Vec<String>, Vec<usize>> = IndexMap::new();

        for a in 0..self.get_n_atoms() {
            if !for_all {
                let atom = self.get_atom(a).unwrap();
                // if the atom is already assigned an element, skip it
                if atom.get_element_name().is_some() || atom.get_element_symbol().is_some() {
                    continue;
                }
            }

            // get all elements matching the atom
            let matched_elements =
                elements
                    .elements
                    .iter()
                    .try_fold(Vec::new(), |mut acc, (name, element)| {
                        if element.matches(a, self).expect(
                        "FATAL GROAN ERROR | System::guess_elements | Query should not be invalid.",
                    ) {
                        acc.push(name.to_string());
                    }
                        Ok(acc)
                    })?;

            if matched_elements.is_empty() {
                no_elements.push(a + 1);
            } else {
                self.set_atom_properties(a, &matched_elements[0], &elements);

                if matched_elements.len() > 1 {
                    multiple_elements
                        .entry(matched_elements)
                        .or_default()
                        .push(a + 1);
                }
            }
        }

        if !no_elements.is_empty() || !multiple_elements.is_empty() {
            Err(ElementError::ElementGuessWarning(Box::new(
                ElementGuessInfo {
                    no_elements,
                    multiple_elements,
                },
            )))
        } else {
            Ok(())
        }
    }

    /// Guess properties of the atoms of the system from their assigned elements.
    /// The guessed properties include mass, van der Waals radius,
    /// and the expected minimal and maximal number of bonds.
    /// The elements assigned to the atoms are not changed.
    ///
    /// ## Returns
    /// - `Ok` if the assignment was fully successful and all properties have been assigned.
    /// - `ElementError::PropertiesGuessWarning` if at least one of the atoms does not have all the properties assigned
    ///   or if any of the atoms does not have an element defined in the `Elements` structure.
    ///   This does not indicate failure of the function. In fact, this is something that will usually happen,
    ///   and is usually not an issue. Nonetheless, the user should be informed about it.
    ///
    /// ## Examples
    /// Assign elements using default `groan_rs` element definitions and
    /// then assign custom properties to the atoms using your own element definitions
    /// provided in `my_elements.yaml` file.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use groan_rs::errors::ElementError;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.guess_elements(Elements::default()).unwrap();
    ///
    /// match system.guess_properties(Elements::from_file("my_elements.yaml").unwrap()) {
    ///     // ignore warnings
    ///     Ok(_) | Err(ElementError::PropertiesGuessWarning(_)) => (),
    ///     Err(e) => eprintln!("{}", e),
    /// }
    /// ```
    /// The `my_elements.yaml` could look like this:
    /// ```yaml
    /// ---
    /// carbon:
    ///    # custom mass and vdw radius for carbon
    ///    mass: 12.0
    ///    vdw: 0.15
    /// nitrogen:
    ///    # custom expected maximal number of bonds for nitrogen
    ///    expected_max_bonds: 5
    /// ...
    /// ```
    ///
    /// ## Notes
    /// - The properties are guessed from the provided `Elements` structure. However, the elements
    ///   are **not assigned** from the `Elements` structure. This function uses information about the elements
    ///   already stored in the atoms.
    /// - If you already called `System::guess_elements`, you generally do not need to call `System::guess_properties`
    ///   as the properties have already been set to the atoms. You would call `System::guess_properties`
    ///   only if you want to force overwrite of the previously assigned properties of the atoms.
    /// - If you have assigned the elements to the atoms in some other way, i.e. without calling
    ///   `System::guess_elements`, use this function to guess the properties of atoms from the assigned elements.
    /// - If the element does not have a specific property assigned, the property of the atom is not changed.
    ///   E.g. if the element has no `mass` set, but the atom corresponding to the element does,
    ///   the mass of the atom is **not** reset.
    /// - When using default element definitions provided by the `groan_rs`, this function will usually
    ///   returns a warning as most elements do not have all properties properly set as they are usually not needed.
    ///   This is intentional and does not indicate failure of the function.
    ///   It only informs the user about potential issues they might not expect.
    pub fn guess_properties(&mut self, elements: Elements) -> Result<(), ElementError> {
        let mut info = PropertiesGuessInfo {
            no_element: Vec::new(),
            not_recognized: Vec::new(),
            no_mass: Vec::new(),
            no_vdw: Vec::new(),
            no_max_bonds: Vec::new(),
            no_min_bonds: Vec::new(),
        };

        for (a, atom) in self.atoms_iter_mut().enumerate() {
            if let Some(elname) = atom.get_element_name() {
                if let Some(element) = elements.elements.get(elname) {
                    match element.mass {
                        None => info.no_mass.push(a + 1),
                        Some(x) => atom.set_mass(x),
                    }

                    match element.vdw {
                        None => info.no_vdw.push(a + 1),
                        Some(x) => atom.set_vdw(x),
                    }

                    match element.expected_max_bonds {
                        None => info.no_max_bonds.push(a + 1),
                        Some(x) => atom.set_expected_max_bonds(x),
                    }

                    match element.expected_min_bonds {
                        None => info.no_min_bonds.push(a + 1),
                        Some(x) => atom.set_expected_min_bonds(x),
                    }
                } else {
                    info.not_recognized.push(a + 1);
                }
            } else {
                info.no_element.push(a + 1);
            }
        }

        if info.is_empty() {
            Ok(())
        } else {
            Err(ElementError::PropertiesGuessWarning(Box::new(info)))
        }
    }

    /// Assign bonds between atoms based on the distances between them, their van der Waals radii
    /// and the provided `radius_factor`.
    /// If the `radius_factor` is not provided, the default value of 0.55 is used.
    ///
    ///
    /// ## Returns
    /// - `Ok` if the bonds were guessed, all atoms have vdw information and no atom has suspicious number of bonds.
    /// - `ElementError::BondGuessWarning` if at least one atom has a suspicious number of
    ///   guessed bonds or if at least one atom is missing vdw radius.
    ///   This does not indicate failure of the function.
    ///
    /// ## Warning
    /// - Currently only works with orthogonal periodic boundary conditions!
    /// - This function is rather unreliable. If you require precise information
    ///   about the topology of your system, load in a TPR file, a PDB file with CONECT information,
    ///   or define the topology manually.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use groan_rs::errors::ElementError;
    /// #
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// match system.guess_elements(Elements::default()) {
    ///     Ok(_) | Err(ElementError::ElementGuessWarning(_)) => (),
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;    
    ///     }
    /// }
    ///
    /// match system.guess_bonds(None) {
    ///     Ok(_) => (),
    ///     Err(ElementError::BondsGuessWarning(e)) => eprintln!("{}", e),
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;
    ///     }
    /// }
    /// ```
    ///
    /// ## How does it work?
    /// A distance between each pair of atoms is calculated. If the distance is
    /// lower than the sum of vdw radii of the atoms multiplied by the provided
    /// `radius_factor`, a bond between the atoms is formed.
    ///
    /// If the `radius_factor` is `None`, the value defined in `DEFAULT_RADIUS_FACTOR`
    /// is used instead (default value is 0.55).
    ///
    /// The function then checks the number of bonds for each atom. If the number of bonds
    /// for any atom is lower than the expected minimal number of bonds, a warning is raised.
    /// Warning is also raised if the number of bonds for any atom is higher than the expected
    /// maximal number of bonds. Warning is also raised if any of the atoms in the system
    /// does not have a van der Waals radius assigned.
    ///
    /// ## Notes
    /// - Atoms which have no assigned van der Waals radius are not assigned any bonds.
    /// - Asymptotic time complexity of this function is O(n^2) where n is the number of atoms in the system.
    /// - It is almost always useful to use the parallelized version of this function: [`System::guess_bonds_parallel`],
    ///   especially if your system is large. To use this function, you have to use the `parallel` feature of `groan_rs`.
    pub fn guess_bonds(&mut self, radius_factor: Option<f32>) -> Result<(), ElementError> {
        let n_atoms = self.get_n_atoms();
        if n_atoms == 0 {
            return Ok(());
        }

        // identify bonds
        let (bonds, no_vdw) = self.identify_bonds(
            radius_factor.unwrap_or(DEFAULT_RADIUS_FACTOR),
            0,
            self.get_n_atoms(),
        );

        // assign bonds
        self.assign_bonds(bonds);

        // check for unexpected number of bonds
        let (too_many_bonds, too_few_bonds) = self.check_unexpected_bonds();

        let info = BondsGuessInfo {
            no_vdw,
            too_many_bonds,
            too_few_bonds,
        };

        if info.is_empty() {
            Ok(())
        } else {
            Err(ElementError::BondsGuessWarning(Box::new(info)))
        }
    }

    /// Assign bonds between atoms based on the distances between them, their van der Waals radii
    /// and the provided `radius_factor`.
    /// If the `radius_factor` is not provided, the default value of 0.55 is used.
    ///
    /// This function works the same as [`System::guess_bonds`] but can employ multiple threads
    /// significantly increasing the speed of the calculation. The number of threads (`n_threads`)
    /// to use is provided by the user.
    ///
    /// ## Panics
    /// Panics if the number of threads (`n_threads`) to be employed equals zero.
    ///
    /// ## Notes
    /// - If the number of threads is higher than the number of atoms (`n_atoms`) in the system,
    ///   only `n_atoms` threads will be used.
    #[cfg(any(feature = "parallel", doc))]
    pub fn guess_bonds_parallel(
        &mut self,
        radius_factor: Option<f32>,
        mut n_threads: usize,
    ) -> Result<(), ElementError> {
        if n_threads == 0 {
            panic!("FATAL GROAN ERROR | System::guess_bonds_parallel | Number of threads should not be zero.");
        }

        let n_atoms = self.get_n_atoms();
        if n_atoms == 0 {
            return Ok(());
        }

        // if the number of atoms is higher than the number of threads, decrease the number of threads
        if n_atoms < n_threads {
            n_threads = n_atoms;
        }

        // identify bonds in parallel
        let (bonds, no_vdw) = self.identify_bonds_parallel(
            radius_factor.unwrap_or(DEFAULT_RADIUS_FACTOR),
            n_atoms,
            n_threads,
        );

        // assign bonds
        self.assign_bonds(bonds);

        // check for unexpected number of bonds
        let (too_many_bonds, too_few_bonds) = self.check_unexpected_bonds();

        let info = BondsGuessInfo {
            no_vdw,
            too_many_bonds,
            too_few_bonds,
        };

        if info.is_empty() {
            Ok(())
        } else {
            Err(ElementError::BondsGuessWarning(Box::new(info)))
        }
    }

    /// Assign bonds to atoms of the system.
    ///
    /// `bonds` must only contain valid indices.
    #[inline]
    fn assign_bonds(&mut self, bonds: Vec<(usize, usize)>) {
        // safety: only safe if bonds vector contains valid indices
        // we apply `System::reset_mol_references` at the end of the function
        unsafe {
            let atoms = self.get_atoms_mut();
            for (index1, index2) in bonds {
                atoms[index1].add_bonded(index2);
                atoms[index2].add_bonded(index1);
            }
        }

        self.reset_mol_references();
    }

    /// Assign bonds to the atoms of the system.
    ///
    /// ## Warning
    /// - `n_threads` must be non-zero and smaller than or equal to `n_atoms`
    /// - `n_atoms` must be non-zero
    #[cfg(feature = "parallel")]
    #[inline]
    fn identify_bonds_parallel(
        &self,
        radius_factor: f32,
        n_atoms: usize,
        n_threads: usize,
    ) -> (Vec<(usize, usize)>, Vec<usize>) {
        // distribute atoms between threads
        let distribution = self.distribute_atoms(n_atoms, n_threads);

        let mut no_vdw = Vec::new();
        let mut bonds = Vec::new();

        std::thread::scope(|s| {
            let mut handles = Vec::new();

            for (start, end) in distribution {
                let handle = s.spawn(move || -> (Vec<(usize, usize)>, Vec<usize>) {
                    self.identify_bonds(radius_factor, start, end)
                });

                handles.push(handle);
            }

            for handle in handles {
                let (mut b, mut vdw) = handle.join().expect(
                    "FATAL GROAN ERROR | System::identify_bonds_parallel | Could not join handle.",
                );
                bonds.append(&mut b);
                no_vdw.append(&mut vdw);
            }
        });

        (bonds, no_vdw)
    }

    /// Identify atoms which should be bonded.
    ///
    /// ## Parameters
    /// - `radius_factor`: factor used to quess bonds between atoms.
    /// - `start` and `end`: range of atom indices for which the bonds should be assigned
    #[inline]
    fn identify_bonds(
        &self,
        radius_factor: f32,
        start: usize,
        end: usize,
    ) -> (Vec<(usize, usize)>, Vec<usize>) {
        let mut no_vdw = Vec::new();
        let mut bonds = Vec::new();

        for a in start..end {
            let atom1 = self
                .get_atom(a)
                .expect("FATAL GROAN ERROR | System::identify_bonds | Atom 1 should exist.");

            let vdw1 = match atom1.get_vdw() {
                Some(x) => x,
                None => {
                    no_vdw.push(a + 1);
                    continue;
                }
            };

            for b in (a + 1)..self.get_n_atoms() {
                let atom2 = self
                    .get_atom(b)
                    .expect("FATAL GROAN ERROR | System::identify_bonds | Atom 2 should exist.");

                let vdw2 = match atom2.get_vdw() {
                    Some(x) => x,
                    None => continue,
                };

                let distance = self
                    .atoms_distance(a, b, Dimension::XYZ)
                    .expect("FATAL GROAN ERROR | System::identify_bonds | Atoms should exist.");
                let limit = (vdw1 + vdw2) * radius_factor;

                if distance < limit {
                    bonds.push((a, b));
                }
            }
        }

        (bonds, no_vdw)
    }

    /// Check for unexpected number of bonds in a system.
    #[inline]
    #[allow(clippy::type_complexity)]
    fn check_unexpected_bonds(
        &self,
    ) -> (IndexMap<usize, (usize, u8)>, IndexMap<usize, (usize, u8)>) {
        let mut too_many_bonds = IndexMap::new();
        let mut too_few_bonds = IndexMap::new();

        for atom in self.atoms_iter() {
            // check limit for maximal number of bonds
            if let Some(limit) = atom.get_expected_max_bonds() {
                if atom.get_n_bonded() > limit as usize
                    && too_many_bonds
                        .insert(atom.get_index() + 1, (atom.get_n_bonded(), limit))
                        .is_some()
                {
                    panic!("FATAL GROAN ERROR | System::guess_bonds | Atom should not be in the `too_many_bonds` map.")
                }
            }

            // check limit for minimal number of bonds
            if let Some(limit) = atom.get_expected_min_bonds() {
                if atom.get_n_bonded() < limit as usize
                    && too_few_bonds
                        .insert(atom.get_index() + 1, (atom.get_n_bonded(), limit))
                        .is_some()
                {
                    panic!("FATAL GROAN ERROR | System::guess_bonds | Atom should not be in the `too_few_bonds` map.")
                }
            }
        }

        (too_many_bonds, too_few_bonds)
    }

    /// Set properties of the atom based on properties of the element.
    /// Panics if the atom or the element does not exist.
    fn set_atom_properties(&mut self, atom_index: usize, element_name: &str, elements: &Elements) {
        let atom = self
            .get_atom_mut(atom_index)
            .expect("FATAL GROAN ERROR | System::set_atom_properties | Atom should exist.");

        let element = elements
            .elements
            .get(element_name)
            .expect("FATAL GROAN ERROR | System::set_atom_properties | Element should exist.");

        atom.set_element_name(element_name);

        if let Some(ref symbol) = element.symbol {
            atom.set_element_symbol(symbol);
        }

        // set mass, vdw, expected bonds if not already set
        if let Some(x) = atom.get_mass().or(element.mass) {
            atom.set_mass(x)
        };
        if let Some(x) = atom.get_vdw().or(element.vdw) {
            atom.set_vdw(x)
        };
        if let Some(x) = atom.get_expected_max_bonds().or(element.expected_max_bonds) {
            atom.set_expected_max_bonds(x)
        };
        if let Some(x) = atom.get_expected_min_bonds().or(element.expected_min_bonds) {
            atom.set_expected_min_bonds(x)
        }
    }

    /// Check validity of all select trees in the elements structure.
    /// Returns `Ok` if the queries are valid, otherwise returns an `ElementError`.
    fn validate_queries(&self, elements: &Elements) -> Result<(), ElementError> {
        for (_, element) in &elements.elements {
            element.validate_select(self)?;
        }

        Ok(())
    }
}

/******************************/
/*  ADVANCED ERROR MESSAGES   */
/******************************/

/// Contains information about issues raised while guessing elements in a system.
#[derive(Debug, PartialEq, Eq)]
pub struct ElementGuessInfo {
    /// Atoms which were not assigned an element.
    no_elements: Vec<usize>,
    /// Atoms which were assigned multiple elements.
    multiple_elements: IndexMap<Vec<String>, Vec<usize>>,
}

impl fmt::Display for ElementGuessInfo {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let no_elements_len = self.no_elements.len();
        if no_elements_len > 0 {
            writeln!(f, "\n> following atoms match no element:")?;
        }

        for (a, atom) in self.no_elements.iter().enumerate() {
            write!(f, "{} ", atom.to_string().yellow())?;
            if a >= 9 && no_elements_len != a + 1 {
                writeln!(f, "and {} other atoms", no_elements_len - a - 1)?;
                break;
            }
        }

        for (matches, atoms) in self.multiple_elements.iter() {
            let assigned = matches
                .first()
                .expect("FATAL GROAN ERROR | ElementGuessInfo::fmt | Element should exist.");
            write!(
                f,
                "\n> following atoms have been identified as {} but also match",
                assigned.to_string().yellow()
            )?;
            for (e, element) in matches.iter().skip(1).enumerate() {
                write!(f, " {}", element.yellow())?;

                if e >= 2 && matches.len() != e + 2 {
                    writeln!(f, " and {} other elements:", matches.len() - e - 2)?;
                    break;
                }
            }

            if matches.len() <= 4 {
                writeln!(f, ":")?;
            }

            for (a, atom) in atoms.iter().enumerate() {
                write!(f, "{} ", atom.to_string().yellow())?;

                if a >= 9 && atoms.len() != a + 1 {
                    writeln!(f, "and {} other atoms ", atoms.len() - a - 1)?;
                    break;
                }
            }

            if atoms.len() <= 10 {
                writeln!(f)?;
            }
        }

        Ok(())
    }
}

/// Contains information about issues raised while guessing properties of atoms in a system.
#[derive(Debug, PartialEq, Eq)]
pub struct PropertiesGuessInfo {
    /// Atoms which do not have an element assigned.
    no_element: Vec<usize>,
    /// Atoms which do have an element assigned but it can't be recognized.
    not_recognized: Vec<usize>,
    /// Atoms which are missing mass information.
    no_mass: Vec<usize>,
    /// Atoms which are missing vdw information.
    no_vdw: Vec<usize>,
    /// Atoms which are missing expected_max_bonds information.
    no_max_bonds: Vec<usize>,
    /// Atoms which are missing expected_min_bonds information.
    no_min_bonds: Vec<usize>,
}

impl PropertiesGuessInfo {
    /// Return `true` if all of the fields of `PropertiesGuessInfo` are empty.
    fn is_empty(&self) -> bool {
        self.no_element.is_empty()
            && self.no_mass.is_empty()
            && self.no_vdw.is_empty()
            && self.no_max_bonds.is_empty()
            && self.no_min_bonds.is_empty()
    }
}

impl fmt::Display for PropertiesGuessInfo {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let fields = [
            &self.not_recognized,
            &self.no_element,
            &self.no_mass,
            &self.no_vdw,
            &self.no_max_bonds,
            &self.no_min_bonds,
        ];
        let strings = [
            "\n> following atoms are assigned an unknown element - no properties were set for them:",
            "\n> following atoms are not assigned an element - no properties were set for them:",
            "\n> following atoms are assigned an element which does not contain information about mass:",
            "\n> following atoms are assigned an element which does not contain information about van der Waals radius:",
            "\n> following atoms are assigned an element which does not contain information about expected maximal number of bonds:",
            "\n> following atoms are assigned an element which does not contain information about expected minimal number of bonds:",
        ];

        for (field, string) in fields.iter().zip(strings.iter()) {
            if !field.is_empty() {
                writeln!(f, "{}", string)?;
            } else {
                continue;
            }

            for (a, atom) in field.iter().enumerate() {
                write!(f, "{} ", atom.to_string().yellow())?;
                if a >= 9 && field.len() != a + 1 {
                    writeln!(f, "and {} other atoms", field.len() - a - 1)?;
                    break;
                }
            }

            if field.len() <= 10 {
                writeln!(f)?;
            }
        }

        Ok(())
    }
}

/// Contains information about issues raised while guessing bonds in a system.
#[derive(Debug, PartialEq, Eq)]
pub struct BondsGuessInfo {
    /// Atoms which have no vdw information assigned.
    no_vdw: Vec<usize>,
    /// Atoms which have too many bonds assigned.
    too_many_bonds: IndexMap<usize, (usize, u8)>,
    /// Atoms which have too few bonds assigned.
    too_few_bonds: IndexMap<usize, (usize, u8)>,
}

impl BondsGuessInfo {
    /// Return `true` if all of the fields of `BondsGuessInfo` are empty.
    fn is_empty(&self) -> bool {
        self.no_vdw.is_empty() && self.too_many_bonds.is_empty() && self.too_few_bonds.is_empty()
    }
}

impl fmt::Display for BondsGuessInfo {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if !self.no_vdw.is_empty() {
            writeln!(f, "\n> following atoms do not have a van der Waals radius assigned - no bonds could be guessed for them:")?;
        }

        for (a, atom) in self.no_vdw.iter().enumerate() {
            write!(f, "{} ", atom.to_string().yellow())?;
            if a >= 9 && self.no_vdw.len() != a + 1 {
                writeln!(f, "and {} other atoms", self.no_vdw.len() - a - 1)?;
                break;
            }
        }

        if self.no_vdw.len() <= 10 {
            writeln!(f)?;
        }

        let fields = [&self.too_many_bonds, &self.too_few_bonds];
        let strings = [
            "\n> following atoms have a suspiciously high number of bonds:",
            "\n> following atoms have a suspiciously low number of bonds:",
        ];

        for (i, (field, string)) in fields.iter().zip(strings.iter()).enumerate() {
            if !field.is_empty() {
                writeln!(f, "{}", string)?;
            } else {
                continue;
            }

            for (a, (atom, (guessed, expected))) in field.iter().enumerate() {
                if i == 0 {
                    writeln!(
                        f,
                        "{} (guessed {}, expected at most {})",
                        atom.to_string().yellow(),
                        guessed.to_string().yellow(),
                        expected.to_string().yellow()
                    )?;
                } else {
                    writeln!(
                        f,
                        "{} (guessed {}, expected at least {})",
                        atom.to_string().yellow(),
                        guessed.to_string().yellow(),
                        expected.to_string().yellow()
                    )?;
                }

                if a >= 9 && field.len() != a + 1 {
                    writeln!(f, "and {} other atoms", field.len() - a - 1)?;
                    break;
                }
            }
        }

        Ok(())
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use crate::{errors::SelectError, structures::simbox::SimBox};

    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn guess_elements() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();

        system.guess_elements(Elements::default()).unwrap();

        for atom in system.atoms_iter() {
            assert!(atom.get_element_name().is_some());
            assert!(atom.get_element_symbol().is_some());
            assert!(atom.get_mass().is_some());
        }

        // N in SER1
        let atom = system.get_atom(0).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "nitrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "N");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 14.0067);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.1625);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 4);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // H1 in SER1
        let atom = system.get_atom(1).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "hydrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "H");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 1.0079);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.10);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 1);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // C in SER23
        let atom = system.get_atom(360).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "carbon");
        assert_eq!(atom.get_element_symbol().unwrap(), "C");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 12.0107);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.17);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 4);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 2);

        // O31 in POPC44
        let atom = system.get_atom(3081).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "oxygen");
        assert_eq!(atom.get_element_symbol().unwrap(), "O");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 15.9994);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.15);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 2);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // P in POPC127
        let atom = system.get_atom(14184).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "phosphorus");
        assert_eq!(atom.get_element_symbol().unwrap(), "P");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 30.9738);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.1871);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 5);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 2);

        // HW1 in SOL4827
        let atom = system.get_atom(31541).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "hydrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "H");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 1.0079);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.10);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 1);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // OW in SOL177
        let atom = system.get_atom(17590).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "oxygen");
        assert_eq!(atom.get_element_symbol().unwrap(), "O");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 15.9994);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.15);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 2);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // NA in NA5250
        let atom = system.get_atom(32795).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "sodium");
        assert_eq!(atom.get_element_symbol().unwrap(), "Na");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 22.9897);
        assert_eq!(atom.get_vdw(), None);
        assert_eq!(atom.get_expected_max_bonds(), None);
        assert_eq!(atom.get_expected_min_bonds(), None);

        // CL in CL5271
        let atom = system.get_atom(32816).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "chlorine");
        assert_eq!(atom.get_element_symbol().unwrap(), "Cl");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 35.453);
        assert_eq!(atom.get_vdw(), None);
        assert_eq!(atom.get_expected_max_bonds(), None);
        assert_eq!(atom.get_expected_min_bonds(), None);
    }

    #[test]
    fn guess_elements_prefilled() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();

        system.get_atom_mut(0).unwrap().set_mass(19.1);
        system.get_atom_mut(0).unwrap().set_element_symbol("Uk");
        system.get_atom_mut(0).unwrap().set_vdw(0.24);
        system.get_atom_mut(360).unwrap().set_expected_max_bonds(7);
        system.get_atom_mut(14184).unwrap().set_vdw(0.20);
        system.get_atom_mut(32795).unwrap().set_mass(19.1);
        system
            .get_atom_mut(32795)
            .unwrap()
            .set_element_name("Unknown");

        system.guess_elements(Elements::default()).unwrap();

        for atom in system.atoms_iter() {
            assert!(atom.get_element_name().is_some());
            assert!(atom.get_element_symbol().is_some());
            assert!(atom.get_mass().is_some());
        }

        // N in SER1
        let atom = system.get_atom(0).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "nitrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "N");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 19.1);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.24);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 4);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // H1 in SER1
        let atom = system.get_atom(1).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "hydrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "H");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 1.0079);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.1);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 1);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // C in SER23
        let atom = system.get_atom(360).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "carbon");
        assert_eq!(atom.get_element_symbol().unwrap(), "C");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 12.0107);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.17);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 7);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 2);

        // O31 in POPC44
        let atom = system.get_atom(3081).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "oxygen");
        assert_eq!(atom.get_element_symbol().unwrap(), "O");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 15.9994);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.15);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 2);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // P in POPC127
        let atom = system.get_atom(14184).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "phosphorus");
        assert_eq!(atom.get_element_symbol().unwrap(), "P");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 30.9738);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.20);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 5);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 2);

        // HW1 in SOL4827
        let atom = system.get_atom(31541).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "hydrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "H");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 1.0079);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.1);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 1);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // OW in SOL177
        let atom = system.get_atom(17590).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "oxygen");
        assert_eq!(atom.get_element_symbol().unwrap(), "O");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 15.9994);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.15);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 2);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // NA in NA5250
        let atom = system.get_atom(32795).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "sodium");
        assert_eq!(atom.get_element_symbol().unwrap(), "Na");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 19.1);
        assert_eq!(atom.get_vdw(), None);
        assert_eq!(atom.get_expected_max_bonds(), None);
        assert_eq!(atom.get_expected_min_bonds(), None);

        // CL in CL5271
        let atom = system.get_atom(32816).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "chlorine");
        assert_eq!(atom.get_element_symbol().unwrap(), "Cl");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 35.453);
        assert_eq!(atom.get_vdw(), None);
        assert_eq!(atom.get_expected_max_bonds(), None);
        assert_eq!(atom.get_expected_min_bonds(), None);
    }

    #[test]
    fn guess_elements_unknown() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();

        system.get_atom_mut(0).unwrap().set_mass(19.1);
        system.get_atom_mut(0).unwrap().set_element_symbol("Uk");
        system.get_atom_mut(0).unwrap().set_vdw(0.24);
        system.get_atom_mut(360).unwrap().set_expected_max_bonds(7);
        system.get_atom_mut(14184).unwrap().set_vdw(0.20);
        system.get_atom_mut(32795).unwrap().set_mass(19.1);
        system
            .get_atom_mut(32795)
            .unwrap()
            .set_element_name("Unknown");

        system.guess_elements_unknown(Elements::default()).unwrap();

        // N in SER1
        let atom = system.get_atom(0).unwrap();
        println!("{:?}", atom);
        assert!(atom.get_element_name().is_none());
        assert_eq!(atom.get_element_symbol().unwrap(), "Uk");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 19.1);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.24);
        assert!(atom.get_expected_max_bonds().is_none());
        assert!(atom.get_expected_min_bonds().is_none());

        // H1 in SER1
        let atom = system.get_atom(1).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "hydrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "H");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 1.0079);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.1);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 1);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // C in SER23
        let atom = system.get_atom(360).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "carbon");
        assert_eq!(atom.get_element_symbol().unwrap(), "C");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 12.0107);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.17);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 7);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 2);

        // O31 in POPC44
        let atom = system.get_atom(3081).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "oxygen");
        assert_eq!(atom.get_element_symbol().unwrap(), "O");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 15.9994);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.15);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 2);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // P in POPC127
        let atom = system.get_atom(14184).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "phosphorus");
        assert_eq!(atom.get_element_symbol().unwrap(), "P");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 30.9738);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.20);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 5);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 2);

        // HW1 in SOL4827
        let atom = system.get_atom(31541).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "hydrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "H");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 1.0079);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.1);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 1);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // OW in SOL177
        let atom = system.get_atom(17590).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "oxygen");
        assert_eq!(atom.get_element_symbol().unwrap(), "O");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 15.9994);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.15);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 2);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // NA in NA5250
        let atom = system.get_atom(32795).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "Unknown");
        assert!(atom.get_element_symbol().is_none());
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 19.1);
        assert!(atom.get_vdw().is_none());
        assert!(atom.get_expected_max_bonds().is_none());
        assert!(atom.get_expected_min_bonds().is_none());

        // CL in CL5271
        let atom = system.get_atom(32816).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "chlorine");
        assert_eq!(atom.get_element_symbol().unwrap(), "Cl");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 35.453);
        assert_eq!(atom.get_vdw(), None);
        assert_eq!(atom.get_expected_max_bonds(), None);
        assert_eq!(atom.get_expected_min_bonds(), None);
    }

    #[test]
    fn guess_elements_with_warnings() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();

        let expected_no: Vec<usize> = vec![
            383, 517, 651, 785, 919, 1053, 1187, 1321, 1455, 1589, 1723, 1857, 1991, 2125, 2259,
            2393, 2527, 2661, 2795, 2929, 3063, 3197, 3331, 3465, 3599, 3733, 3867, 4001, 4135,
            4269, 4403, 4537, 4671, 4805, 4939, 5073, 5207, 5341, 5475, 5609, 5743, 5877, 6011,
            6145, 6279, 6413, 6547, 6681, 6815, 6949, 7083, 7217, 7351, 7485, 7619, 7753, 7887,
            8021, 8155, 8289, 8423, 8557, 8691, 8825, 8959, 9093, 9227, 9361, 9495, 9629, 9763,
            9897, 10031, 10165, 10299, 10433, 10567, 10701, 10835, 10969, 11103, 11237, 11371,
            11505, 11639, 11773, 11907, 12041, 12175, 12309, 12443, 12577, 12711, 12845, 12979,
            13113, 13247, 13381, 13515, 13649, 13783, 13917, 14051, 14185, 14319, 14453, 14587,
            14721, 14855, 14989, 15123, 15257, 15391, 15525, 15659, 15793, 15927, 16061, 16195,
            16329, 16463, 16597, 16731, 16865, 16999, 17133, 17267, 17401,
        ];

        let expected_multiple1: Vec<usize> = vec![
            32803, 32808, 32809, 32810, 32811, 32812, 32813, 32814, 32815, 32816, 32817,
        ];
        let expected_multiple2: Vec<usize> = vec![32804, 32805, 32806, 32807];

        match system
            .guess_elements(Elements::from_file("test_files/elements_incomplete.yaml").unwrap())
        {
            Ok(_) => panic!("Function should have failed."),
            Err(ElementError::ElementGuessWarning(x)) => {
                assert_eq!(x.no_elements, expected_no);
                let atoms = x
                    .multiple_elements
                    .get(&vec!["carbon".to_string(), "chlorine".to_string()])
                    .unwrap();
                assert_eq!(*atoms, expected_multiple1);
                let atoms = x
                    .multiple_elements
                    .get(&vec![
                        "carbon".to_string(),
                        "chlorine".to_string(),
                        "unknown".to_string(),
                    ])
                    .unwrap();
                assert_eq!(*atoms, expected_multiple2);
            }
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }

        for (a, atom) in system.atoms_iter().enumerate() {
            if expected_no.contains(&(a + 1)) {
                assert!(atom.get_element_name().is_none());
                assert!(atom.get_element_symbol().is_none());
                assert!(atom.get_mass().is_none());
            } else {
                assert!(atom.get_element_name().is_some());
                assert!(atom.get_element_symbol().is_some());
                assert!(atom.get_mass().is_some());
            }
        }

        // N in SER1
        let atom = system.get_atom(0).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "nitrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "N");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 14.0067);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.155);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 3);

        // H1 in SER1
        let atom = system.get_atom(1).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "hydrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "H");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 1.0079);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.12);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 1);

        // C in SER23
        let atom = system.get_atom(360).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "carbon");
        assert_eq!(atom.get_element_symbol().unwrap(), "C");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 12.0107);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.17);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 4);

        // O31 in POPC44
        let atom = system.get_atom(3081).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "oxygen");
        assert_eq!(atom.get_element_symbol().unwrap(), "O");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 15.9994);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.152);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 2);

        // P in POPC127
        let atom = system.get_atom(14184).unwrap();
        assert_eq!(atom.get_element_name(), None);
        assert_eq!(atom.get_element_symbol(), None);
        assert_eq!(atom.get_mass(), None);
        assert_eq!(atom.get_vdw(), None);
        assert_eq!(atom.get_expected_max_bonds(), None);

        // HW1 in SOL4827
        let atom = system.get_atom(31541).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "hydrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "H");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 1.0079);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.12);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 1);

        // OW in SOL177
        let atom = system.get_atom(17590).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "oxygen");
        assert_eq!(atom.get_element_symbol().unwrap(), "O");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 15.9994);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.152);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 2);

        // NA in NA5250
        let atom = system.get_atom(32795).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "sodium");
        assert_eq!(atom.get_element_symbol().unwrap(), "Na");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 22.9897);
        assert_eq!(atom.get_vdw(), None);
        assert_eq!(atom.get_expected_max_bonds(), None);

        // CL in CL5259
        let atom = system.get_atom(32804).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "carbon");
        assert_eq!(atom.get_element_symbol().unwrap(), "C");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 12.0107);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.17);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 4);

        // CL in CL5271
        let atom = system.get_atom(32816).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "carbon");
        assert_eq!(atom.get_element_symbol().unwrap(), "C");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 12.0107);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.17);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 4);
    }

    #[test]
    fn guess_elements_invalid_query() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();

        match system
            .guess_elements(Elements::from_file("test_files/elements_invalid_group.yaml").unwrap())
        {
            Ok(_) => panic!("Function should have failed."),
            Err(ElementError::InvalidQuery(SelectError::GroupNotFound(x))) => {
                assert_eq!(x, "Membrane");
            }
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }

        // check that the system is not changed
        for atom in system.atoms_iter() {
            assert!(atom.get_element_name().is_none());
            assert!(atom.get_element_symbol().is_none());
            assert!(atom.get_mass().is_none());
            assert!(atom.get_vdw().is_none());
            assert!(atom.get_expected_max_bonds().is_none());
            assert!(atom.get_position().is_some());
        }
    }

    #[test]
    fn guess_elements_complicated_groups() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system
            .guess_elements(
                Elements::from_file("test_files/elements_complicated_group.yaml").unwrap(),
            )
            .unwrap();

        for (a, atom) in system.atoms_iter().enumerate() {
            if a < 61 {
                assert_eq!(atom.get_element_name().unwrap(), "protein element");
                assert_eq!(atom.get_element_symbol().unwrap(), "P");
                assert!(atom.get_mass().is_none());
                assert!(atom.get_vdw().is_none());
                assert!(atom.get_expected_max_bonds().is_none());
            } else {
                assert_eq!(atom.get_element_name().unwrap(), "other");
                assert_eq!(atom.get_element_symbol().unwrap(), "O");
                assert!(atom.get_mass().is_none());
                assert!(atom.get_vdw().is_none());
                assert!(atom.get_expected_max_bonds().is_none());
            }
        }
    }

    #[test]
    fn guess_properties_1() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();

        for atom in system.atoms_iter_mut() {
            atom.set_element_name("carbon");
        }

        system
            .guess_properties(
                Elements::from_file("test_files/elements_properties_complete.yaml").unwrap(),
            )
            .unwrap();

        for atom in system.atoms_iter() {
            assert_approx_eq!(f32, atom.get_mass().unwrap(), 16.0107);
            assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.21);
            assert_eq!(atom.get_expected_max_bonds().unwrap(), 3);
            assert_eq!(atom.get_expected_min_bonds().unwrap(), 3);
        }
    }

    #[test]
    fn guess_properties_2() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.guess_elements(Elements::default()).unwrap();
        system
            .guess_properties(
                Elements::from_file("test_files/elements_properties_complete.yaml").unwrap(),
            )
            .unwrap();

        for atom in system.atoms_iter() {
            assert!(atom.get_element_name().is_some());
            assert!(atom.get_element_symbol().is_some());
            assert!(atom.get_mass().is_some());
            assert!(atom.get_vdw().is_some());
            assert!(atom.get_expected_max_bonds().is_some());
            assert!(atom.get_expected_min_bonds().is_some());
        }

        // N in SER1
        let atom = system.get_atom(0).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "nitrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "N");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 17.0067);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.255);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 5);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 4);

        // H1 in SER1
        let atom = system.get_atom(1).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "hydrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "H");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 1.5079);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.15);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 2);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 2);

        // C in SER23
        let atom = system.get_atom(360).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "carbon");
        assert_eq!(atom.get_element_symbol().unwrap(), "C");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 16.0107);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.21);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 3);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 3);

        // O31 in POPC44
        let atom = system.get_atom(3081).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "oxygen");
        assert_eq!(atom.get_element_symbol().unwrap(), "O");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 19.9994);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.08);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 4);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 3);

        // P in POPC127
        let atom = system.get_atom(14184).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "phosphorus");
        assert_eq!(atom.get_element_symbol().unwrap(), "P");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 32.9738);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.32);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 6);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 5);

        // HW1 in SOL4827
        let atom = system.get_atom(31541).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "hydrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "H");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 1.5079);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.15);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 2);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 2);

        // OW in SOL177
        let atom = system.get_atom(17590).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "oxygen");
        assert_eq!(atom.get_element_symbol().unwrap(), "O");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 19.9994);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.08);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 4);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 3);

        // NA in NA5250
        let atom = system.get_atom(32795).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "sodium");
        assert_eq!(atom.get_element_symbol().unwrap(), "Na");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 25.9897);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.21);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 0);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 0);

        // CL in CL5271
        let atom = system.get_atom(32816).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "chlorine");
        assert_eq!(atom.get_element_symbol().unwrap(), "Cl");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 37.453);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.20);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 0);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 0);
    }

    #[test]
    fn guess_properties_with_warnings() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.guess_elements(Elements::default()).unwrap();

        system.get_atom_mut(1).unwrap().reset_element_name();

        let no_element: Vec<usize> = vec![2];
        let not_recognized: Vec<usize> = vec![
            32789, 32790, 32791, 32792, 32793, 32794, 32795, 32796, 32797, 32798, 32799, 32800,
            32801, 32802,
        ];
        let no_mass: Vec<usize> = vec![
            32803, 32804, 32805, 32806, 32807, 32808, 32809, 32810, 32811, 32812, 32813, 32814,
            32815, 32816, 32817,
        ];
        let no_vdw: Vec<usize> = vec![
            383, 517, 651, 785, 919, 1053, 1187, 1321, 1455, 1589, 1723, 1857, 1991, 2125, 2259,
            2393, 2527, 2661, 2795, 2929, 3063, 3197, 3331, 3465, 3599, 3733, 3867, 4001, 4135,
            4269, 4403, 4537, 4671, 4805, 4939, 5073, 5207, 5341, 5475, 5609, 5743, 5877, 6011,
            6145, 6279, 6413, 6547, 6681, 6815, 6949, 7083, 7217, 7351, 7485, 7619, 7753, 7887,
            8021, 8155, 8289, 8423, 8557, 8691, 8825, 8959, 9093, 9227, 9361, 9495, 9629, 9763,
            9897, 10031, 10165, 10299, 10433, 10567, 10701, 10835, 10969, 11103, 11237, 11371,
            11505, 11639, 11773, 11907, 12041, 12175, 12309, 12443, 12577, 12711, 12845, 12979,
            13113, 13247, 13381, 13515, 13649, 13783, 13917, 14051, 14185, 14319, 14453, 14587,
            14721, 14855, 14989, 15123, 15257, 15391, 15525, 15659, 15793, 15927, 16061, 16195,
            16329, 16463, 16597, 16731, 16865, 16999, 17133, 17267, 17401, 32803, 32804, 32805,
            32806, 32807, 32808, 32809, 32810, 32811, 32812, 32813, 32814, 32815, 32816, 32817,
        ];
        let no_max_bonds: Vec<usize> = vec![
            32803, 32804, 32805, 32806, 32807, 32808, 32809, 32810, 32811, 32812, 32813, 32814,
            32815, 32816, 32817,
        ];
        let no_min_bonds: Vec<usize> = vec![
            383, 517, 651, 785, 919, 1053, 1187, 1321, 1455, 1589, 1723, 1857, 1991, 2125, 2259,
            2393, 2527, 2661, 2795, 2929, 3063, 3197, 3331, 3465, 3599, 3733, 3867, 4001, 4135,
            4269, 4403, 4537, 4671, 4805, 4939, 5073, 5207, 5341, 5475, 5609, 5743, 5877, 6011,
            6145, 6279, 6413, 6547, 6681, 6815, 6949, 7083, 7217, 7351, 7485, 7619, 7753, 7887,
            8021, 8155, 8289, 8423, 8557, 8691, 8825, 8959, 9093, 9227, 9361, 9495, 9629, 9763,
            9897, 10031, 10165, 10299, 10433, 10567, 10701, 10835, 10969, 11103, 11237, 11371,
            11505, 11639, 11773, 11907, 12041, 12175, 12309, 12443, 12577, 12711, 12845, 12979,
            13113, 13247, 13381, 13515, 13649, 13783, 13917, 14051, 14185, 14319, 14453, 14587,
            14721, 14855, 14989, 15123, 15257, 15391, 15525, 15659, 15793, 15927, 16061, 16195,
            16329, 16463, 16597, 16731, 16865, 16999, 17133, 17267, 17401, 32803, 32804, 32805,
            32806, 32807, 32808, 32809, 32810, 32811, 32812, 32813, 32814, 32815, 32816, 32817,
        ];

        match system.guess_properties(
            Elements::from_file("test_files/elements_properties_incomplete.yaml").unwrap(),
        ) {
            Ok(_) => panic!("Function should have failed."),
            Err(ElementError::PropertiesGuessWarning(x)) => {
                assert_eq!(x.no_element, no_element);
                assert_eq!(x.not_recognized, not_recognized);
                assert_eq!(x.no_mass, no_mass);
                assert_eq!(x.no_vdw, no_vdw);
                assert_eq!(x.no_max_bonds, no_max_bonds);
                assert_eq!(x.no_min_bonds, no_min_bonds);
            }
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }

        // N in SER1
        let atom = system.get_atom(0).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "nitrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "N");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 17.0067);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.255);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 5);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 5);

        // H1 in SER1
        let atom = system.get_atom(1).unwrap();
        assert_eq!(atom.get_element_name(), None);
        assert_eq!(atom.get_element_symbol().unwrap(), "H");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 1.0079);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.1);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 1);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // C in SER23
        let atom = system.get_atom(360).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "carbon");
        assert_eq!(atom.get_element_symbol().unwrap(), "C");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 16.0107);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.21);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 3);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 2);

        // O31 in POPC44
        let atom = system.get_atom(3081).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "oxygen");
        assert_eq!(atom.get_element_symbol().unwrap(), "O");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 19.9994);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.08);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 4);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 2);

        // P in POPC127
        let atom = system.get_atom(14184).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "phosphorus");
        assert_eq!(atom.get_element_symbol().unwrap(), "P");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 32.9738);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.1871);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 6);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 2);

        // HW1 in SOL4827
        let atom = system.get_atom(31541).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "hydrogen");
        assert_eq!(atom.get_element_symbol().unwrap(), "H");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 1.5079);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.15);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 2);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 1);

        // OW in SOL177
        let atom = system.get_atom(17590).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "oxygen");
        assert_eq!(atom.get_element_symbol().unwrap(), "O");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 19.9994);
        assert_approx_eq!(f32, atom.get_vdw().unwrap(), 0.08);
        assert_eq!(atom.get_expected_max_bonds().unwrap(), 4);
        assert_eq!(atom.get_expected_min_bonds().unwrap(), 2);

        // NA in NA5250
        let atom = system.get_atom(32795).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "sodium");
        assert_eq!(atom.get_element_symbol().unwrap(), "Na");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 22.9897);
        assert_eq!(atom.get_vdw(), None);
        assert_eq!(atom.get_expected_max_bonds(), None);
        assert_eq!(atom.get_expected_min_bonds(), None);

        // CL in CL5271
        let atom = system.get_atom(32816).unwrap();
        assert_eq!(atom.get_element_name().unwrap(), "chlorine");
        assert_eq!(atom.get_element_symbol().unwrap(), "Cl");
        assert_approx_eq!(f32, atom.get_mass().unwrap(), 35.453);
        assert_eq!(atom.get_vdw(), None);
        assert_eq!(atom.get_expected_max_bonds(), None);
        assert_eq!(atom.get_expected_min_bonds(), None);
    }

    #[test]
    fn guess_bonds() {
        let mut system = System::from_file("test_files/aa_peptide.pdb").unwrap();

        system.guess_elements(Elements::default()).unwrap();
        system.guess_bonds(None).unwrap();

        let mut system_from_pdb = System::from_file("test_files/aa_peptide.pdb").unwrap();
        system_from_pdb
            .add_bonds_from_pdb("test_files/aa_peptide.pdb")
            .unwrap();

        for (atom1, atom2) in system.atoms_iter().zip(system_from_pdb.atoms_iter()) {
            assert_eq!(atom1.get_bonded(), atom2.get_bonded());
        }
    }

    #[test]
    #[cfg(feature = "parallel")]
    fn guess_bonds_parallel() {
        let mut system = System::from_file("test_files/aa_peptide.pdb").unwrap();

        system.guess_elements(Elements::default()).unwrap();
        system.guess_bonds_parallel(None, 4).unwrap();

        let mut system_from_pdb = System::from_file("test_files/aa_peptide.pdb").unwrap();
        system_from_pdb
            .add_bonds_from_pdb("test_files/aa_peptide.pdb")
            .unwrap();

        for (atom1, atom2) in system.atoms_iter().zip(system_from_pdb.atoms_iter()) {
            assert_eq!(atom1.get_bonded(), atom2.get_bonded());
        }
    }

    #[test]
    fn guess_bonds_warnings() {
        let mut system = System::from_file("test_files/aa_peptide.pdb").unwrap();

        system.guess_elements(Elements::default()).unwrap();

        let mut elements_for_bonds = Elements::default();
        elements_for_bonds.update(
            Elements::from_file("test_files/elements_update_guess_bonds_warning.yaml").unwrap(),
        );

        match system.guess_properties(elements_for_bonds) {
            Ok(_) | Err(_) => (),
        }

        system.get_atom_mut(1).unwrap().reset_vdw();

        let no_vdw = vec![2];
        let too_few_bonds = vec![
            2, 12, 31, 50, 61, 72, 91, 110, 121, 132, 151, 170, 192, 211, 230, 241, 252, 271, 290,
            301, 312, 331, 350, 361,
        ];
        let too_many_bonds = vec![
            1, 14, 33, 52, 63, 74, 93, 112, 123, 134, 153, 172, 188, 194, 213, 232, 243, 254, 273,
            292, 303, 314, 333, 352,
        ];

        match system.guess_bonds(None) {
            Ok(_) => panic!("Function should have returned a warning."),
            Err(ElementError::BondsGuessWarning(e)) => {
                assert_eq!(e.no_vdw, no_vdw);
                assert_eq!(
                    e.too_few_bonds.keys().cloned().collect::<Vec<usize>>(),
                    too_few_bonds
                );
                assert_eq!(
                    e.too_many_bonds.keys().cloned().collect::<Vec<usize>>(),
                    too_many_bonds
                );
            }
            Err(e) => panic!("Incorrect warning type `{:?}` returned.", e),
        }

        let mut system_from_pdb = System::from_file("test_files/aa_peptide.pdb").unwrap();
        system_from_pdb
            .add_bonds_from_pdb("test_files/aa_peptide.pdb")
            .unwrap();

        for (atom1, atom2) in system.atoms_iter().zip(system_from_pdb.atoms_iter()) {
            if atom1.get_atom_number() <= 2 {
                continue;
            }

            assert_eq!(atom1.get_bonded(), atom2.get_bonded());
        }
    }

    #[test]
    #[cfg(feature = "parallel")]
    fn guess_bonds_multithreaded() {
        let no_vdw = vec![2];
        let too_few_bonds = vec![
            2, 12, 31, 50, 61, 72, 91, 110, 121, 132, 151, 170, 192, 211, 230, 241, 252, 271, 290,
            301, 312, 331, 350, 361,
        ];
        let too_many_bonds = vec![
            1, 14, 33, 52, 63, 74, 93, 112, 123, 134, 153, 172, 188, 194, 213, 232, 243, 254, 273,
            292, 303, 314, 333, 352,
        ];

        for n_threads in [1, 2, 3, 4, 5, 6, 7, 8, 16, 32, 512, 1024] {
            let mut system = System::from_file("test_files/aa_peptide.pdb").unwrap();

            system.guess_elements(Elements::default()).unwrap();

            let mut elements_for_bonds = Elements::default();
            elements_for_bonds.update(
                Elements::from_file("test_files/elements_update_guess_bonds_warning.yaml").unwrap(),
            );

            match system.guess_properties(elements_for_bonds) {
                Ok(_) | Err(_) => (),
            }

            system.get_atom_mut(1).unwrap().reset_vdw();

            match system.guess_bonds_parallel(None, n_threads) {
                Ok(_) => panic!("Function should have returned a warning."),
                Err(ElementError::BondsGuessWarning(e)) => {
                    assert_eq!(e.no_vdw, no_vdw);
                    assert_eq!(
                        e.too_few_bonds.keys().cloned().collect::<Vec<usize>>(),
                        too_few_bonds
                    );
                    assert_eq!(
                        e.too_many_bonds.keys().cloned().collect::<Vec<usize>>(),
                        too_many_bonds
                    );
                }
                Err(e) => panic!("Incorrect warning type `{:?}` returned.", e),
            }

            let mut system_from_pdb = System::from_file("test_files/aa_peptide.pdb").unwrap();
            system_from_pdb
                .add_bonds_from_pdb("test_files/aa_peptide.pdb")
                .unwrap();

            for (atom1, atom2) in system.atoms_iter().zip(system_from_pdb.atoms_iter()) {
                if atom1.get_atom_number() <= 2 {
                    continue;
                }

                assert_eq!(atom1.get_bonded(), atom2.get_bonded());
            }
        }
    }

    #[test]
    fn guess_bonds_empty_system() {
        let mut system = System::new(
            "Empty system",
            vec![],
            Some(SimBox::from([10.0, 10.0, 10.0])),
        );
        system.guess_bonds(None).unwrap();
        assert!(!system.has_bonds());
    }

    #[test]
    #[cfg(feature = "parallel")]
    fn guess_bonds_empty_system_parallel() {
        let mut system = System::new(
            "Empty system",
            vec![],
            Some(SimBox::from([10.0, 10.0, 10.0])),
        );
        system.guess_bonds_parallel(None, 7).unwrap();
        assert!(!system.has_bonds());
    }

    #[test]
    #[should_panic(
        expected = "FATAL GROAN ERROR | System::guess_bonds_parallel | Number of threads should not be zero."
    )]
    #[cfg(feature = "parallel")]
    fn guess_bonds_no_threads() {
        let mut system = System::from_file("test_files/aa_peptide.pdb").unwrap();

        system.guess_bonds_parallel(None, 0).unwrap();
    }
}
