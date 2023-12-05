// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of the Element structure and its methods.

use serde::Deserialize;
use std::fs::File;
use std::io::Read;
use std::{collections::HashMap, path::Path};

use crate::{
    errors::ParseElementError,
    selections::select::{parse_query, Select},
};

/*/// ## Methods for working with elements in the System.
impl System {
    /// Guess elements of the atoms in the system based on the provided `SupportedElements` file.
    pub fn guess_elements(&mut self, elements: SupportedElements) -> Result<(), ElementError> {
        // only elements that actually exist in the system 
        // will be included in the `SupportedElements` structure associated with the system.
        let mut new_elements = Vec::new();

        for atom in self.atoms_iter_mut() {

        }

        
    }
}*/

/// Contains information about all elements that can occur in the system.
#[derive(Debug, Clone)]
#[allow(unused)]
pub struct SupportedElements {
    /// All supported elements
    elements: Vec<Element>,
}

/// Contains information about specific element.
#[derive(Debug, Clone)]
#[allow(unused)]
pub struct Element {
    /// Name of the element
    name: String,
    /// Symbol of the element
    symbol: Option<String>,
    /// `Select` structure identifying the atoms of this element
    select: Option<Select>,
    /// Atomic mass of the element in amu (daltons)
    mass: Option<f32>,
    /// Van der Waals radius of the atom in nm
    vdw: Option<f32>,
    /// Expected maximal number of bonds
    expected_max_bonds: Option<u8>,
}

impl Element {
    /// Get name of the element.
    pub fn get_name(&self) -> &str {
        &self.name
    }

    /// Get symbol of the element.
    pub fn get_symbol(&self) -> Option<&str> {
        self.symbol.as_deref()
    }

    /// Get select tree of the element.
    pub fn get_select(&self) -> Option<&Select> {
        self.select.as_ref()
    }

    /// Get mass of the element.
    pub fn get_mass(&self) -> Option<f32> {
        self.mass
    }

    /// Get van der Waals radius of the element.
    pub fn get_vdw(&self) -> Option<f32> {
        self.vdw
    }

    /// Get the expected maximal number of bonds.
    pub fn get_expected_max_bonds(&self) -> Option<u8> {
        self.expected_max_bonds
    }

}

impl Default for SupportedElements {
    /// Construct a default `SupportedElements` structure.
    /// The structure will contain default information about elements
    /// that are supported and recognized by the `groan_rs` library.
    ///
    /// ## Notes
    /// - This function parses YAML content from `src/config/elements.yaml`
    /// which is included in the `groan_rs` library at compile time.
    /// - This is a relatively slow operation and there is no reason to call it multiple times in a program!
    /// If you need to use the `SupportedElements` structure for multiple systems, clone it.
    fn default() -> Self {
        let yaml = include_str!("../config/elements.yaml");

        SupportedElements::new_from_string(yaml)
            .expect("FATAL GROAN ERROR | SupportedElements::default | Default `elements.yaml` file could not be read or parsed.")
    }
}

impl SupportedElements {
    /// Construct a new `SupportedElements` structure from the provided YAML file.
    ///
    /// ## Returns
    /// `SupportedElements` structure if parsing was successful.
    /// `ParseElementError` otherwise.
    ///
    /// ## Notes
    /// - For an example of the 'elements yaml file', see `src/config/elements.yaml`.
    pub fn new_from_file(filename: impl AsRef<Path>) -> Result<Self, ParseElementError> {
        SupportedElements::new_from_string(&SupportedElements::load_yaml_to_string(filename)?)
    }

    /// Parse yaml string into `SupportedElements` structure.
    fn new_from_string(yaml: &str) -> Result<Self, ParseElementError> {
        let elements = parse_yaml_raw(yaml)?
            .into_iter()
            .map(|(name, temp)| temp.to_element(name))
            .collect::<Result<Vec<Element>, ParseElementError>>()?;

        Ok(SupportedElements { elements })
    }

    /// Update `SupportedElements` structure using data from the provided YAML file.
    ///
    /// ## Returns
    /// `Ok` if the parsing was successful.
    /// `ParseElementError` otherwise. If an error occurs, the `SupportedElements` structure is not changed.
    ///
    /// ## Example
    /// Let's suppose that in the default `SupportedElements` structure provided by `groan_rs` library,
    /// you are missing information about polonium
    /// and also you don't like the default van der Waals radius for carbon.
    /// You can construct a new yaml file containing the missing and modified information.
    /// The yaml file, `my_elements.yaml`, should look like this:
    /// ```yaml
    /// ---
    /// carbon:
    ///   vdw: 0.20
    /// polonium:
    ///   symbol: Po
    ///   query: name r'^[Pp][Oo].*$'
    ///   mass: 209.0
    /// ...
    /// ```
    ///
    /// The modified `SupportedElements` structure can be obtained using:
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// // load default parameters for the elements
    /// let mut elements = SupportedElements::default();
    /// // update the `elements` structure with your custom information
    /// elements.update_from_file("my_elements.yaml").unwrap();
    /// ```
    ///
    /// In the `elements` structure, carbon will now have `vdw` of 0.2,
    /// and a new element, polonium, will be added.
    /// No other information about the elements will be changed.
    pub fn update_from_file(
        &mut self,
        filename: impl AsRef<Path>,
    ) -> Result<(), ParseElementError> {
        SupportedElements::update_from_string(
            self,
            &SupportedElements::load_yaml_to_string(filename)?,
        )
    }

    /// Opens the specified file and loads its contents into a string.
    fn load_yaml_to_string(filename: impl AsRef<Path>) -> Result<String, ParseElementError> {
        let mut yaml_file = match File::open(&filename) {
            Ok(x) => x,
            Err(_) => {
                return Err(ParseElementError::FileNotFound(Box::from(
                    filename.as_ref(),
                )))
            }
        };

        let mut yaml_string = String::new();
        match yaml_file.read_to_string(&mut yaml_string) {
            Ok(_) => Ok(yaml_string),
            Err(_) => {
                return Err(ParseElementError::FileCouldNotBeRead(Box::from(
                    filename.as_ref(),
                )))
            }
        }
    }

    /// Update `SupportedElements` structure from yaml string.
    /// Previously defined elements are updated. New elements are added.
    /// In case an error occurs, `self` is not changed.
    fn update_from_string(&mut self, yaml: &str) -> Result<(), ParseElementError> {
        let parsed_elements = parse_yaml_raw(yaml)?;

        // first check for the validity of all elements
        // this is necessary because we do not want to disrupt the SupportedElements
        // structure if the parsed yaml is invalid
        for (name, temp) in &parsed_elements {
            temp.validate(name)?;
        }

        // then update the structure
        for (name, temp) in parsed_elements {
            if let Some(element) = self.elements.iter_mut().find(|e| e.name == name) {
                temp.update_element(element);
            } else {
                let new_element = temp.to_element(name)
                    .expect("FATAL GROAN ERROR | SupportedElements::update_from_string | TempElement should be valid.");
                self.elements.push(new_element);
            }
        }

        Ok(())
    }
}

/// Temporary `Element` structure for parsing YAML.
#[derive(Deserialize, Debug)]
#[serde(deny_unknown_fields)]
struct TempElement {
    symbol: Option<String>,
    query: Option<String>,
    mass: Option<f32>,
    vdw: Option<f32>,
    expected_max_bonds: Option<u8>,
}

impl TempElement {
    /// Convert `TempElement` and name to `Element` structure.
    /// Checks for validity of the query, mass, and vdw radius.
    #[allow(clippy::wrong_self_convention)]
    fn to_element(self, name: String) -> Result<Element, ParseElementError> {
        let select = TempElement::query2select(self.query.as_deref())?;

        TempElement::check_mass(self.mass, &name)?;
        TempElement::check_vdw(self.vdw, &name)?;

        Ok(Element {
            name,
            symbol: self.symbol,
            select,
            mass: self.mass,
            vdw: self.vdw,
            expected_max_bonds: self.expected_max_bonds,
        })
    }

    /// Update `Element` structure based on the provided `TempElement` structure.
    ///
    /// ## Panics
    /// Panics if the `TempElement` is invalid (its validity should be checked before using `TempElement::validate`).
    fn update_element(self, element: &mut Element) {
        if self.symbol.is_some() {
            element.symbol = self.symbol;
        }

        if self.query.is_some() {
            element.select = TempElement::query2select(self.query.as_deref())
                .expect("FATAL GROAN ERROR | TempElement::update_element | Query should be valid.");
        }

        if self.mass.is_some() {
            TempElement::check_mass(self.mass, &element.name)
                .expect("FATAL GROAN ERROR | TempElement::update_element | Mass should be valid.");
            element.mass = self.mass;
        }

        if self.vdw.is_some() {
            TempElement::check_vdw(self.vdw, &element.name)
                .expect("FATAL GROAN ERROR | TempElement::update_element | VdW should be valid.");
            element.vdw = self.vdw;
        }

        if self.expected_max_bonds.is_some() {
            element.expected_max_bonds = self.expected_max_bonds;
        }
    }

    /// Checks that the properties of the TempElement are valid.
    fn validate(&self, name: &str) -> Result<(), ParseElementError> {
        TempElement::query2select(self.query.as_deref())?;
        TempElement::check_mass(self.mass, name)?;
        TempElement::check_vdw(self.vdw, name)?;

        Ok(())
    }

    /// Convert query to select tree. If `query` is None, returns `Ok(None)`.
    fn query2select(query: Option<&str>) -> Result<Option<Select>, ParseElementError> {
        match query {
            Some(x) => match parse_query(x) {
                Ok(x) => Ok(Some(*x)),
                Err(e) => Err(ParseElementError::InvalidQuery(e)),
            },
            None => Ok(None),
        }
    }

    /// Check that `mass` is not negative. If `mass` is None, returns `Ok`.
    fn check_mass(mass: Option<f32>, name: &str) -> Result<(), ParseElementError> {
        match mass {
            Some(x) if x < 0.0 => Err(ParseElementError::InvalidMass(name.to_owned())),
            Some(_) | None => Ok(()),
        }
    }

    /// Check that `vdw` is not negative. If `mass` is None, returns `Ok`.
    fn check_vdw(vdw: Option<f32>, name: &str) -> Result<(), ParseElementError> {
        match vdw {
            Some(x) if x < 0.0 => Err(ParseElementError::InvalidVdW(name.to_owned())),
            Some(_) | None => Ok(()),
        }
    }
}

fn parse_yaml_raw(yaml: &str) -> Result<HashMap<String, TempElement>, ParseElementError> {
    // if parsing yaml fails, return location of the error
    let raw_elements: HashMap<String, TempElement> = match serde_yaml::from_str(yaml) {
        Ok(x) => x,
        Err(e) => {
            let (line, col) = match e.location() {
                Some(x) => (x.line(), x.column()),
                None => (0, 0),
            };

            return Err(ParseElementError::CouldNotParseYaml(line, col));
        }
    };

    Ok(raw_elements)
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use crate::errors::SelectError;
    use float_cmp::assert_approx_eq;

    use super::*;

    #[test]
    fn elements_default() {
        let elements = SupportedElements::default();

        assert_eq!(elements.elements.len(), 39);

        let names: Vec<String> = elements.elements.iter().map(|x| x.name.clone()).collect();

        assert!(names.contains(&String::from("dummy")));
        assert!(names.contains(&String::from("hydrogen")));
        assert!(names.contains(&String::from("carbon")));
        assert!(names.contains(&String::from("bismuth")));
        assert!(names.contains(&String::from("lithium")));

        let mut visited = 0;

        for e in &elements.elements {
            // first element
            if e.name == "dummy" {
                assert!(e.symbol.is_none());
                assert!(e.select.is_some());
                assert_approx_eq!(f32, e.mass.unwrap(), 0.0);
                assert!(e.vdw.is_none());
                assert!(e.expected_max_bonds.is_none());
                visited += 1;

                continue;
            }

            // important element
            if e.name == "carbon" {
                assert_eq!(e.symbol.as_ref().unwrap(), "C");
                assert!(e.select.is_some());
                assert_approx_eq!(f32, e.mass.unwrap(), 12.0107);
                assert_approx_eq!(f32, e.vdw.unwrap(), 0.17);
                assert_eq!(e.expected_max_bonds.unwrap(), 4);
                visited += 1;

                continue;
            }

            // last element
            if e.name == "bismuth" {
                assert_eq!(e.symbol.as_ref().unwrap(), "Bi");
                assert!(e.select.is_some());
                assert_approx_eq!(f32, e.mass.unwrap(), 208.98040);
                assert!(e.vdw.is_none());
                assert!(e.expected_max_bonds.is_none());
                visited += 1;

                continue;
            }
        }

        assert_eq!(visited, 3);
    }

    fn elements_match(el1: &mut SupportedElements, el2: &mut SupportedElements) {
        assert_eq!(
            el1.elements.len(),
            el2.elements.len()
        );

        el1
            .elements
            .sort_by(|x, y| x.name.cmp(&y.name));
        el2
            .elements
            .sort_by(|x, y| x.name.cmp(&y.name));

        for (e1, e2) in el1
            .elements
            .iter()
            .zip(el2.elements.iter())
        {
            assert_eq!(e1.name, e2.name);

            match (&e1.symbol, &e2.symbol) {
                (Some(x), Some(y)) => assert_eq!(x, y),
                (None, None) => (),
                _ => panic!("Symbols do not match."),
            }

            match (&e1.select, &e2.select) {
                (Some(x), Some(y)) => assert_eq!(x, y),
                (None, None) => (),
                _ => panic!("Select trees do not match."),
            }

            match (&e1.mass, &e2.mass) {
                (Some(x), Some(y)) => assert_approx_eq!(f32, *x, *y),
                (None, None) => (),
                _ => panic!("Masses do not match."),
            }

            match (&e1.vdw, &e2.vdw) {
                (Some(x), Some(y)) => assert_approx_eq!(f32, *x, *y),
                (None, None) => (),
                _ => panic!("Van der Waals radii do not match."),
            }

            match (&e1.expected_max_bonds, &e2.expected_max_bonds) {
                (Some(x), Some(y)) => assert_eq!(x, y),
                (None, None) => (),
                _ => panic!("Expected max bonds do not match."),
            }
        }
    }

    #[test]
    fn elements_from_file() {
        let mut default_elements = SupportedElements::default();
        let mut elements_from_file =
            SupportedElements::new_from_file("src/config/elements.yaml").unwrap();

        elements_match(&mut default_elements, &mut elements_from_file);
    }

    #[test]
    fn elements_invalid_field() {
        match SupportedElements::new_from_file("test_files/elements_invalid_field.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::CouldNotParseYaml(line, col)) => {
                assert_eq!(line, 10);
                assert_eq!(col, 3);
            }
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn elements_nonexistent_file() {
        match SupportedElements::new_from_file("nonexistent_file.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::FileNotFound(name)) => {
                assert_eq!(name.to_str().unwrap(), "nonexistent_file.yaml")
            }
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn elements_invalid_query() {
        match SupportedElements::new_from_file("test_files/elements_invalid_query.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::InvalidQuery(SelectError::InvalidNumber(query))) => {
                assert_eq!(query, "resid 14 15 17 X 18 19 20")
            }
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn elements_invalid_mass() {
        match SupportedElements::new_from_file("test_files/elements_invalid_mass.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::InvalidMass(name)) => assert_eq!(name, String::from("hydrogen")),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn elements_invalid_vdw() {
        match SupportedElements::new_from_file("test_files/elements_invalid_vdw.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::InvalidVdW(name)) => assert_eq!(name, String::from("carbon")),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn elements_update() {
        let mut elements = SupportedElements::default();
        elements
            .update_from_file("test_files/elements_update.yaml")
            .unwrap();

        let mut visited = 0;

        for e in &elements.elements {
            // first element
            if e.name == "dummy" {
                assert!(e.symbol.is_none());
                assert!(e.select.is_some());
                assert_approx_eq!(f32, e.mass.unwrap(), 0.0);
                assert!(e.vdw.is_none());
                assert!(e.expected_max_bonds.is_none());

                visited += 1;
                continue;
            }

            // changed element #1
            if e.name == "carbon" {
                assert_eq!(e.symbol.as_ref().unwrap(), "C");
                assert!(e.select.is_some());
                assert_approx_eq!(f32, e.mass.unwrap(), 12.0107);
                assert_approx_eq!(f32, e.vdw.unwrap(), 0.20);
                assert_eq!(e.expected_max_bonds.unwrap(), 4);

                visited += 1;
                continue;
            }

            // change element #2
            if e.name == "hydrogen" {
                assert_eq!(e.symbol.as_ref().unwrap(), "HH");
                assert!(e.select.is_some());
                assert_approx_eq!(f32, e.mass.unwrap(), 2.014);
                assert_approx_eq!(f32, e.vdw.unwrap(), 0.12);
                assert_eq!(e.expected_max_bonds.unwrap(), 8);

                visited += 1;
                continue;
            }

            // added element
            if e.name == "polonium" {
                assert_eq!(e.symbol.as_ref().unwrap(), "Po");
                assert!(e.select.is_some());
                assert_approx_eq!(f32, e.mass.unwrap(), 209.0);
                assert!(e.vdw.is_none());
                assert!(e.expected_max_bonds.is_none());

                visited += 1;
                continue;
            }

            // last element
            if e.name == "bismuth" {
                assert_eq!(e.symbol.as_ref().unwrap(), "Bi");
                assert!(e.select.is_some());
                assert_approx_eq!(f32, e.mass.unwrap(), 208.98040);
                assert!(e.vdw.is_none());
                assert!(e.expected_max_bonds.is_none());

                visited += 1;
                continue;
            }
        }

        assert_eq!(visited, 5);
    }

    #[test]
    fn elements_update_invalid_field() {
        let mut elements = SupportedElements::default();

        match elements.update_from_file("test_files/elements_invalid_field.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::CouldNotParseYaml(line, col)) => {
                assert_eq!(line, 10);
                assert_eq!(col, 3);
            }
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }

        let mut default = SupportedElements::default();
        elements_match(&mut elements, &mut default);
    }

    #[test]
    fn elements_update_nonexistent_file() {
        let mut elements = SupportedElements::default();

        match elements.update_from_file("nonexistent_file.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::FileNotFound(name)) => {
                assert_eq!(name.to_str().unwrap(), "nonexistent_file.yaml")
            }
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }

        let mut default = SupportedElements::default();
        elements_match(&mut elements, &mut default);
    }

    #[test]
    fn elements_update_invalid_query() {
        let mut elements = SupportedElements::default();

        match elements.update_from_file("test_files/elements_invalid_query.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::InvalidQuery(SelectError::InvalidNumber(query))) => {
                assert_eq!(query, "resid 14 15 17 X 18 19 20")
            }
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }

        let mut default = SupportedElements::default();
        elements_match(&mut elements, &mut default);
    }

    #[test]
    fn elements_update_invalid_mass() {
        let mut elements = SupportedElements::default();

        match SupportedElements::new_from_file("test_files/elements_invalid_mass.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::InvalidMass(name)) => assert_eq!(name, String::from("hydrogen")),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }

        let mut default = SupportedElements::default();
        elements_match(&mut elements, &mut default);
    }

    #[test]
    fn elements_update_invalid_vdw() {
        let mut elements = SupportedElements::default();

        match SupportedElements::new_from_file("test_files/elements_invalid_vdw.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::InvalidVdW(name)) => assert_eq!(name, String::from("carbon")),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }

        let mut default = SupportedElements::default();
        elements_match(&mut elements, &mut default);
    }

}
