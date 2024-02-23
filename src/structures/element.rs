// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of the Element structure and its methods.

use indexmap::IndexMap;
use serde::de::{MapAccess, Visitor};
use serde::{de, Deserialize, Deserializer};
use std::fmt;
use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::errors::ElementError;
use crate::{
    errors::ParseElementError,
    selections::select::Select,
    structures::group::Group,
    system::general::System,
};

/// Contains information about all elements that can occur in the system.
/// For atomistic systems, construct with `Elements::default`.
#[derive(Debug, Clone)]
pub struct Elements {
    /// All supported elements.
    /// Keys are names of the elements.
    pub(crate) elements: IndexMap<String, Element>,
}

impl Default for Elements {
    /// Construct a default `Elements` structure.
    /// The structure will contain default information about elements
    /// that are supported and recognized by the `groan_rs` library.
    ///
    /// ## Notes
    /// - This function parses YAML content from `src/config/elements.yaml`
    /// which is included in the `groan_rs` library at compile time.
    /// - This is a relatively slow operation and there is no reason to call it multiple times in a program!
    /// If you need to use the `Elements` structure for multiple systems, clone it.
    fn default() -> Self {
        let yaml = include_str!("../config/elements.yaml");

        Elements::new_from_string(yaml)
            .expect("FATAL GROAN ERROR | Elements::default | Default `elements.yaml` file could not be read or parsed.")
    }
}

impl Elements {
    /// Construct a new `Elements` structure from the provided YAML file.
    ///
    /// ## Returns
    /// `Elements` structure if parsing was successful.
    /// `ParseElementError` otherwise.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let elements = match Elements::from_file("my_custom_elements.yaml") {
    ///     Ok(x) => x,
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;
    ///     }
    /// };
    /// ```
    ///
    /// ## Notes
    /// - For an example of the 'elements yaml file', see `src/config/elements.yaml`.
    pub fn from_file(filename: impl AsRef<Path>) -> Result<Self, ParseElementError> {
        Elements::new_from_string(&Elements::load_yaml_to_string(filename)?)
    }

    /// Parse yaml string into `Elements` structure.
    fn new_from_string(yaml: &str) -> Result<Self, ParseElementError> {
        let elements: IndexMap<String, Element> = match serde_yaml::from_str(yaml) {
            Ok(x) => x,
            Err(e) => return Err(ParseElementError::CouldNotParseYaml(e)),
        };

        Ok(Elements { elements })
    }

    /// Update `Elements` structure using data from another `Elements` structure.
    ///
    /// ## Example
    /// Let's suppose that in the default `Elements` structure provided by `groan_rs` library,
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
    /// The modified `Elements` structure can be obtained using:
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// // load default parameters for the elements
    /// let mut elements = Elements::default();
    /// // update the `elements` structure with your custom information
    /// elements.update(Elements::from_file("my_elements.yaml").unwrap());
    /// ```
    ///
    /// In the `elements` structure, carbon will now have `vdw` of 0.2,
    /// and a new element, polonium, will be added.
    /// No other information about the elements will be changed.
    pub fn update(&mut self, update_elements: Elements) {
        for (name, element) in update_elements.elements.into_iter() {
            if let Some(old) = self.elements.get_mut(&name) {
                old.update(element);
            } else {
                self.elements.insert(name, element);
            }
        }
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
}

/// Contains information about specific element.
#[derive(Debug, Clone)]
pub(crate) struct Element {
    /// Symbol of the element
    pub symbol: Option<String>,
    /// `Select` structure identifying the atoms of this element
    pub select: Option<Select>,
    /// Atomic mass of the element in amu (daltons)
    pub mass: Option<f32>,
    /// Van der Waals radius of the atom in nm
    pub vdw: Option<f32>,
    /// Expected maximal number of bonds
    pub expected_max_bonds: Option<u8>,
    /// Expected minimal number of bonds
    pub expected_min_bonds: Option<u8>,
}

impl Element {
    /// Check whether the select tree associated with the element is valid.
    /// This actually only checks for the existence of the groups in the `Select` structure
    /// as the other sanity checks should be performed before even constructing the `Element`.
    pub(crate) fn validate_select(&self, system: &System) -> Result<(), ElementError> {
        if let Some(select) = &self.select {
            match select.validate_groups(system) {
                Ok(_) => return Ok(()),
                Err(e) => return Err(ElementError::InvalidQuery(e)),
            }
        }

        Ok(())
    }

    /// Check whether the select tree associated with the element matches properties of the atom.
    /// Returns `true` if it does, returns `false` if not or if the select tree is not available.
    /// Returns `ElementError` if a nonexistent group is present in the select tree.
    pub(crate) fn matches(&self, atom_index: usize, system: &System) -> Result<bool, ElementError> {
        match &self.select {
            Some(select) => match Group::matches_select(atom_index, select, system) {
                Ok(x) => Ok(x),
                Err(e) => Err(ElementError::InvalidQuery(e)),
            },
            None => Ok(false),
        }
    }

    /// Update fields of `self` based on another `Element` structure.
    fn update(&mut self, element: Element) {
        if element.symbol.is_some() {
            self.symbol = element.symbol;
        }

        if element.select.is_some() {
            self.select = element.select;
        }

        if element.mass.is_some() {
            self.mass = element.mass;
        }

        if element.vdw.is_some() {
            self.vdw = element.vdw;
        }

        if element.expected_max_bonds.is_some() {
            self.expected_max_bonds = element.expected_max_bonds;
        }

        if element.expected_min_bonds.is_some() {
            self.expected_min_bonds = element.expected_min_bonds;
        }
    }
}

/// Handles parsing of the yaml file and basic sanity checks.
struct ElementVisitor;

impl<'de> Visitor<'de> for ElementVisitor {
    type Value = Element;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("struct Element")
    }

    fn visit_map<V>(self, mut map: V) -> Result<Element, V::Error>
    where
        V: MapAccess<'de>,
    {
        let mut symbol = None;
        let mut select = None;
        let mut mass = None;
        let mut vdw = None;
        let mut expected_max_bonds = None;
        let mut expected_min_bonds = None;

        while let Some(key) = map.next_key::<String>()? {
            match key.as_str() {
                "symbol" => symbol = map.next_value()?,
                "query" => {
                    let query: Option<String> = map.next_value()?;
                    select = match query {
                        Some(q) => match Select::parse_query(&q) {
                            Ok(sel) => Some(*sel),
                            Err(_) => return Err(de::Error::custom("invalid query")),
                        },
                        None => None,
                    };
                }
                "mass" => {
                    let value: Option<f32> = map.next_value()?;
                    mass = match value {
                        Some(v) if v < 0.0 => return Err(de::Error::custom("mass is negative")),
                        _ => value,
                    };
                }
                "vdw" => {
                    let value: Option<f32> = map.next_value()?;
                    vdw = match value {
                        Some(v) if v < 0.0 => return Err(de::Error::custom("vdw is negative")),
                        _ => value,
                    };
                }
                "expected_max_bonds" => expected_max_bonds = map.next_value()?,
                "expected_min_bonds" => expected_min_bonds = map.next_value()?,
                _ => return Err(de::Error::unknown_field(&key, FIELDS)),
            }
        }

        Ok(Element {
            symbol,
            select,
            mass,
            vdw,
            expected_max_bonds,
            expected_min_bonds,
        })
    }
}

const FIELDS: &[&str] = &[
    "symbol",
    "query",
    "mass",
    "vdw",
    "expected_max_bonds",
    "expected_min_bonds",
];

impl<'de> Deserialize<'de> for Element {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_struct("Element", FIELDS, ElementVisitor)
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use float_cmp::assert_approx_eq;

    use super::*;

    #[test]
    fn elements_default() {
        let elements = Elements::default();

        assert_eq!(elements.elements.len(), 39);

        // first element
        let e = elements.elements.get("dummy").unwrap();
        assert!(e.symbol.is_none());
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 0.0);
        assert!(e.vdw.is_none());
        assert!(e.expected_max_bonds.is_none());
        assert!(e.expected_min_bonds.is_none());

        // important element
        let e = elements.elements.get("carbon").unwrap();
        assert_eq!(e.symbol.as_ref().unwrap(), "C");
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 12.0107);
        assert_approx_eq!(f32, e.vdw.unwrap(), 0.17);
        assert_eq!(e.expected_max_bonds.unwrap(), 4);
        assert_eq!(e.expected_min_bonds.unwrap(), 2);

        // last element
        let e = elements.elements.get("bismuth").unwrap();
        assert_eq!(e.symbol.as_ref().unwrap(), "Bi");
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 208.98040);
        assert!(e.vdw.is_none());
        assert!(e.expected_max_bonds.is_none());
        assert!(e.expected_min_bonds.is_none());
    }

    fn elements_match(el1: &mut Elements, el2: &mut Elements) {
        assert_eq!(el1.elements.len(), el2.elements.len());

        for (name1, e1) in el1.elements.iter() {
            let e2 = el2.elements.get(name1).unwrap();

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
        let mut default_elements = Elements::default();
        let mut elements_from_file = Elements::from_file("src/config/elements.yaml").unwrap();

        elements_match(&mut default_elements, &mut elements_from_file);
    }

    #[test]
    fn elements_invalid_field() {
        match Elements::from_file("test_files/elements_invalid_field.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::CouldNotParseYaml(e)) => {
                let string = e.to_string();
                assert!(string.contains("unknown field"));
            }
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn elements_nonexistent_file() {
        match Elements::from_file("nonexistent_file.yaml") {
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
        match Elements::from_file("test_files/elements_invalid_query.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::CouldNotParseYaml(e)) => {
                let string = e.to_string();
                assert!(string.contains("invalid query"));
            }
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn elements_invalid_mass() {
        match Elements::from_file("test_files/elements_invalid_mass.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::CouldNotParseYaml(e)) => {
                let string = e.to_string();
                assert!(string.contains("mass is negative"));
            }
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn elements_invalid_vdw() {
        match Elements::from_file("test_files/elements_invalid_vdw.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::CouldNotParseYaml(e)) => {
                let string = e.to_string();
                assert!(string.contains("vdw is negative"));
            }
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn elements_update() {
        let mut elements = Elements::default();
        elements.update(Elements::from_file("test_files/elements_update.yaml").unwrap());

        assert_eq!(elements.elements.len(), 40);

        // first element
        let e = elements.elements.get("dummy").unwrap();
        assert!(e.symbol.is_none());
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 0.0);
        assert!(e.vdw.is_none());
        assert!(e.expected_max_bonds.is_none());
        assert!(e.expected_min_bonds.is_none());

        // changed element #1
        let e = elements.elements.get("carbon").unwrap();
        assert_eq!(e.symbol.as_ref().unwrap(), "C");
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 12.0107);
        assert_approx_eq!(f32, e.vdw.unwrap(), 0.20);
        assert_eq!(e.expected_max_bonds.unwrap(), 4);
        assert_eq!(e.expected_min_bonds.unwrap(), 2);

        // changed element #2
        let e = elements.elements.get("hydrogen").unwrap();
        assert_eq!(e.symbol.as_ref().unwrap(), "HH");
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 2.014);
        assert_approx_eq!(f32, e.vdw.unwrap(), 0.10);
        assert_eq!(e.expected_max_bonds.unwrap(), 8);
        assert_eq!(e.expected_min_bonds.unwrap(), 2);

        // added element
        let e = elements.elements.get("polonium").unwrap();
        assert_eq!(e.symbol.as_ref().unwrap(), "Po");
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 209.0);
        assert!(e.vdw.is_none());
        assert!(e.expected_max_bonds.is_none());
        assert!(e.expected_min_bonds.is_none());

        // last element
        let e = elements.elements.get("bismuth").unwrap();
        assert_eq!(e.symbol.as_ref().unwrap(), "Bi");
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 208.98040);
        assert!(e.vdw.is_none());
        assert!(e.expected_max_bonds.is_none());
        assert!(e.expected_max_bonds.is_none());
    }
}
