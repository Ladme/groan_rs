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

/// Contains information about all elements that can occur in the system.
#[derive(Debug)]
#[allow(unused)]
pub struct ElementsProperties {
    /// all supported elements
    elements: Vec<Element>,
}

/// Contains information about specific element.
#[derive(Debug, Clone)]
#[allow(unused)]
pub struct Element {
    /// name of the element
    name: String,
    /// symbol of the element
    symbol: Option<String>,
    /// `Select` structure identifying the atoms of this element
    select: Option<Select>,
    /// atomic mass of the element in amu (daltons)
    mass: Option<f32>,
    /// van der Waals radius of the atom in nm
    vdw: Option<f32>,
    /// expected maximal number of bonds
    expected_max_bonds: Option<u8>,
}

impl Default for ElementsProperties {
    /// Construct a default `ElementsProperties` structure.
    ///
    /// This function parses YAML content from `src/config/elements.yaml`
    /// which is included in the `groan_rs` library at compile time.
    ///
    /// This is a relatively slow operation and there is no reason to call it multiple times in a program!
    fn default() -> Self {
        let yaml = include_str!("../config/elements.yaml");

        ElementsProperties::new_from_string(yaml)
            .expect("FATAL GROAN ERROR | ElementsProperties::default | Default `elements.yaml` file could not be read or parsed.")
    }
}

impl ElementsProperties {
    /// Construct a new `ElementsProperties` structure from the provided YAML file.
    ///
    /// ## Returns
    /// `ElementsProperties` structure if parsing was successful.
    /// `ParseElementError` otherwise.
    ///
    /// ## Notes
    /// - For an example of the 'elements yaml file', see `src/config/elements.yaml`.
    pub fn new_from_file(filename: impl AsRef<Path>) -> Result<Self, ParseElementError> {
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
            Ok(_) => ElementsProperties::new_from_string(&yaml_string),
            Err(_) => {
                return Err(ParseElementError::FileCouldNotBeRead(Box::from(
                    filename.as_ref(),
                )))
            }
        }
    }

    /// Parse yaml string into `ElementsProperties` structure.
    fn new_from_string(yaml: &str) -> Result<Self, ParseElementError> {
        #[derive(Deserialize, Debug)]
        #[serde(deny_unknown_fields)]
        struct TempElement {
            symbol: Option<String>,
            query: Option<String>,
            mass: Option<f32>,
            vdw: Option<f32>,
            expected_max_bonds: Option<u8>,
        }

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

        // convert temporary elements into final elements
        // parse selection queries to select trees
        // check for validity of the mass and vdw
        let elements = raw_elements
            .into_iter()
            .map(|(name, element)| {
                let select = match element.query {
                    Some(x) => match parse_query(&x) {
                        Ok(x) => Some(*x),
                        Err(e) => return Err(ParseElementError::InvalidQuery(e)),
                    },
                    None => None,
                };

                match element.mass {
                    Some(x) if x < 0.0 => return Err(ParseElementError::InvalidMass(name)),
                    Some(_) | None => (),
                }

                match element.vdw {
                    Some(x) if x < 0.0 => return Err(ParseElementError::InvalidVdW(name)),
                    Some(_) | None => (),
                }

                Ok(Element {
                    name,
                    symbol: element.symbol,
                    select,
                    mass: element.mass,
                    vdw: element.vdw,
                    expected_max_bonds: element.expected_max_bonds,
                })
            })
            .collect::<Result<Vec<Element>, ParseElementError>>()?;

        Ok(ElementsProperties { elements })
    }
}

/******************************/
/*   UNIT TESTS FOR ELEMENT   */
/******************************/

#[cfg(test)]
mod tests {
    use crate::errors::SelectError;
    use float_cmp::assert_approx_eq;

    use super::*;

    #[test]
    fn elements_default() {
        let elements = ElementsProperties::default();

        assert_eq!(elements.elements.len(), 39);

        let names: Vec<String> = elements.elements.iter().map(|x| x.name.clone()).collect();

        assert!(names.contains(&String::from("dummy")));
        assert!(names.contains(&String::from("hydrogen")));
        assert!(names.contains(&String::from("carbon")));
        assert!(names.contains(&String::from("bismuth")));
        assert!(names.contains(&String::from("lithium")));

        for e in &elements.elements {
            // first element
            if e.name == "dummy" {
                assert!(e.symbol.is_none());
                assert!(e.select.is_some());
                assert_approx_eq!(f32, e.mass.unwrap(), 0.0);
                assert!(e.vdw.is_none());
                assert!(e.expected_max_bonds.is_none());
                continue;
            }

            // important element
            if e.name == "carbon" {
                assert_eq!(e.symbol.as_ref().unwrap(), "C");
                assert!(e.select.is_some());
                assert_approx_eq!(f32, e.mass.unwrap(), 12.0107);
                assert_approx_eq!(f32, e.vdw.unwrap(), 0.17);
                assert_eq!(e.expected_max_bonds.unwrap(), 4);
                continue;
            }

            // last element
            if e.name == "bismuth" {
                assert_eq!(e.symbol.as_ref().unwrap(), "Bi");
                assert!(e.select.is_some());
                assert_approx_eq!(f32, e.mass.unwrap(), 208.98040);
                assert!(e.vdw.is_none());
                assert!(e.expected_max_bonds.is_none());
                continue;
            }
        }
    }

    #[test]
    fn elements_from_file() {
        let mut default_elements = ElementsProperties::default();
        let mut elements_from_file =
            ElementsProperties::new_from_file("src/config/elements.yaml").unwrap();

        assert_eq!(
            default_elements.elements.len(),
            elements_from_file.elements.len()
        );

        default_elements
            .elements
            .sort_by(|x, y| x.name.cmp(&y.name));
        elements_from_file
            .elements
            .sort_by(|x, y| x.name.cmp(&y.name));

        for (e1, e2) in default_elements
            .elements
            .iter()
            .zip(elements_from_file.elements.iter())
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
    fn elements_invalid_field() {
        match ElementsProperties::new_from_file("test_files/elements_invalid_field.yaml") {
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
        match ElementsProperties::new_from_file("nonexistent_file.yaml") {
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
        match ElementsProperties::new_from_file("test_files/elements_invalid_query.yaml") {
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
        match ElementsProperties::new_from_file("test_files/elements_invalid_mass.yaml") {
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
        match ElementsProperties::new_from_file("test_files/elements_invalid_vdw.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::InvalidVdW(name)) => assert_eq!(name, String::from("carbon")),
            Err(e) => panic!(
                "Function failed successfully but incorrect error type `{}` was returned.",
                e
            ),
        }
    }
}
