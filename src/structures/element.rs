// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of the Element structure and its methods.

use serde::{Deserialize, Deserializer, de};
use serde::de::{Visitor, MapAccess};
use std::fmt;
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
    /// All supported elements.
    /// Keys are names of the elements.
    elements: HashMap<String, Element>,
    /// HashMap converting element symbol to element name.
    symbols2names: HashMap<String, String>,
}

/// Contains information about specific element.
#[derive(Debug, Clone)]
#[allow(unused)]
struct Element {
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
        let elements: HashMap<String, Element> = match serde_yaml::from_str(yaml) {
            Ok(x) => x,
            Err(e) => return Err(ParseElementError::CouldNotParseYaml(e)),
        };
        
        let symbols2names = SupportedElements::make_symbols2names(&elements)?;

        Ok(SupportedElements { elements, symbols2names })
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
    /// # use groan_rs::prelude::*;
    /// #
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
        let parsed_elements = SupportedElements::new_from_file(filename)?;

        // check that the merged elements do not have duplicate element symbols
        for (name1, element1) in parsed_elements.elements.iter() {
            for (name2, element2) in self.elements.iter() {
                match (element1.symbol.as_ref(), element2.symbol.as_ref()) {
                    (Some(sym1), Some(sym2)) if sym1 == sym2 && name1 != name2 => {
                        // check for switched
                        let exchange = match parsed_elements.elements.get(name2) {
                            Some(element3) => element3.symbol.as_ref(),
                            None => return Err(ParseElementError::DuplicateSymbol(sym1.to_string(), name1.to_string(), name2.to_string())),
                        };

                        match exchange {
                            Some(sym3) if sym3 != sym1 => (),
                            _ => return Err(ParseElementError::DuplicateSymbol(sym1.to_string(), name1.to_string(), name2.to_string())),
                        }
                    }
                    _ => continue,
                }
            }
        }

        // update the elements
        for (name, element) in parsed_elements.elements.into_iter() {
            match self.elements.get_mut(&name) {
                Some(old) => old.update(element),
                None => {
                    match self.elements.insert(name, element) {
                        None => (),
                        Some(_) => panic!("FATAL GROAN ERROR | SupportedElements::update_from_file | Element should not exist.")
                    }
                }
            }
        }

        // generate new symbols2names map
        self.symbols2names = SupportedElements::make_symbols2names(&self.elements)
            .expect("FATAL GROAN ERROR | SupportedElements::update_from_file | `symbols2names` hashmap should have been created.");

        Ok(())
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

    /// Creates new `symbols2names` hashmap from the provided `elements`.
    fn make_symbols2names(elements: &HashMap<String, Element>) -> Result<HashMap<String, String>, ParseElementError> {
        let mut symbols2names = HashMap::new();

        for (name, element) in elements.iter() {
            if let Some(ref symbol) = element.symbol {
                match symbols2names.insert(symbol.to_string(), name.to_string()) {
                    Some(x) => return Err(ParseElementError::DuplicateSymbol(symbol.to_string(), name.to_string(), x)),
                    None => (),
                }
            }
        }

        Ok(symbols2names)
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

        while let Some(key) = map.next_key::<String>()? {
            match key.as_str() {
                "symbol" => symbol = map.next_value()?,
                "query" => {
                    let query: Option<String> = map.next_value()?;
                    select = match query {
                        Some(q) => match parse_query(&q) {
                            Ok(sel) => Some(*sel),
                            Err(_) => return Err(de::Error::custom("invalid query")),
                        },
                        None => None,
                    };
                },
                "mass" => {
                    let value: Option<f32> = map.next_value()?;
                    mass = match value {
                        Some(v) if v < 0.0 => return Err(de::Error::custom("mass is negative")),
                        _ => value,
                    };
                },
                "vdw" => {
                    let value: Option<f32> = map.next_value()?;
                    vdw = match value {
                        Some(v) if v < 0.0 => return Err(de::Error::custom("vdw is negative")),
                        _ => value,
                    };
                },
                "expected_max_bonds" => expected_max_bonds = map.next_value()?,
                _ => return Err(de::Error::unknown_field(&key, FIELDS)),
            }
        }

        Ok(Element {
            symbol,
            select,
            mass,
            vdw,
            expected_max_bonds,
        })
    }
}

const FIELDS: &'static [&'static str] = &["symbol", "query", "mass", "vdw", "expected_max_bonds"];

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
        let elements = SupportedElements::default();

        assert_eq!(elements.elements.len(), 39);
        assert_eq!(elements.symbols2names.len(), 38); // dummy does not have a symbol

        for name in elements.symbols2names.values() {
            assert!(elements.elements.get(name).is_some());
        }

        // first element
        let e = elements.elements.get("dummy").unwrap();
        assert!(e.symbol.is_none());
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 0.0);
        assert!(e.vdw.is_none());
        assert!(e.expected_max_bonds.is_none());

        // important element
        let e = elements.elements.get("carbon").unwrap();
        assert_eq!(e.symbol.as_ref().unwrap(), "C");
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 12.0107);
        assert_approx_eq!(f32, e.vdw.unwrap(), 0.17);
        assert_eq!(e.expected_max_bonds.unwrap(), 4);

        // last element
        let e = elements.elements.get("bismuth").unwrap();
        assert_eq!(e.symbol.as_ref().unwrap(), "Bi");
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 208.98040);
        assert!(e.vdw.is_none());
        assert!(e.expected_max_bonds.is_none());
    }

    fn elements_match(el1: &mut SupportedElements, el2: &mut SupportedElements) {
        assert_eq!(
            el1.elements.len(),
            el2.elements.len()
        );

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
        let mut default_elements = SupportedElements::default();
        let mut elements_from_file =
            SupportedElements::new_from_file("src/config/elements.yaml").unwrap();

        elements_match(&mut default_elements, &mut elements_from_file);
    }

    #[test]
    fn elements_invalid_field() {
        match SupportedElements::new_from_file("test_files/elements_invalid_field.yaml") {
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
        match SupportedElements::new_from_file("test_files/elements_invalid_mass.yaml") {
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
        match SupportedElements::new_from_file("test_files/elements_invalid_vdw.yaml") {
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
    fn elements_duplicate_symbol() {
        match SupportedElements::new_from_file("test_files/elements_duplicate_symbol.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::DuplicateSymbol(sym, n1, n2)) => {
                assert_eq!(sym, String::from("H"));
                assert!(n1 == String::from("helium") || n1 == String::from("hydrogen"));
                assert!(n2 == String::from("helium") || n2 == String::from("hydrogen"));
                assert_ne!(n1, n2);
            }
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
        
        assert_eq!(elements.elements.len(), 40);
        assert_eq!(elements.symbols2names.len(), 39); // dummy does not have a symbol

        for name in elements.symbols2names.values() {
            assert!(elements.elements.get(name).is_some());
        }

        // first element
        let e = elements.elements.get("dummy").unwrap();
        assert!(e.symbol.is_none());
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 0.0);
        assert!(e.vdw.is_none());
        assert!(e.expected_max_bonds.is_none());

        // changed element #1
        let e = elements.elements.get("carbon").unwrap();
        assert_eq!(e.symbol.as_ref().unwrap(), "C");
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 12.0107);
        assert_approx_eq!(f32, e.vdw.unwrap(), 0.20);
        assert_eq!(e.expected_max_bonds.unwrap(), 4);

        // changed element #2
        let e = elements.elements.get("hydrogen").unwrap();
        assert_eq!(e.symbol.as_ref().unwrap(), "HH");
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 2.014);
        assert_approx_eq!(f32, e.vdw.unwrap(), 0.12);
        assert_eq!(e.expected_max_bonds.unwrap(), 8);

        // added element
        let e = elements.elements.get("polonium").unwrap();
        assert_eq!(e.symbol.as_ref().unwrap(), "Po");
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 209.0);
        assert!(e.vdw.is_none());
        assert!(e.expected_max_bonds.is_none());

        // last element
        let e = elements.elements.get("bismuth").unwrap();
        assert_eq!(e.symbol.as_ref().unwrap(), "Bi");
        assert!(e.select.is_some());
        assert_approx_eq!(f32, e.mass.unwrap(), 208.98040);
        assert!(e.vdw.is_none());
        assert!(e.expected_max_bonds.is_none());

    }

    #[test]
    fn elements_update_invalid_field() {
        let mut elements = SupportedElements::default();

        match elements.update_from_file("test_files/elements_invalid_field.yaml") {
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
            Err(ParseElementError::CouldNotParseYaml(e)) => {
                let string = e.to_string();
                assert!(string.contains("invalid query"));
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

        match elements.update_from_file("test_files/elements_invalid_mass.yaml") {
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

        let mut default = SupportedElements::default();
        elements_match(&mut elements, &mut default);
    }

    #[test]
    fn elements_update_invalid_vdw() {
        let mut elements = SupportedElements::default();

        match elements.update_from_file("test_files/elements_invalid_vdw.yaml") {
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

        let mut default = SupportedElements::default();
        elements_match(&mut elements, &mut default);
    }

    #[test]
    fn elements_update_duplicate_symbol() {
        let mut elements = SupportedElements::default();

        match elements.update_from_file("test_files/elements_duplicate_symbol.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::DuplicateSymbol(sym, n1, n2)) => {
                assert_eq!(sym, String::from("H"));
                assert!(n1 == String::from("helium") || n1 == String::from("hydrogen"));
                assert!(n2 == String::from("helium") || n2 == String::from("hydrogen"));
                assert_ne!(n1, n2);
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
    fn elements_update_duplicate_symbol_new() {
        let mut elements = SupportedElements::default();

        match elements.update_from_file("test_files/elements_duplicate_symbol3.yaml") {
            Ok(_) => panic!("Function should have failed but it succeeded."),
            Err(ParseElementError::DuplicateSymbol(sym, n1, n2)) => {
                assert_eq!(sym, String::from("H"));
                assert!(n1 == String::from("polonium") || n1 == String::from("hydrogen"));
                assert!(n2 == String::from("polonium") || n2 == String::from("hydrogen"));
                assert_ne!(n1, n2);
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
    fn elements_update_switched_symbols() {
        let mut elements = SupportedElements::default();

        elements.update_from_file("test_files/elements_symbols_switched.yaml").unwrap();

        let e = elements.elements.get("carbon").unwrap();
        assert_eq!(e.symbol.as_deref().unwrap(), "H");
        let e = elements.elements.get("hydrogen").unwrap();
        assert_eq!(e.symbol.as_deref().unwrap(), "C");
        // unchanged
        let e = elements.elements.get("bismuth").unwrap();
        assert_eq!(e.symbol.as_deref().unwrap(), "Bi");
    }

    #[test]
    fn elements_update_complex_switched_symbols() {
        let mut elements = SupportedElements::default();

        elements.update_from_file("test_files/elements_complex_switches.yaml").unwrap();

        let e = elements.elements.get("carbon").unwrap();
        assert_eq!(e.symbol.as_deref().unwrap(), "N");
        let e = elements.elements.get("hydrogen").unwrap();
        assert_eq!(e.symbol.as_deref().unwrap(), "C");
        let e = elements.elements.get("nitrogen").unwrap();
        assert_eq!(e.symbol.as_deref().unwrap(), "P");
        let e = elements.elements.get("phosphorus").unwrap();
        assert_eq!(e.symbol.as_deref().unwrap(), "H");
        // unchanged
        let e = elements.elements.get("bismuth").unwrap();
        assert_eq!(e.symbol.as_deref().unwrap(), "Bi");
    }
}