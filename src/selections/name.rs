// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of the Name enum for the Groan selection language.

use crate::errors::SelectError;
use crate::system::general::System;
use regex::Regex;

#[derive(Debug)]
pub enum Name {
    String(String),
    Regex(Regex),
}

impl Name {
    /// Create new `Name` enum. `Name` enum contains either a String or a Regex.
    pub fn new(string: &str) -> Result<Self, SelectError> {
        if string.starts_with("r'") {
            let regex = match Regex::new(&(string[2..string.len() - 1])) {
                Ok(reg) => reg,
                Err(_) => return Err(SelectError::InvalidRegex(string.to_owned())),
            };

            Ok(Name::Regex(regex))
        } else {
            Ok(Name::String(string.to_owned()))
        }
    }
    
    /// Check whether `atom_index` is part of the specified group.
    pub fn match_groups(&self, system: &System, atom_index: usize) -> Result<bool, SelectError> {
        match self {
            Name::String(s) => match system.group_isin(s, atom_index) {
                Err(e) => Err(SelectError::GroupNotFound(e.to_string())),
                Ok(result) => Ok(result),
            },

            // this is inefficient; it's better to expand the regular expression into a vector of group names
            // matching the regex; this is what happens when `Group::from_select` is called;
            // this match arm is here just in case this is forgotten
            Name::Regex(r) => {
                let group_names = system.get_groups_as_ref().keys();

                for name in group_names {
                    if r.is_match(name) {
                        return Ok(system
                            .group_isin(name, atom_index)
                            .expect("Groan error. `match_groups` regex: group must exist."));
                    }
                }
                Ok(false)
            }
        }
    }

    /// Converts `Name` to string without consuming `Name`.
    pub fn to_string(&self) -> String {
        match self {
            Name::String(s) => s.clone(),
            Name::Regex(r) => r.to_string(),
        }
    }
}

impl PartialEq<Name> for Name {
    /// Compare `Name` enum with another `Name` enum.
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Name::String(s), Name::String(t)) => s == t,
            (Name::String(s), Name::Regex(t)) => s == t.as_str(),
            (Name::Regex(s), Name::String(t)) => s.as_str() == t,
            (Name::Regex(s), Name::Regex(t)) => s.as_str() == t.as_str(),
        }
    }
}

impl PartialEq<str> for Name {
    /// Compare `Name` enum with `&str`.
    fn eq(&self, other: &str) -> bool {
        match self {
            Name::String(s) => s == other,
            Name::Regex(s) => s.is_match(other),
        }
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use super::*;
    use crate::selections::select::Select;

    #[test]
    fn match_groups() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let selection = Select::GroupName(vec![
            Name::new("Membrane").unwrap(),
            Name::new("r'membrane'").unwrap(),
        ]);

        match selection {
            Select::GroupName(v) => {
                let membrane = &v[0];
                for atom_index in 0..system.get_n_atoms() {
                    let result = membrane.match_groups(&system, atom_index).unwrap();
                    if atom_index >= 61 && atom_index <= 6204 {
                        assert!(result);
                    } else {
                        assert!(!result);
                    }
                }

                let regex = &v[1];
                for atom_index in 0..system.get_n_atoms() {
                    let result = regex.match_groups(&system, atom_index).unwrap();
                    if atom_index < 61 {
                        assert!(result);
                    } else {
                        assert!(!result);
                    }
                }
            }
            _ => panic!("This can't happen."),
        }
    }
}
