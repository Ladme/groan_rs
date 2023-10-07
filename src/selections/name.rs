// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of the Name enum for the Groan selection language.

use crate::errors::SelectError;
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
            let regex = match Regex::new(&(string[2..string.len() - 1].to_string())) {
                Ok(reg) => reg,
                Err(_) => return Err(SelectError::InvalidRegex(string.to_owned())),
            };

            Ok(Name::Regex(regex))
        } else {
            Ok(Name::String(string.to_owned()))
        }
    }

    pub fn to_str(&self) -> &str {
        match self {
            Name::String(s) => s,
            Name::Regex(_) => todo!("Regex is currently not implemented for Groups."),
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
