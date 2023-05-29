// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

use thiserror::Error;
use std::path::Path;

/// Errors that can occur when reading and parsing gro file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseGroError {
    #[error("File `{0}` was not found.")]
    FileNotFoundErr(Box<Path>),
    #[error("File `{0}` ended unexpectedly.")]
    LineNotFoundErr(Box<Path>),
    #[error("Could not parse line `{0}`.")]
    ParseLineErr(String),
    #[error("Could not parse line `{0}` as atom.")]
    ParseAtomLineErr(String),
    #[error("Could not parse line `{0}` as box dimensions.")]
    ParseBoxLineErr(String),
}

#[derive(Error, Debug, PartialEq, Eq)]
pub enum GroupError {
    #[error("Group `{0}` does not exist.")]
    NotFound(String),
    #[error("Group `{0}` already exists.")]
    AlreadyExists(String),
}
