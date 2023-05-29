// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

use thiserror::Error;
use std::path::Path;

/// Errors that can occur when reading and parsing gro file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseGroError {
    #[error("File `{0}` was not found.")]
    FileNotFound(Box<Path>),
    #[error("File `{0}` ended unexpectedly.")]
    LineNotFound(Box<Path>),
    #[error("Could not parse line `{0}`.")]
    ParseLineErr(String),
    #[error("Could not parse line `{0}` as atom.")]
    ParseAtomLineErr(String),
    #[error("Could not parse line `{0}` as box dimensions.")]
    ParseBoxLineErr(String),
}


/// Errors that can occur when working with atom Groups.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum GroupError {
    #[error("Group `{0}` does not exist.")]
    NotFound(String),
    #[error("Group `{0}` already exists.")]
    AlreadyExists(String),
}

/// Errors that can occur when reading and parsing ndx file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseNdxError {
    #[error("File `{0} was not found.")]
    FileNotFound(Box<Path>),
    #[error("File `{0}` ended unexpectedly.")]
    LineNotFound(Box<Path>),
    #[error("Could not parse line `{0}` as group name.")]
    ParseGroupNameErr(String),
    #[error("Could not parse line `{0}`.")]
    ParseLineErr(String),
    #[error("Group `{0}` already exists in the system.")]
    GroupAlreadyExists(String),
    #[error("The ndx file contains multiple groups names `{0}`.")]
    GroupsShareName(String),
    #[error("Atom index `{0}` does not exist in the system.")]
    InvalidAtomIndex(u64),
}

#[derive(Error, Debug, PartialEq, Eq)]
pub enum WriteGroError {
    
}