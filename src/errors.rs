// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of errors that can be returned by the library.

use colored::{ColoredString, Colorize};
use std::collections::HashSet;
use std::path::Path;
use thiserror::Error;

fn path_to_yellow(path: &Box<Path>) -> ColoredString {
    path.to_str().unwrap().yellow()
}

fn unpack_set(set: &Box<HashSet<String>>) -> ColoredString {
    let mut output = String::new();
    let len = set.len();
    output.push_str("\n");

    for (i, key) in set.iter().enumerate() {
        output.push_str(&key);
        output.push_str("\n");

        if i >= 9 && len != i + 1 {
            let and_more = format!("...and {} more...", len - i - 1);
            output.push_str(&and_more);
            break;
        }
    }

    output.yellow()
}

/// Errors that can occur when reading a file of unknown type.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseFileError {
    #[error("{} file '{}' has an unknown file extension", "error:".red().bold(), path_to_yellow(.0))]
    UnknownExtension(Box<Path>),
}

/// Errors that can occur when reading and parsing gro file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseGroError {
    #[error("{} file '{}' was not found", "error:".red().bold(), path_to_yellow(.0))]
    FileNotFound(Box<Path>),
    #[error("{} file '{}' ended unexpectedly", "error:".red().bold(), path_to_yellow(.0))]
    LineNotFound(Box<Path>),
    #[error("{} could not parse line '{}'", "error:".red().bold(), .0.yellow())]
    ParseLineErr(String),
    #[error("{} could not parse line '{}' as atom", "error:".red().bold(), .0.yellow())]
    ParseAtomLineErr(String),
    #[error("{} could not parse line '{}' as box dimensions", "error:".red().bold(), .0.yellow())]
    ParseBoxLineErr(String),
}

/// Errors that can occur when reading and parsing pdb file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParsePdbError {
    #[error("{} file '{}' was not found", "error:".red().bold(), path_to_yellow(.0))]
    FileNotFound(Box<Path>),
    #[error("{} file '{}' ended unexpectedly", "error:".red().bold(), path_to_yellow(.0))]
    LineNotFound(Box<Path>),
    #[error("{} could not parse line '{}'", "error:".red().bold(), .0.yellow())]
    ParseLineErr(String),
    #[error("{} could not parse line '{}' as atom", "error:".red().bold(), .0.yellow())]
    ParseAtomLineErr(String),
    #[error("{} could not parse line '{}' as box dimensions", "error:".red().bold(), .0.yellow())]
    ParseBoxLineErr(String),
    #[error("{} could not parse line '{}' as title", "error:".red().bold(), .0.yellow())]
    ParseTitleLineErr(String),
    #[error("{} simulation box specified on line '{}' is not orthogonal", "error:".red().bold(), .0.yellow())]
    NonOrthogonalBox(String),
}

/// Errors that can occur when writing a gro file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum WriteGroError {
    #[error("{} file '{}' could not be created", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotCreate(Box<Path>),
    #[error("{} could not write line into file", "error:".red().bold())]
    CouldNotWrite(),
    #[error("{} group '{}' does not exist", "error:".red().bold(), .0.yellow())]
    GroupNotFound(String),
}

/// Errors that can occur when working with Groups of atoms.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum GroupError {
    #[error("{} group '{}' does not exist", "error:".red().bold(), .0.yellow())]
    NotFound(String),
    #[error("{} group '{}' already exists", "warning:".yellow().bold(), .0.yellow())]
    AlreadyExistsWarning(String),
    #[error("{} following groups already exist and have been overwritten: {}", "warning:".yellow().bold(), unpack_set(.0))]
    MultipleAlreadyExistWarning(Box<HashSet<String>>),
    #[error("{} name '{}' contains invalid characters", "error:".red().bold(), .0.yellow())]
    InvalidName(String),
}

/// Errors that can occur when reading and parsing ndx file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseNdxError {
    #[error("{} file '{}' was not found", "error:".red().bold(), path_to_yellow(.0))]
    FileNotFound(Box<Path>),
    #[error("{} file '{}' ended unexpectedly", "error:".red().bold(), path_to_yellow(.0))]
    LineNotFound(Box<Path>),
    #[error("{} could not parse line '{}' as group name", "error:".red().bold(), .0.yellow())]
    ParseGroupNameErr(String),
    #[error("{} could not parse line '{}'", "error:".red().bold(), .0.yellow())]
    ParseLineErr(String),
    #[error("{} atom index '{}' does not exist in the system", "error:".red().bold(), .0.to_string().yellow())]
    InvalidAtomIndex(usize),
    #[error("{} duplicate groups detected: {}", "warning:".yellow().bold(), unpack_set(.0))]
    DuplicateGroupsWarning(Box<HashSet<String>>),
    #[error("{} ignored the following groups with invalid names: {}", "warning:".yellow().bold(), unpack_set(.0))]
    InvalidNamesWarning(Box<HashSet<String>>),
}

/// Errors that can occur when writing an ndx file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum WriteNdxError {
    #[error("{} file '{}' could not be created", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotCreate(Box<Path>),
    #[error("{} could not write line into file", "error:".red().bold())]
    CouldNotWrite,
}

/// Errors that can occur when working with xtc file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum XtcError {
    #[error("{} unable to work with path '{}'", "error:".red().bold(), path_to_yellow(.0))]
    InvalidPath(Box<Path>),
    #[error("{} file '{}' could not be opened or created", "error:".red().bold(), path_to_yellow(.0))]
    FileNotFound(Box<Path>),
}

/// Errors that can occur when reading an xtc file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ReadXtcError {
    #[error("{} unable to work with path '{}'", "error:".red().bold(), path_to_yellow(.0))]
    InvalidPath(Box<Path>),
    #[error("{} file '{}' was not found or could not be read as an xtc file", "error:".red().bold(), path_to_yellow(.0))]
    FileNotFound(Box<Path>),
    #[error("{} could not read frame in an xtc file", "error:".red().bold())]
    FrameNotFound,
    #[error("{} number of atoms in the xtc file '{}' does not match the number of atoms in the system", "error:".red().bold(), path_to_yellow(.0))]
    AtomsNumberMismatch(Box<Path>),
}

/// Errors that can occur when writing an xtc file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum WriteXtcError {
    #[error("{} unable to work with path '{}'", "error:".red().bold(), path_to_yellow(.0))]
    InvalidPath(Box<Path>),
    #[error("{} file '{}' could not be created", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotCreate(Box<Path>),
    #[error("{} could not write frame to an xtc file", "error:".red().bold())]
    CouldNotWrite,
}

/// Errors that can occur when parsing atom selection query.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum SelectError {
    #[error("{} the provided query is empty", "error:".red().bold())]
    EmptyQuery,
    #[error("{} invalid operator detected in query '{}'", "error:".red().bold(), .0.to_string().yellow())]
    InvalidOperator(String),
    #[error("{} missing argument in query '{}'", "error:".red().bold(), .0.to_string().yellow())]
    MissingArgument(String),
    #[error("{} empty argument in query '{}'", "error:".red().bold(), .0.to_string().yellow())]
    EmptyArgument(String),
    #[error("{} unmatching parentheses in query '{}'", "error:".red().bold(), .0.to_string().yellow())]
    InvalidParentheses(String),
    #[error("{} unmatching number of quotes in query '{}'", "error:".red().bold(), .0.to_string().yellow())]
    InvalidQuotes(String),
    #[error("{} could not understand the residue/atom numbers in query '{}'", "error:".red().bold(), .0.to_string().yellow())]
    InvalidNumber(String),
    #[error("{} group '{}' does not exist", "error:".red().bold(), .0.to_string().yellow())]
    GroupNotFound(String),
    #[error("{} the provided query '{}' could not be understood for unknown reason", "error:".red().bold(), .0.to_string().yellow())]
    UnknownError(String),
}
