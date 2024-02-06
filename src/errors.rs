// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of errors that can be returned by the library.

use colored::{ColoredString, Colorize};
use std::collections::HashSet;
use std::path::Path;
use thiserror::Error;

use crate::system::guess::{BondsGuessInfo, ElementGuessInfo, PropertiesGuessInfo};

fn path_to_yellow(path: &Path) -> ColoredString {
    path.to_str().unwrap().yellow()
}

fn unpack_set(set: &HashSet<String>) -> ColoredString {
    let mut output = String::new();
    let len = set.len();
    output.push('\n');

    for (i, key) in set.iter().enumerate() {
        output.push_str(key);
        output.push('\n');

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
    /// Used when a file has an unknown or unsupported file extension or is missing the file extension entirely
    /// and therefore the file type/format can not be identified.
    #[error("{} file '{}' has an unknown or unsupported file extension", "error:".red().bold(), path_to_yellow(.0))]
    UnknownExtension(Box<Path>),
}

/// Errors that can occur when reading and parsing gro file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseGroError {
    /// Used when the gro file was not found (i.e. does not exist).
    #[error("{} file '{}' was not found", "error:".red().bold(), path_to_yellow(.0))]
    FileNotFound(Box<Path>),
    /// Used when the gro file ended unexpectedly.
    #[error("{} file '{}' ended unexpectedly", "error:".red().bold(), path_to_yellow(.0))]
    LineNotFound(Box<Path>),
    /// Used when a generic line in the gro file could not be parsed.
    #[error("{} could not parse line '{}'", "error:".red().bold(), .0.yellow())]
    ParseLineErr(String),
    /// Used when an "atom line" in the gro file could not be parsed.
    #[error("{} could not parse line '{}' as atom", "error:".red().bold(), .0.yellow())]
    ParseAtomLineErr(String),
    /// Used when a "box line" in the gro file could not be parsed.
    #[error("{} could not parse line '{}' as box dimensions", "error:".red().bold(), .0.yellow())]
    ParseBoxLineErr(String),
}

/// Errors that can occur when reading and parsing pdb file.
/// Does not include errors that can occur when reading connectivity section of a PDB file.
/// For these errors, see `ParsePdbConnectivityError`.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParsePdbError {
    /// Used when the pdb file was not found (i.e. does not exist).
    #[error("{} file '{}' was not found", "error:".red().bold(), path_to_yellow(.0))]
    FileNotFound(Box<Path>),
    /// Used when the pdb file ended unexpectedly.
    #[error("{} file '{}' ended unexpectedly", "error:".red().bold(), path_to_yellow(.0))]
    LineNotFound(Box<Path>),
    /// Used when a generic line in the pdb file could not be parsed.
    #[error("{} could not parse line '{}'", "error:".red().bold(), .0.yellow())]
    ParseLineErr(String),
    /// Used when an "ATOM" or "HETATM" line in the pdb file could not be parsed.
    #[error("{} could not parse line '{}' as atom", "error:".red().bold(), .0.yellow())]
    ParseAtomLineErr(String),
    /// Used when a "CRYST1" line in the pdb file could not be parsed.
    #[error("{} could not parse line '{}' as box dimensions", "error:".red().bold(), .0.yellow())]
    ParseBoxLineErr(String),
    /// Used when a "TITLE" line in the pdb file could not be parsed.
    #[error("{} could not parse line '{}' as title", "error:".red().bold(), .0.yellow())]
    ParseTitleLineErr(String),
}

/// Errors that can occur when reading the connectivity section of a PDB file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParsePdbConnectivityError {
    /// Used when the pdb file was not found (i.e. does not exist).
    #[error("{} file '{}' was not found", "error:".red().bold(), path_to_yellow(.0))]
    FileNotFound(Box<Path>),
    /// Used when the pdb file ended unexpectedly.
    #[error("{} file '{}' ended unexpectedly", "error:".red().bold(), path_to_yellow(.0))]
    LineNotFound(Box<Path>),
    /// Used when "CONECT" line in the pdb file could not be parsed.
    #[error("{} could not parse line '{}' as connectivity information", "error:".red().bold(), .0.yellow())]
    ParseConectLineErr(String),
    /// Used when an non-existent atom number is found on a "CONECT" line.
    #[error("{} atom number '{}' mentioned on line '{}' does not exist", "error:".red().bold(), .0.to_string().yellow(), .1.yellow())]
    AtomNotFound(usize, String),
    /// Used when there are multiple atoms with the same number in the PDB file.
    #[error("{} multiple atoms have the same number in the system and connectivity is thus ambiguous", "error:".red().bold())]
    DuplicateAtomNumbers,
    /// Used when an atom claims to be bonded to itself.
    #[error("{} atom '{}' claims to be bonded to itself which does not make sense", "error:".red().bold(), .0.to_string().yellow())]
    SelfBonding(usize),
    /// Used when no connectivity information has been read from the PDB file.
    #[error("{} no bonds have been found in the PDB file '{}'", "warning:".yellow().bold(), path_to_yellow(.0))]
    NoBondsWarning(Box<Path>),
}

/// Errors that can occur when writing a gro file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum WriteGroError {
    /// Used when the gro file could not be opened for writing (i.e. path is invalid).
    #[error("{} file '{}' could not be created", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotCreate(Box<Path>),
    /// Used when writing into the gro file failed for any reason.
    #[error("{} could not write line into file", "error:".red().bold())]
    CouldNotWrite,
    /// Used when the group of atoms selected to be written into the gro file does not exist.
    #[error("{} group '{}' does not exist", "error:".red().bold(), .0.yellow())]
    GroupNotFound(String),
}

/// Errors that can occur when writing a pdb file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum WritePdbError {
    /// Used when the pdb file could not be opened for writing (i.e. path is invalid).
    #[error("{} file '{}' could not be created", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotCreate(Box<Path>),
    /// Used when writing into the pdb file failed for any reason.
    #[error("{} could not write line into file", "error:".red().bold())]
    CouldNotWrite,
    /// Used when the group of atoms selected to be written into the pdb file does not exist.
    #[error("{} group '{}' does not exist", "error:".red().bold(), .0.yellow())]
    GroupNotFound(String),
    /// Used when connectivity printing is requested and the system is too large for the PDB file.
    #[error("{} system is too large ('{}' atoms) for PDB connectivity section", "error:".red().bold(), .0.to_string().yellow())]
    ConectTooLarge(usize),
    /// Used when there are multiple atoms with the same number in the system and connectivity thus can't be printed.
    #[error("{} multiple atoms have the same number in the system and connectivity is thus ambiguous", "error:".red().bold())]
    ConectDuplicateAtomNumbers,
    /// Used when the atom number to be printed in the connectivity section is higher than 99,999.
    #[error("{} atom number '{}' is too high for PDB connectivity section and can not be wrapped", "error:".red().bold(), .0.to_string().yellow())]
    ConectInvalidNumber(usize),
}

/// Errors that can occur when reading and parsing pqr file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParsePqrError {
    /// Used when the pqr file was not found (i.e. does not exist).
    #[error("{} file '{}' was not found", "error:".red().bold(), path_to_yellow(.0))]
    FileNotFound(Box<Path>),
    /// Used when the pqr file ended unexpectedly.
    #[error("{} file '{}' ended unexpectedly", "error:".red().bold(), path_to_yellow(.0))]
    LineNotFound(Box<Path>),
    /// Used when a generic line in the pqr file could not be parsed.
    #[error("{} could not parse line '{}'", "error:".red().bold(), .0.yellow())]
    ParseLineErr(String),
    /// Used when an "ATOM" or "HETATM" line in the pqr file could not be parsed.
    #[error("{} could not parse line '{}' as atom", "error:".red().bold(), .0.yellow())]
    ParseAtomLineErr(String),
}

/// Errors that can occur when working with Groups of atoms.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum GroupError {
    /// Used when the specified group of atoms does not exist.
    #[error("{} group '{}' does not exist", "error:".red().bold(), .0.yellow())]
    NotFound(String),
    /// Used when the group with the same name already exists. This is a warning and does not indicate failure.
    #[error("{} group '{}' already existed and has been overwritten", "warning:".yellow().bold(), .0.yellow())]
    AlreadyExistsWarning(String),
    /// Used when at least one group shares its name with any of the newly created groups. This is a warning and does not indicate failure.
    /// This warning is used when creating multiple groups of atoms in a single function. Otherwise, `GroupError::AlreadyExistsWarning` is used.
    #[error("{} following groups already existed and have been overwritten: {}", "warning:".yellow().bold(), unpack_set(.0))]
    MultipleAlreadyExistWarning(Box<HashSet<String>>),
    /// Used when the name of the group contains invalid character(s).
    #[error("{} name '{}' contains invalid characters", "error:".red().bold(), .0.yellow())]
    InvalidName(String),
    /// Used when the groan selection language query provided to create the group is invalid.
    /// Encapsulates the `SelectError` providing more information about the type of error.
    #[error("{}", .0)]
    InvalidQuery(SelectError),
    /// Used when there is an issue with simulation box of the system.
    #[error("{}", .0)]
    InvalidSimBox(SimBoxError),
    /// Used when there is an issue with positions of atoms in the system.
    #[error("{}", .0)]
    InvalidPosition(PositionError),
    /// Used when there is an issue with masses of atoms in the system.
    #[error("{}", .0)]
    InvalidMass(MassError),
    /// Used when the group is empty and is expected not to be.
    #[error("{} group '{}' contains no atoms", "error:".red().bold(), .0.yellow())]
    EmptyGroup(String),
}

/// Errors that can occur when working with atoms in a system.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum AtomError {
    /// Used when selecting an atom from the system with invalid atom index.
    #[error("{} atom index '{}' is out of range", "error:".red().bold(), .0.to_string().yellow())]
    OutOfRange(usize),
    /// Used when attempting to create a bond that is invalid.
    #[error("{} bond could not be created between atoms {} and {}", "error:".red().bold(), .0.to_string().yellow(), .1.to_string().yellow())]
    InvalidBond(usize, usize),
    /// Used when there is an issue with simulation box of the system.
    #[error("{}", .0)]
    InvalidSimBox(SimBoxError),
    /// Used when there is an issue with position of an atom.
    #[error("{}", .0)]
    InvalidPosition(PositionError),
}

/// Errors that can occur when reading and parsing ndx file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ParseNdxError {
    /// Used when the ndx file was not found (i.e. does not exist).
    #[error("{} file '{}' was not found", "error:".red().bold(), path_to_yellow(.0))]
    FileNotFound(Box<Path>),
    /// Used when the ndx file ended unexpectedly.
    #[error("{} file '{}' ended unexpectedly", "error:".red().bold(), path_to_yellow(.0))]
    LineNotFound(Box<Path>),
    /// Used when line containing the name of the ndx group could not be parsed.
    #[error("{} could not parse line '{}' as group name", "error:".red().bold(), .0.yellow())]
    ParseGroupNameErr(String),
    /// Used when a generic line in the ndx file could not be parsed.
    #[error("{} could not parse line '{}'", "error:".red().bold(), .0.yellow())]
    ParseLineErr(String),
    /// Used when an atom index from the ndx file does not correspond to any atom in the system.
    #[error("{} atom index '{}' does not exist in the system", "error:".red().bold(), .0.to_string().yellow())]
    InvalidAtomIndex(usize),
    /// Used when at least one group shares its name with any group read from the ndx file. This is a warning and does not indicate failure.
    #[error("{} duplicate groups detected: {}", "warning:".yellow().bold(), unpack_set(.0))]
    DuplicateGroupsWarning(Box<HashSet<String>>),
    /// Used when at least one group read from the ndx file has an invalid name. Note that even though this is a warning, it indicates
    /// partial failure of the function as groups with invalid names are NOT CREATED.
    #[error("{} ignored the following groups with invalid names: {}", "warning:".yellow().bold(), unpack_set(.0))]
    InvalidNamesWarning(Box<HashSet<String>>),
}

/// Errors that can occur when writing an ndx file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum WriteNdxError {
    /// Used when the ndx file could not be opened for writing (i.e. path is invalid).
    #[error("{} file '{}' could not be created", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotCreate(Box<Path>),
    /// Used when writing into the ndx file failed for any reason.
    #[error("{} could not write line into file", "error:".red().bold())]
    CouldNotWrite,
}

/// Errors that can occur when working with a trajectory file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum TrajError {
    /// Used when the path to the trajectory file is invalid (i.e. contains invalid characters).
    #[error("{} unable to work with path '{}'", "error:".red().bold(), path_to_yellow(.0))]
    InvalidPath(Box<Path>),
    /// Used when the specified trajectory file could not be reached.
    #[error("{} file '{}' could not be opened or created", "error:".red().bold(), path_to_yellow(.0))]
    FileNotFound(Box<Path>),
}

/// Errors that can occur when reading a trajectory file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ReadTrajError {
    /// Used when the path to the trajectory file is invalid (i.e. contains invalid characters).
    #[error("{} unable to work with path '{}'", "error:".red().bold(), path_to_yellow(.0))]
    InvalidPath(Box<Path>),
    /// Used when the trajectory file does not exist, could not be read or is not a valid trajectory file.
    #[error("{} file '{}' was not found or could not be read as a trajectory file", "error:".red().bold(), path_to_yellow(.0))]
    FileNotFound(Box<Path>),
    /// Used when a frame could not be read from a trajectory file.
    #[error("{} could not read frame in a trajectory file", "error:".red().bold())]
    FrameNotFound,
    /// Used when the number of atoms in the trajectory file does not match the number of atoms in the corresponding `System` structure.
    #[error("{} number of atoms in the trajectory file '{}' does not match the number of atoms in the system", "error:".red().bold(), path_to_yellow(.0))]
    AtomsNumberMismatch(Box<Path>),
    /// Used when the time provided as the start of the time range is higher than the time provided as the end of the time range.
    #[error("{} invalid time range (starting time '{}' ps is higher than the ending time '{}' ps)", "error:".red().bold(), .0.yellow(), .1.yellow())]
    InvalidTimeRange(String, String),
    /// Used when the start time is higher than the time of all the frames in the trajectory file
    #[error("{} start time ('{}' ps) exceeds the time of all frames in the trajectory file", "error:".red().bold(), .0.yellow())]
    StartNotFound(String),
    /// Used when the time provided as the start/end of the time range is negative.
    #[error("{} negative time ('{}' ps) is not allowed in a time range", "error:".red().bold(), .0.yellow())]
    TimeRangeNegative(String),
    /// Used when a frame of the trajectory file could not be skipped over.
    #[error("{} could not skip a frame in a trajectory file", "error:".red().bold())]
    SkipFailed,
    /// Used when the step of the iteration is invalid, usually zero.
    #[error("{} unsupported iteration step '{}' - must be > 0", "error:".red().bold(), .0.to_string().yellow())]
    InvalidStep(usize),
    /// Used when simulation box read from the trajectory is invalid.
    #[error("{} simulation box is invalid", "error:".red().bold())]
    InvalidSimBox,
}

/// Errors that can occur when writing a trajectory file.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum WriteTrajError {
    /// Used when the path to the trajectory file is invalid (i.e. contains invalid characters).
    #[error("{} unable to work with path '{}'", "error:".red().bold(), path_to_yellow(.0))]
    InvalidPath(Box<Path>),
    /// Used when the trajectory file could not be opened for writing (i.e. path is invalid).
    #[error("{} file '{}' could not be created", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotCreate(Box<Path>),
    /// Used when a frame could not be written into a trajectory file for any reason.
    #[error("{} could not write frame to a trajectory file", "error:".red().bold())]
    CouldNotWrite,
    /// Used when the group of atoms selected to be written into the trajectory file does not exist.
    #[error("{} group '{}' does not exist", "error:".red().bold(), .0.yellow())]
    GroupNotFound(String),
}

/// Errors that can occur when parsing atom selection query.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum SelectError {
    /// Used when the provided groan selection language query is empty.
    #[error("{} the provided query is empty", "error:".red().bold())]
    EmptyQuery,
    /// Used when an unknown operator is detected in the groan selection language query.
    #[error("{} invalid operator detected in query '{}'", "error:".red().bold(), .0.to_string().yellow())]
    InvalidOperator(String),
    /// Used when an operator is missing an argument in the groan selection language query.
    #[error("{} missing argument in query '{}'", "error:".red().bold(), .0.to_string().yellow())]
    MissingArgument(String),
    /// Used when a keyword is missing an argument in the groan selection language query.
    #[error("{} empty argument in query '{}'", "error:".red().bold(), .0.to_string().yellow())]
    EmptyArgument(String),
    /// Used when the parentheses are incorrectly used in the groan selection language query.
    #[error("{} unmatching parentheses in query '{}'", "error:".red().bold(), .0.to_string().yellow())]
    InvalidParentheses(String),
    /// Used when the quotes are incorrectly used in the groan selection language query.
    #[error("{} unmatching number of quotes in query '{}'", "error:".red().bold(), .0.to_string().yellow())]
    InvalidQuotes(String),
    /// Used when a `()` expression is not followed or preceeded by a binary operator.
    /// (e.g. things like `(name CA CB) Protein` or `element (name CA)`)
    #[error("{} invalid token following or preceeding parentheses in query '{}'", "error:".red().bold(), .0.to_string().yellow())]
    InvalidTokenParentheses(String),
    /// Used when any error occurs while parsing atom/residue numbers or ranges in the groan selection language query.
    #[error("{} could not understand the residue/atom numbers in query '{}'", "error:".red().bold(), .0.to_string().yellow())]
    InvalidNumber(String),
    /// Used when a group specified in the groan selection language query does not exist.
    #[error("{} group '{}' does not exist", "error:".red().bold(), .0.to_string().yellow())]
    GroupNotFound(String),
    /// Used when an invalid identifier of a chain (i.e. longer than one character) is used in the groan selection language query.
    #[error("{} invalid chain identifier(s) in query '{}'", "error:".red().bold(), .0.to_string().yellow())]
    InvalidChainId(String),
    /// Used when the groan selection language query contains a regular expression that is invalid.
    #[error("{} string '{}' is not a valid regular expression", "error:".red().bold(), .0.to_string().yellow())]
    InvalidRegex(String),
    /// Used when the regular expression is used to select groups but corresponds to no groups in the system.
    /// This is currently only used when no regular expression in the entire subquery corresponds to any groups of atoms.
    #[error("{} regular expression '{}' matches no atom groups in the system", "error:".red().bold(), .0.to_string().yellow())]
    NoRegexMatch(String),
    /// Used when an unknown error which does not have a specific `SelectError` variant occurs while parsing the groan selection language query.
    #[error("{} the provided query '{}' could not be understood for unknown reason", "error:".red().bold(), .0.to_string().yellow())]
    UnknownError(String),
}

/// Errors that can occur when reading element data.
#[derive(Error, Debug)]
pub enum ParseElementError {
    /// Used when the yaml file was not found (i.e. does not exist).
    #[error("{} file '{}' was not found", "error:".red().bold(), path_to_yellow(.0))]
    FileNotFound(Box<Path>),
    /// Used when the content of the yaml file could not be read.
    #[error("{} file '{}' could not be read", "error:".red().bold(), path_to_yellow(.0))]
    FileCouldNotBeRead(Box<Path>),
    /// Used when the yaml string containing element data can't be parsed.
    #[error("{} could not parse yaml input as elements ({})", "error:".red().bold(), .0.to_string().yellow())]
    CouldNotParseYaml(serde_yaml::Error),
    /// Used when the same element symbol corresponds to multiple elements.
    #[error("{} element symbol '{}' corresponds to both '{}' and '{}'", "error:".red().bold(), .0.yellow(), .1.yellow(), .2.yellow())]
    DuplicateSymbol(String, String, String),
}

/// Errors that can occur when working with elements.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum ElementError {
    /// Used in `System::guess_elements`.
    /// Used when a select tree associated with the given element is invalid.
    /// This should only happen when the query uses group names.
    #[error("{} (this error occured when guessing elements)", .0)]
    InvalidQuery(SelectError),

    /// Used in `System::guess_elements`.
    /// Used when there is at least one atom which was not assigned an element
    /// or when there is at least one atom which matches multiple elements.
    /// This is a warning and does not indicate failure of the function.
    #[error("{} when guessing elements, following concerns have been raised:\n{}", "warning:".yellow().bold(), .0.to_string())]
    ElementGuessWarning(Box<ElementGuessInfo>),

    /// Used in `System::guess_properties`.
    /// Used when there is at least one atom which was not assigned all the properties.
    /// This is a warning and does not indicate failure of the function.
    #[error("{} when guessing properties, following concerns have been raised:\n{}", "warning:".yellow().bold(), .0.to_string())]
    PropertiesGuessWarning(Box<PropertiesGuessInfo>),

    /// Used in `System::guess_bonds`.
    /// Used when there is at least one atom with suspicious number of bonds.
    /// This is a warning and does not indicate failure of the function.
    #[error("{} when guessing bonds, following concerns have been raised:\n{}", "warning:".yellow().bold(), .0.to_string())]
    BondsGuessWarning(Box<BondsGuessInfo>),
}

/// Errors that can occur when working with simulation box.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum SimBoxError {
    /// Used when the system has no defined simulation box.
    #[error("{} system has no simulation box", "error:".red().bold())]
    DoesNotExist,
    /// Used when the simulation box is not orthogonal but is required to be.
    #[error("{} simulation box is not orthogonal but is required to be", "error:".red().bold())]
    NotOrthogonal,
}

/// Errors that can occur when working with positions of atoms.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum PositionError {
    /// Used when the atom has no position but it is required.
    #[error("{} atom with atom number '{}' has undefined position", "error:".red().bold(), .0.to_string().yellow())]
    NoPosition(usize),
}

/// Errors that can occur when working with masses of atoms.
#[derive(Error, Debug, PartialEq, Eq)]
pub enum MassError {
    /// Used when the atom has no mass but it is required.
    #[error("{} atom with atom number '{}' has undefined mass", "error:".red().bold(), .0.to_string().yellow())]
    NoMass(usize),
}
