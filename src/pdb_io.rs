// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of functions for reading and writing pdb files.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::atom::Atom;
use crate::errors::ParsePdbError;
use crate::simbox::SimBox;
use crate::system::System;
use crate::vector3d::Vector3D;

/// Read a pdb file and construct a System structure.
///
/// ## Warning
/// Currently only supports orthogonal simulation boxes!
///
/// ## Supported keywords
/// This function can handle lines starting with ATOM, HETATM, TITLE and CRYST1.
/// All other lines are ignored.
///
/// ## Notes
/// - In case multiple TITLE lines are provided, the **last one** is used as the
/// name of the system. If not TITLE line is provided, "Unknown" is used as the name.
///
/// - In case multiple CRYST1 lines are provided, information from the **last one** is used.
/// If no CRYST1 line is provided, the simulation box size is set to 0 in all dimensions.
/// Note again that only orthogonal simulation boxes are supported.
///
pub fn read_pdb(filename: impl AsRef<Path>) -> Result<System, ParsePdbError> {
    let file = match File::open(filename.as_ref()) {
        Ok(x) => x,
        Err(_) => return Err(ParsePdbError::FileNotFound(Box::from(filename.as_ref()))),
    };

    let reader = BufReader::new(file);

    let mut atoms: Vec<Atom> = Vec::new();
    let mut title = "Unknown".to_string();
    let mut simbox = SimBox::from([0.0, 0.0, 0.0]);

    for raw_line in reader.lines() {
        let line = match raw_line {
            Ok(x) => x,
            Err(_) => return Err(ParsePdbError::LineNotFound(Box::from(filename.as_ref()))),
        };

        // parse ATOM/HETATM line
        if (line.len() >= 4 && line[0..4] == "ATOM".to_string())
            || (line.len() >= 6 && line[0..6] == "HETATM".to_string())
        {
            atoms.push(line_as_atom(&line)?);
        }
        // parse TITLE line
        else if line.len() >= 5 && line[0..5] == "TITLE".to_string() {
            title = line_as_title(&line)?;
        }
        // parse CRYST1 line
        else if line.len() >= 6 && line[0..6] == "CRYST1".to_string() {
            simbox = line_as_box(&line)?;
        }
    }

    Ok(System::new(&title, atoms, simbox))
}

/// Parse a single line from a pdb file as an atom.
///
/// ## Notes
/// - Parses lines starting with ATOM or HETATM.
fn line_as_atom(line: &str) -> Result<Atom, ParsePdbError> {
    // check line length
    if line.len() < 54 {
        return Err(ParsePdbError::ParseAtomLineErr(line.to_string()));
    }

    // parsing atom number
    let atom_number = line[6..11]
        .trim()
        .parse::<usize>()
        .map_err(|_| ParsePdbError::ParseAtomLineErr(line.to_string()))?;

    // parsing atom name
    let atom_name = line[12..16].trim().to_string();
    if atom_name.is_empty() {
        return Err(ParsePdbError::ParseAtomLineErr(line.to_string()));
    }

    // parsing residue name
    let residue_name = line[17..21].trim().to_string();
    if residue_name.is_empty() {
        return Err(ParsePdbError::ParseAtomLineErr(line.to_string()));
    }

    // parsing residue number
    let residue_number = line[22..26]
        .trim()
        .parse::<usize>()
        .map_err(|_| ParsePdbError::ParseAtomLineErr(line.to_string()))?;

    // parsing position
    let mut curr = 30usize;
    let mut position = [0.0, 0.0, 0.0];
    for i in 0..3 {
        position[i] = line[curr..curr + 8]
            .trim()
            .parse::<f32>()
            .map(|x| x / 10.0)
            .map_err(|_| ParsePdbError::ParseAtomLineErr(line.to_string()))?;

        curr += 8;
    }

    Ok(Atom::new(
        residue_number,
        &residue_name,
        atom_number,
        &atom_name,
        position.into(),
        Vector3D::default(),
        Vector3D::default(),
    ))
}

/// Parse a single line as a simulation box.
///
/// ## Notes
/// - Parses a line starting with CRYST1.
fn line_as_box(line: &str) -> Result<SimBox, ParsePdbError> {
    // check line length
    if line.len() < 33 {
        return Err(ParsePdbError::ParseBoxLineErr(line.to_string()));
    }

    let mut boxsize = [0.0, 0.0, 0.0];
    let mut curr = 6usize;
    for i in 0..3 {
        boxsize[i] = line[curr..curr + 9]
            .trim()
            .parse::<f32>()
            .map(|x| x / 10.0)
            .map_err(|_| ParsePdbError::ParseBoxLineErr(line.to_string()))?;

        curr += 9;
    }

    Ok(SimBox::from(boxsize))
}

/// Parse a single line as a title.
///
/// ## Notes
/// - Parses a line starting with TITLE.
fn line_as_title(line: &str) -> Result<String, ParsePdbError> {
    let title = line[5..].trim().to_string();
    if title.is_empty() {
        return Err(ParsePdbError::ParseTitleLineErr(line.to_string()));
    }

    Ok(title)
}
