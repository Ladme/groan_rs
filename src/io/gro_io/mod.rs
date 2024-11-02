// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of functions for reading and writing gro files.

pub mod structure;
pub mod trajectory;

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub use structure::read_gro;
pub use trajectory::GroReader;

use crate::errors::ParseGroError;
use crate::prelude::SimBox;

/// Read the next line in the provided buffer and parse it as a title.
fn get_title(
    buffer: &mut BufReader<File>,
    filename: impl AsRef<Path>,
) -> Result<String, ParseGroError> {
    let mut title = String::new();
    match buffer.read_line(&mut title) {
        Ok(0) | Err(_) => return Err(ParseGroError::LineNotFound(Box::from(filename.as_ref()))),
        Ok(_) => return Ok(title.trim().to_string()),
    };
}

/// Read the next line in the provided buffer and parse it as the number of atoms.
fn get_natoms(
    buffer: &mut BufReader<File>,
    filename: impl AsRef<Path>,
) -> Result<usize, ParseGroError> {
    let mut line = String::new();
    match buffer.read_line(&mut line) {
        Ok(0) | Err(_) => Err(ParseGroError::LineNotFound(Box::from(filename.as_ref()))),
        Ok(_) => match line.trim().parse::<usize>() {
            Ok(x) => Ok(x),
            Err(_) => Err(ParseGroError::ParseLineErr(line.trim().to_string())),
        },
    }
}

/// Parse a line as simulation box dimensions.
fn line_as_box(line: &str) -> Result<SimBox, ParseGroError> {
    let mut simulation_box = [0.0f32; 9];
    let mut i = 0usize;
    for split in line.split_whitespace() {
        simulation_box[i] = split
            .trim()
            .parse::<f32>()
            .map_err(|_| ParseGroError::ParseBoxLineErr(line.to_string()))?;
        i += 1;
    }

    if i != 3 && i != 9 {
        Err(ParseGroError::ParseBoxLineErr(line.to_string()))?;
    }

    // check that the simulation box is valid
    if simulation_box[3] != 0.0 || simulation_box[4] != 0.0 || simulation_box[6] != 0.0 {
        return Err(ParseGroError::UnsupportedBox(line.to_string()));
    }

    Ok(simulation_box.into())
}
