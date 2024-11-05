// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of functions for reading and writing gro files.

pub mod structure;
pub mod trajectory;

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

pub use structure::read_gro;
pub use trajectory::GroReader;
pub use trajectory::GroWriter;

use crate::errors::{ParseGroError, WriteGroError};
use crate::prelude::{AtomIterator, SimBox};
use crate::system::System;

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

/// Write box dimensions into an open gro file.
fn write_box(simbox: Option<&SimBox>, writer: &mut BufWriter<File>) -> Result<(), WriteGroError> {
    match simbox {
        Some(simbox) if simbox.is_orthogonal() => {
            writeln!(
                writer,
                " {:9.5} {:9.5} {:9.5}",
                simbox.x, simbox.y, simbox.z
            )
            .map_err(|_| WriteGroError::CouldNotWrite)?;
        }
        Some(simbox) => {
            writeln!(
                writer,
                " {:9.5} {:9.5} {:9.5} {:9.5} {:9.5} {:9.5} {:9.5} {:9.5} {:9.5}",
                simbox.x,
                simbox.y,
                simbox.z,
                simbox.v1y,
                simbox.v1z,
                simbox.v2x,
                simbox.v2z,
                simbox.v3x,
                simbox.v3y
            )
            .map_err(|_| WriteGroError::CouldNotWrite)?;
        }
        None => {
            let x = 0.0;
            writeln!(writer, " {x:9.5} {x:9.5} {x:9.5}",)
                .map_err(|_| WriteGroError::CouldNotWrite)?;
        }
    }
    Ok(())
}

/// Write gro file header into an open gro file.
#[inline(always)]
fn write_header(
    writer: &mut BufWriter<File>,
    title: &str,
    n_atoms: usize,
) -> Result<(), WriteGroError> {
    writeln!(writer, "{}", title).map_err(|_| WriteGroError::CouldNotWrite)?;
    writeln!(writer, "{:>5}", n_atoms).map_err(|_| WriteGroError::CouldNotWrite)?;

    Ok(())
}

/// Determine the title to use for the gro frame.
#[inline(always)]
fn determine_title(system: &System, group: &str, is_trajectory: bool) -> String {
    let title = match group {
        "all" => system.get_name().to_owned(),
        _ => format!("Group `{}` from {}", group, system.get_name()),
    };

    if is_trajectory {
        format!(
            "{} t={} step={}",
            title,
            system.get_simulation_time(),
            system.get_simulation_step()
        )
    } else {
        title
    }
}

/// Write a single simulation frame in gro format.
/// Does not check coordinate sizes.
fn write_frame(
    system: &System,
    writer: &mut BufWriter<File>,
    group_name: &str,
    iterator: AtomIterator,
    n_atoms: usize,
    write_velocities: bool,
    is_trajectory: bool,
) -> Result<(), WriteGroError> {
    // write gro file header
    write_header(
        writer,
        &determine_title(system, group_name, is_trajectory),
        n_atoms,
    )?;

    // write atoms
    for atom in iterator {
        atom.write_gro(writer, write_velocities)?;
    }

    // write simulation box
    write_box(system.get_box(), writer)?;

    writer.flush().map_err(|_| WriteGroError::CouldNotWrite)?;

    Ok(())
}
