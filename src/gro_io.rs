// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of functions for reading and writing gro files.

use std::fs::File;
use std::io::{BufReader, BufRead};
use std::path::Path;

use crate::atom::Atom;
use crate::errors::ParseGroError;
use crate::system::System;
use crate::simbox::SimBox;

impl System {

    /// Create a new System from file. 
    /// ## Returns
    /// System structure if successful or ParseGroError if parsing fails. 
    /// ## Example
    /// ```no_run
    /// use groan_rs::system::System;
    /// 
    /// let system = match System::from_file("system.gro") {
    ///     Ok(x) => x,
    ///     Err(e) => {
    ///         eprintln!("{}", e);
    ///         return;
    ///     }
    /// };
    /// ```
    /// ## Notes
    /// - The returned System structure will contain two default groups "all" and "All"
    /// consisting of all the atoms in the system.
    pub fn from_file(filename: impl AsRef<Path>) -> Result<Self, ParseGroError> { 
        // TODO: add more supported file types
        read_gro(filename)
    }


    /*pub fn write_gro(&self, filename: impl AsRef<Path>) -> Result<(), WriteGroError> {

    }*/
}

/// Read a gro file and construct a System structure.
fn read_gro(filename: impl AsRef<Path>) -> Result<System, ParseGroError> {
    
    let file = match File::open(filename.as_ref()) {
        Ok(x) => x,
        Err(_) => return Err(ParseGroError::FileNotFound(Box::from(filename.as_ref()))),
    };

    let mut buffer = BufReader::new(file);

    // get title and number of atoms
    let title = get_title(&mut buffer, filename.as_ref())?;
    let n_atoms = get_natoms(&mut buffer, filename.as_ref())?;
    let mut simulation_box = [0.0f32; 9].into();

    let mut atoms: Vec<Atom> = Vec::with_capacity(n_atoms as usize);

    // parse all remaining lines
    for (gmx_index, raw_line) in buffer.lines().enumerate() {
        let line = match raw_line {
            Ok(x) => x,
            Err(_) => return Err(ParseGroError::LineNotFound(Box::from(filename.as_ref()))),
        };

        if gmx_index == n_atoms as usize { 
            simulation_box = line_as_box(&line)?;
        } else {
            let atom = line_as_atom(&line, gmx_index as u64 + 1)?;
            atoms.push(atom);
        }
    }

    if atoms.len() != n_atoms as usize {
        return Err(ParseGroError::LineNotFound(Box::from(filename.as_ref())));
    }

    Ok( System::new(&title, atoms, simulation_box) )

}

/// Read the next line in the provided buffer and parse it as a title.
fn get_title(buffer: &mut BufReader<File>, filename: impl AsRef<Path>) -> Result<String, ParseGroError> {

    let mut title = String::new();
    match buffer.read_line(&mut title) {
        Ok(0) | Err(_) => return Err(ParseGroError::LineNotFound(Box::from(filename.as_ref()))),
        Ok(_) => return Ok(title.trim().to_string()),
    };
}

/// Read the next line in the provided buffer and parse it as the number of atoms.
fn get_natoms(buffer: &mut BufReader<File>, filename: impl AsRef<Path>) -> Result<u64, ParseGroError> {

    let mut line = String::new();
    match buffer.read_line(&mut line) {
        Ok(0) | Err(_) => return Err(ParseGroError::LineNotFound(Box::from(filename.as_ref()))),
        Ok(_) => {
            match line.trim().parse::<u64>() {
                Ok(x) => return Ok(x),
                Err(_) => return Err(ParseGroError::ParseLineErr(line.trim().to_string())),
            };
        }
    };
}

/// Parse a line as atom.
fn line_as_atom(line: &str, gmx_index: u64) -> Result<Atom, ParseGroError> {

    if line.len() < 44 {
        return Err(ParseGroError::ParseAtomLineErr(line.to_string()));
    }

    // parse residue number
    let resid = line[0..5].trim().parse::<u32>().map_err(|_| ParseGroError::ParseAtomLineErr(line.to_string()))?;

    // parse residue name
    let resname = line[5..10].trim().to_string();
    if resname.is_empty() {
        return Err(ParseGroError::ParseAtomLineErr(line.to_string()));
    }

    // parse atom name
    let atomname = line[10..15].trim().to_string();
    if atomname.is_empty() {
        return Err(ParseGroError::ParseAtomLineErr(line.to_string()));
    }

    // parse atom number
    let atomid = line[15..20].trim().parse::<u32>().map_err(|_| ParseGroError::ParseAtomLineErr(line.to_string()))?;

    // parse position
    let mut position = [0.0; 3];
    for (i, item) in position.iter_mut().enumerate() {
        let curr = 20 + i * 8;
        *item = line[curr..curr + 8].trim().parse::<f32>().map_err(|_| ParseGroError::ParseAtomLineErr(line.to_string()))?;
    }

    // parse velocity, if present
    let mut velocity = [0.0; 3];
    if line.len() >= 68 {
        for (i, item) in velocity.iter_mut().enumerate() {
            let curr = 44 + i * 8;
            *item = line[curr..curr + 8].trim().parse::<f32>().map_err(|_| ParseGroError::ParseAtomLineErr(line.to_string()))?;
        }
    } 

    Ok(Atom::new(resid, 
                &resname, 
                atomid, 
                gmx_index, 
                &atomname, 
                position.into(), 
                velocity.into(),
                [0.0; 3].into()))
}

/// Parse a line as simulation box dimensions.
fn line_as_box(line: &str) -> Result<SimBox, ParseGroError> {

    let mut simulation_box = [0.0f32; 9];
    let mut i = 0usize;
    for split in line.split_whitespace() {
        simulation_box[i] = split.trim().parse::<f32>().map_err(|_| ParseGroError::ParseBoxLineErr(line.to_string()))?;
        i += 1;
    }

    if i != 3 && i != 9 {
        return Err(ParseGroError::ParseBoxLineErr(line.to_string()))?;
    }

    Ok(simulation_box.into())
}

// NOTE: Writing box dimensions: %10.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f


/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::approx_eq;

    #[test]
    fn test_read_gro() {

        let system = read_gro("test_files/example.gro").unwrap();

        // check the system's name
        assert_eq!(system.get_name(), "INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1");

        // check the number of atoms
        assert_eq!(system.get_n_atoms(), 16844);

        // check box size
        let simbox = system.get_box_as_ref();
        assert!(approx_eq!(f32, simbox.x, 13.01331));
        assert!(approx_eq!(f32, simbox.y, 13.01331));
        assert!(approx_eq!(f32, simbox.z, 11.25347));

        assert!(approx_eq!(f32, simbox.v1x, 13.01331));
        assert!(approx_eq!(f32, simbox.v2y, 13.01331));
        assert!(approx_eq!(f32, simbox.v3z, 11.25347));

        assert_eq!(simbox.v1y, 0.0f32);
        assert_eq!(simbox.v1z, 0.0f32);
        assert_eq!(simbox.v2x, 0.0f32);


        assert_eq!(simbox.v2z, 0.0f32);
        assert_eq!(simbox.v3x, 0.0f32);
        assert_eq!(simbox.v3y, 0.0f32);

        let atoms = system.get_atoms_as_ref();

        // check the first atom
        let first = &atoms[0];
        assert_eq!(first.get_residue_number(), 1);
        assert_eq!(first.get_residue_name(), "GLY");
        assert_eq!(first.get_atom_name(), "BB");
        assert_eq!(first.get_atom_number(), 1);
        assert_eq!(first.get_gmx_atom_number(), 1);
        
        assert!(approx_eq!(f32, first.get_position().x, 9.497));
        assert!(approx_eq!(f32, first.get_position().y, 1.989));
        assert!(approx_eq!(f32, first.get_position().z, 7.498));

        assert!(approx_eq!(f32, first.get_velocity().x, -0.0683));
        assert!(approx_eq!(f32, first.get_velocity().y, 0.1133));
        assert!(approx_eq!(f32, first.get_velocity().z, 0.0005));

        assert_eq!(first.get_force().x, 0.0f32);
        assert_eq!(first.get_force().y, 0.0f32);
        assert_eq!(first.get_force().z, 0.0f32);

        // check atom somewhere in the middle
        let middle = &atoms[4932];
        assert_eq!(middle.get_residue_number(), 435);
        assert_eq!(middle.get_residue_name(), "POPC");
        assert_eq!(middle.get_atom_name(), "C4B");
        assert_eq!(middle.get_atom_number(), 4933);
        assert_eq!(middle.get_gmx_atom_number(), 4933);
        
        assert!(approx_eq!(f32, middle.get_position().x, 6.384));
        assert!(approx_eq!(f32, middle.get_position().y, 11.908));
        assert!(approx_eq!(f32, middle.get_position().z, 5.471));

        assert!(approx_eq!(f32, middle.get_velocity().x, -0.2271));
        assert!(approx_eq!(f32, middle.get_velocity().y, 0.1287));
        assert!(approx_eq!(f32, middle.get_velocity().z, 0.1784));

        assert_eq!(middle.get_force().x, 0.0f32);
        assert_eq!(middle.get_force().y, 0.0f32);
        assert_eq!(middle.get_force().z, 0.0f32);

        // check the last atom
        let last = &atoms[16843];
        assert_eq!(last.get_residue_number(), 11180);
        assert_eq!(last.get_residue_name(), "ION");
        assert_eq!(last.get_atom_name(), "CL");
        assert_eq!(last.get_atom_number(), 16844);
        assert_eq!(last.get_gmx_atom_number(), 16844);
        
        assert!(approx_eq!(f32, last.get_position().x, 8.829));
        assert!(approx_eq!(f32, last.get_position().y, 11.186));
        assert!(approx_eq!(f32, last.get_position().z, 2.075));

        assert!(approx_eq!(f32, last.get_velocity().x, 0.0712));
        assert!(approx_eq!(f32, last.get_velocity().y, 0.2294));
        assert!(approx_eq!(f32, last.get_velocity().z, -0.1673));

        assert_eq!(last.get_force().x, 0.0f32);
        assert_eq!(last.get_force().y, 0.0f32);
        assert_eq!(last.get_force().z, 0.0f32);
    }

    #[test]
    fn test_read_gro_novelocities() {

        let system = read_gro("test_files/example_novelocities.gro").unwrap();

        assert_eq!(system.get_name(), "Buforin II peptide P11L");
        assert_eq!(system.get_n_atoms(), 50);

        // check box size
        let simbox = system.get_box_as_ref();
        assert!(approx_eq!(f32, simbox.x, 6.08608));
        assert!(approx_eq!(f32, simbox.y, 6.08608));
        assert!(approx_eq!(f32, simbox.z, 6.08608));

        assert_eq!(simbox.v1y, 0.0f32);
        assert_eq!(simbox.v1z, 0.0f32);
        assert_eq!(simbox.v2x, 0.0f32);

        assert_eq!(simbox.v2z, 0.0f32);
        assert_eq!(simbox.v3x, 0.0f32);
        assert_eq!(simbox.v3y, 0.0f32);


        let atoms = system.get_atoms_as_ref();

        // check the first atom
        let first = &atoms[0];
        assert_eq!(first.get_residue_number(), 1);
        assert_eq!(first.get_residue_name(), "THR");
        assert_eq!(first.get_atom_name(), "BB");
        assert_eq!(first.get_atom_number(), 1);
        assert_eq!(first.get_gmx_atom_number(), 1);
        
        assert!(approx_eq!(f32, first.get_position().x, 1.660));
        assert!(approx_eq!(f32, first.get_position().y, 2.061));
        assert!(approx_eq!(f32, first.get_position().z, 3.153));


        // check atom somewhere in the middle
        let middle = &atoms[24];
        assert_eq!(middle.get_residue_number(), 11);
        assert_eq!(middle.get_residue_name(), "LEU");
        assert_eq!(middle.get_atom_name(), "SC1");
        assert_eq!(middle.get_atom_number(), 25);
        assert_eq!(middle.get_gmx_atom_number(), 25);
        
        assert!(approx_eq!(f32, middle.get_position().x, 3.161));
        assert!(approx_eq!(f32, middle.get_position().y, 2.868));
        assert!(approx_eq!(f32, middle.get_position().z, 2.797));

        // check the last atom
        let last = &atoms[49];
        assert_eq!(last.get_residue_number(), 21);
        assert_eq!(last.get_residue_name(), "LYS");
        assert_eq!(last.get_atom_name(), "SC2");
        assert_eq!(last.get_atom_number(), 50);
        assert_eq!(last.get_gmx_atom_number(), 50);
        
        assert!(approx_eq!(f32, last.get_position().x, 4.706));
        assert!(approx_eq!(f32, last.get_position().y, 4.447));
        assert!(approx_eq!(f32, last.get_position().z, 2.813));

        // check that the velocity and force of all atoms is zero
        for atom in atoms.iter() {
            assert_eq!(atom.get_velocity().x, 0.0f32);
            assert_eq!(atom.get_velocity().y, 0.0f32);
            assert_eq!(atom.get_velocity().z, 0.0f32);

            assert_eq!(atom.get_force().x, 0.0f32);
            assert_eq!(atom.get_force().y, 0.0f32);
            assert_eq!(atom.get_force().z, 0.0f32);
        }
    }

    #[test]
    fn test_read_gro_box9() {

        let system = read_gro("test_files/example_box9.gro").unwrap();

        let simbox = system.get_box_as_ref();
        assert!(approx_eq!(f32, simbox.x, 6.08608));
        assert!(approx_eq!(f32, simbox.y, 6.08608));
        assert!(approx_eq!(f32, simbox.z, 6.08608));

        assert!(approx_eq!(f32, simbox.v1y, 1.4));
        assert!(approx_eq!(f32, simbox.v1z, 0.16));
        assert!(approx_eq!(f32, simbox.v2x, 2.2));

        assert!(approx_eq!(f32, simbox.v2z, 0.0));
        assert!(approx_eq!(f32, simbox.v3x, 1.4));
        assert!(approx_eq!(f32, simbox.v3y, 3.856));
    }

    #[test]
    fn test_read_gro_nonexistent() {

        if let Ok(_) = read_gro("test_files/nonexistent.gro") {
            panic!("Nonexistent file seems to exist.");
        }
    }

    macro_rules! read_gro_fails {
        ($name:ident, $file:expr, $variant:path, $expected:expr) => {
            #[test]
            fn $name() {
                let file = $file;
                match read_gro(file) {
                    Err($variant(e)) => assert_eq!(e, $expected),
                    Ok(_) => panic!("Parsing should have failed, but it succeeded."),
                    Err(e) => panic!("Parsing successfully failed but incorrect error type `{:?}` was returned.", e),
                }
            }
        };
    }

    read_gro_fails!(test_read_gro_incomplete_line, "test_files/example_incomplete_line.gro", 
    ParseGroError::ParseAtomLineErr, "   16HIS    SC1   35   3.458   3.653   ");

    read_gro_fails!(test_read_gro_empty_file, "test_files/example_empty.gro", 
    ParseGroError::LineNotFound, Box::from(Path::new("test_files/example_empty.gro")));

    read_gro_fails!(test_read_gro_only_title, "test_files/example_only_title.gro", 
    ParseGroError::LineNotFound, Box::from(Path::new("test_files/example_only_title.gro")));

    read_gro_fails!(test_read_gro_missing_natoms, "test_files/example_missing_natoms.gro", 
    ParseGroError::ParseLineErr, "");

    read_gro_fails!(test_read_gro_unparsable_natoms, "test_files/example_unparsable_natoms.gro", 
    ParseGroError::ParseLineErr, "6A5F");

    read_gro_fails!(test_read_gro_missing_atom, "test_files/example_missing_atom.gro", 
    ParseGroError::ParseAtomLineErr, "   6.08608   6.08608   6.08608");

    read_gro_fails!(test_read_gro_invalid_resid, "test_files/example_invalid_resid.gro", 
    ParseGroError::ParseAtomLineErr, "   1APHE    SC2   22   2.519   3.025   3.387");

    read_gro_fails!(test_read_gro_invalid_atomid, "test_files/example_invalid_atomid.gro", 
    ParseGroError::ParseAtomLineErr, "   21LYS     BB        4.362   4.008   3.161");

    read_gro_fails!(test_read_gro_invalid_position, "test_files/example_invalid_position.gro", 
    ParseGroError::ParseAtomLineErr, "    2ARG    SC1    4   1.877   1. 73   3.023");

    read_gro_fails!(test_read_gro_invalid_velocity, "test_files/example_invalid_velocity.gro", 
    ParseGroError::ParseAtomLineErr, "   15LEU    SC1   31   9.638   2.052   5.595  0.0685  O.0634  0.1453");

    read_gro_fails!(test_read_gro_shifted_line, "test_files/example_shifted_line.gro", 
    ParseGroError::ParseAtomLineErr, "    20ARG     BB   45   4.265   3.832   2.925");

    read_gro_fails!(test_read_gro_empty_box_line, "test_files/example_empty_box_line.gro", 
    ParseGroError::ParseBoxLineErr, "");

    read_gro_fails!(test_read_gro_short_box, "test_files/example_short_box.gro",
    ParseGroError::ParseBoxLineErr, "   6.08608   6.08608");

    read_gro_fails!(test_read_gro_long_box, "test_files/example_long_box.gro",
    ParseGroError::ParseBoxLineErr, "   6.08608   6.08608   6.08608   6.08608   6.08608");

    read_gro_fails!(test_read_gro_unparsable_box, "test_files/example_unparsable_box.gro",
    ParseGroError::ParseBoxLineErr, "   6.08608   6.08608   6,08608");

    #[test]
    fn test_system_from_file() {

        let system = System::from_file("test_files/example.gro").expect("File not found.");

        assert_eq!(system.get_name(), "INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1");
        assert_eq!(system.get_n_atoms(), 16844);

        assert!(system.group_exists("all"));
    }

    #[test]
    fn test_system_from_file_fails() {
        
        match System::from_file("test_files/example_invalid_position.gro") {
            Ok(_) => panic!("Parsing should have failed, but it succeeded."),
            Err(_) => (),
        }
    }

}