// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! # groan_rs: Gromacs Analysis Library for Rust
//!
//! Rust library for analyzing Gromacs simulations.
//! Currently in very early stage of development:
//! anything can break, change or stop working at any time.
//!
//! ## Usage
//!
//! Run
//!
//! ```bash
//! $ cargo add groan_rs
//! ```
//!
//! Import the crate in your Rust code:
//! ```
//! use groan_rs::prelude::*;
//! ```
//!
//! ## Examples
//!
//! #### Analyzing a structure file
//!
//! Read a gro file and an ndx file and calculate the center of geometry of a protein.
//!
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     // read a gro file
//!     let mut system = System::from_file("system.gro")?;
//!     // `groan_rs` also supports pdb files which can be read as:
//!     // let mut system = System::from_file("system.pdb")?;
//!
//!     // read an ndx file
//!     system.read_ndx("index.ndx")?;
//!
//!     // calculate the center of geometry of a protein
//!     // note that the ndx file `index.ndx` must contain a group `Protein` for this to work
//!     let center = system.group_get_center("Protein")?;
//!
//!     // print the result
//!     println!("{:?}", center);
//!
//!     Ok(())
//! }
//!
//! ```
//!
//! #### Selecting atoms and creating groups
//!
//! Read a gro file and select atoms belonging to specific parts of the system.
//! `groan_rs` uses groan selection language (GSL) to select atoms.
//! GSL is similar to VMD query language, see below.
//!
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     // read a gro file
//!     let mut system = System::from_file("system.gro")?;
//!
//!     // select atoms using groan selection language creating a group 'My Group'
//!     system.group_create("My Group", "serial 25-28 or (resname POPC and name P)")?;
//!     // the group 'My Group' is now tightly associated with the system
//!
//!     // we can now use the previously created group to construct another group
//!     // note that
//!     // a) each atom can be in any number of groups
//!     // b) the group names are case-sensitive
//!     system.group_create("my group", "'My Group' || resid 87 to 124")?;
//!
//!     // we can then perform operations with the groups, e.g. write them into separate pdb files
//!     system.group_write_pdb("My Group", "My_Group.pdb")?;
//!     system.group_write_pdb("my group", "my_group.pdb")?;
//!
//!     // each system also be default contains two groups consisting of all atoms in the system
//!     // this groups are called 'All' and 'all'
//!     assert!(system.group_exists("all"));
//!     assert!(system.group_exists("All"));
//!
//!     // if you read an ndx file into the system like this:
//!     system.read_ndx("index.ndx")?;
//!     // the groups defined in the ndx file will be associated
//!     // with the system in the same way as the manually created groups
//!
//!     Ok(())
//! }
//! ```
//!
//! #### Analyzing a trajectory file
//!
//! Read an xtc file and calculate distance between two groups of atoms for each frame.
//! _(`groan_rs` supports procedural as well as functional approaches.)_
//!
//! Procedural approach:
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     // read a gro file
//!     let mut system = System::from_file("structure.gro")?;
//!
//!     // create the groups we are interested in
//!     // `groan_rs` uses VMD-like selection language for specifying groups of atoms
//!     system.group_create("group 1", "serial 1 to 5")?;
//!     system.group_create("group 2", "resid 45")?;
//!
//!     // prepare a vector for the calculated distances
//!     let mut distances = Vec::new();
//!
//!     // read the trajectory calculating the distance between the groups for each frame
//!     for frame in system.xtc_iter("files/md_short.xtc")? {
//!         // check that the xtc frame has been read correctly
//!         let frame = frame?;
//!         // calculate the distance and put it into the vector
//!         distances.push(
//!             frame
//!             .group_distance("group 1", "group 2", Dimension::XYZ)
//!             .expect("Groups do not exist but they should."),
//!         );
//!     }
//!
//!     // print the calculated distances
//!     println!("{:?}", distances);
//!
//!     Ok(())
//! }
//! ```
//!
//! Functional approach:
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     // read a gro file
//!     let mut system = System::from_file("structure.gro")?;
//!
//!     // create the groups we are interested in
//!     // `groan_rs` uses VMD-like selection language for specifying groups of atoms
//!     system.group_create("group 1", "serial 1 to 5")?;
//!     system.group_create("group 2", "resid 45")?;
//!
//!     // read the trajectory calculating distance between the groups for each frame
//!     // and collect the results into a vector
//!     let distances: Vec<f32> = system
//!         // read an xtc trajectory
//!         .xtc_iter("trajectory.xtc")?
//!         // calculate distance between the groups for each frame
//!         .map(|frame| {
//!             // check that the xtc frame has been read correctly
//!             let frame = frame?;
//!             // calculate the distance
//!             Ok(frame
//!                 .group_distance("group 1", "group 2", Dimension::XYZ)
//!                 .expect("Groups do not exist but they should."))
//!         })
//!         // collect the calculated distances
//!         // if any error occured while reading the trajectory, propagate it
//!         .collect::<Result<Vec<f32>, Box<dyn Error>>>()?;
//!
//!     // print the calculated distances
//!     println!("{:?}", distances);
//!
//!     Ok(())
//! }
//! ```
//!
//! #### Converting between trr and xtc files
//!
//! Read a trr file and write the velocities of particles into a new xtc file.
//!
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     // read a gro file
//!     let mut system = System::from_file("structure.gro")?;
//!
//!     // create an output xtc file
//!     // the `XtcWriter` structure is tightly coupled with the corresponding `System` structure
//!     // if the `System` is updated, the `XtcWriter` reflects this change
//!     let mut writer = XtcWriter::new(&system, "output.xtc")?;
//!
//!     // iterate through the trr trajectory
//!     // the `trr_iter` (as well as `xtc_iter`) function updates the `System` structure
//!     // with which it is associated with the properties of the system in the current frame
//!     for frame in system.trr_iter("input.trr")? {
//!         // check that the trr frame has been read correctly
//!         frame?;
//!         // write the current frame into the output xtc file
//!         writer.write_frame()?;
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! #### Iterating over atoms
//!
//! Iterate over atoms in two different groups, collecting information about them or modifying them.
//!
//! Procedural approach:
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     // read a gro file and an ndx file
//!     let mut system = System::from_file("system.gro")?;
//!     system.read_ndx("index.ndx")?;
//!
//!     // immutably iterate over atoms of the group 'Protein' collecting atom names
//!     let mut names = Vec::new();
//!     for atom in system.group_iter("Protein")? {
//!         names.push(atom.get_atom_name().to_owned());
//!     }
//!     // (you can also iterate over all the atoms in the system using `System::atoms_iter`)
//!
//!     // print the names
//!     println!("{:?}", names);
//!
//!     // mutably iterate over atoms of the group 'Membrane' changing their residue names
//!     for atom in system.group_iter_mut("Membrane")? {
//!         atom.set_residue_name("MEMB");
//!     }
//!     // (you can also iterate over all the atoms in the system using `System::atoms_iter_mut`)
//!
//!     Ok(())
//! }
//! ```
//!
//! Functional approach:
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     // read a gro file and an ndx file
//!     let mut system = System::from_file("system.gro")?;
//!     system.read_ndx("index.ndx")?;
//!
//!     // immutably iterate over atoms of the group 'Protein' collecting atom names
//!     let names: Vec<String> = system
//!         .group_iter("Protein")?
//!         .map(|atom| atom.get_atom_name().to_owned())
//!         .collect();
//!     // (you can also iterate over all the atoms in the system using `System::atoms_iter`)
//!
//!     // print the names
//!     println!("{:?}", names);
//!
//!     // mutably iterate over atoms of the group 'Membrane' changing their residue names
//!     system
//!         .group_iter_mut("Membrane")?
//!         .for_each(|atom| atom.set_residue_name("MEMB"));
//!     // (you can also iterate over all the atoms in the system using `System::atoms_iter_mut`)
//!
//!     Ok(())
//! }
//! ```
//!
//! #### Filtering and extracting atoms based on geometric criteria
//!
//! Filter atoms located inside a cylinder with specified position and dimensions and extract them into a separate vector.
//!
//! Procedural approach:
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     // read a gro file
//!     let system = System::from_file("system.gro")?;
//!
//!     // extract atoms that are located inside a specified cylinder
//!     let mut inside_cylinder = Vec::new();
//!     // the cylinder has its center of the base at x = 1.5, y = 2.5, z = 3.5 (in nm)
//!     // the cylinder has radius of 2.1 nm and height of 4.3 nm
//!     // the cylinder has its main axis oriented along the z-axis of the system
//!     let cylinder = Cylinder::new([1.5, 2.5, 3.5].into(), 2.1, 4.3, Dimension::Z);
//!
//!     for atom in system.atoms_iter() {
//!         if cylinder.inside(atom.get_position(), system.get_box_as_ref()) {
//!             inside_cylinder.push(atom.clone());
//!         }
//!     }
//!
//!     // atoms in the `inside_cylinder` vector are fully independent copies
//!     // of the original atoms in the `System` structure
//!
//!     // print the number of extracted atoms
//!     println!("{:?}", inside_cylinder.len());
//!
//!     Ok(())
//! }
//! ```
//!
//! Functional approach:
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error>> {
//!     // read a gro file
//!     let system = System::from_file("system.gro")?;
//!
//!     // extract atoms that are located inside a specified cylinder
//!     let inside_cylinder: Vec<Atom> = system
//!         .atoms_iter()
//!         // the cylinder has its center of the base at x = 1.5, y = 2.5, z = 3.5 (in nm)
//!         // the cylinder has radius of 2.1 nm and height of 4.3 nm
//!         // the cylinder has its main axis oriented along the z-axis of the system
//!         .filter_geometry(Cylinder::new(
//!             [1.5, 2.5, 3.5].into(),
//!             2.1,
//!             4.3,
//!             Dimension::Z,
//!         ))
//!         .cloned()
//!         .collect();
//!
//!     // atoms in the `inside_cylinder` vector are fully independent copies
//!     // of the original atoms in the `System` structure
//!
//!     // print the number of extracted atoms
//!     println!("{}", inside_cylinder.len());
//!
//!     Ok(())
//! }
//! ```
//!
//! ## Groan selection language
//!
//! Groan selection language (GSL) is a query language for selecting atoms in the `groan_rs` library. In general, GSL is similar to the selection language used by VMD.
//!
//! **Basic queries**
//!
//! You can select atoms based on:
//! - Their **residue names** using `resname XYZ`. For instance, `resname POPE` will select all atoms of the system corresponding to residues named POPE.
//! - Their **residue numbers** using `resid XYZ` or `resnum XYZ`. For instance, `resid 17` will select all atoms of the system corresponding to residues with number 17.
//! - Their **atom names** using `name XYZ` or `atomname XYZ`. For instance, `name P` will select all atoms of the system which name is P.
//! - Their **real atom numbers** using `serial XYZ`. For instance, `serial 256` will select the atom with the number 256 (this is guaranteed to be a single atom). Note that in this case, the atoms are selected based on their "real" atom numbers, as understood by gromacs, **not** by the atom numbers specified in the `gro` or `pdb` file the system originates from.
//! - Their **gro atom numbers** using `atomid XYZ` or `atomnum XYZ`. For instance, `atomid 124` will select all atoms which have the atom number 124 in the `gro` or `pdb` file the system originates from. This can be multiple atoms.
//! - Their **chain identifiers** using `chain X`. For instance `chain A` will select all atoms belonging to the chain 'A'. The information about chains is usually present in pdb files. Note that if no chain information is present, using this keyword will select no atoms.
//!
//! **Multiple identifiers**
//!
//! You can specify multiple identifiers in the query. For instance, by using `resname POPE POPG`, you will select all atoms of the system corresponding to residues named POPE as well as atoms corresponding to residues named POPG. See examples of similar queries below:
//!
//! - `resid 13 15 16 17` will select all atoms corresponding to residues with numbers 13, 15, 16, or 17.
//! - `name P CA HA` will select all atoms with atom names P, CA, or HA.
//! - `serial 245 267 269 271` will select atoms numbered 245, 267, 269, or 271.
//! - `chain A B C` will select atoms belonging to the chains 'A', 'B', or 'C'.
//!
//! **Selecting atoms using groups**
//!
//! You can also select atoms using previously created groups of atoms. For instance, if you have previously created groups named "Protein" and "Membrane", you can use query `group Protein Membrane` or just `Protein Membrane` to select atoms of these two groups. In case any of the groups does not exist, an error will be raised.
//!
//! If you load an `ndx file` for your system using `System::read_ndx()`, you can also use groups from the `ndx file`.
//!
//! In case your group consists of multiple words, you have to enclose it into quotes (' or "). For instance `Protein 'My Group'` will select all atoms of the group "Protein" as well as all atoms of the group "My Group".
//!
//! **Selecting atoms by autodetection**
//!
//! You can select atoms using internally defined macros. Currently, `groan_rs` library provides seven of such macros:
//! - `@protein` will select all atoms corresponding to amino acids (supports some 140 different amino acids).
//! - `@water` will select all atoms of water.
//! - `@ion` will select all atoms of ions.
//! - `@membrane` will select all atoms corresponding to membrane lipids (supports over 200 membrane lipid types).
//! - `@dna` will select all atoms belonging to a DNA molecule.
//! - `@rna` will select all atoms belonging to an RNA molecule.
//! - `@hydrogen` will select all hydrogen atoms.
//!
//! Please be aware that while groan macros are generally dependable, there is no absolute guarantee that they will unfailingly identify all atoms correctly. Be careful when using them.
//!
//! **Selecting all atoms**
//!
//! You can select all atoms of the system by using `all`.
//!
//! **Ranges**
//!
//! Instead of writing residue or atom numbers explicitly, you can use keyword `to` or `-` to specify a range. For example, instead of writing `resid 14 15 16 17 18 19 20`, you can use `resid 14 to 20` or `resid 14-20`. This will select all atoms corresponding to residues with residue numbers 14, 15, 16, 17, 18, 19, or 20.
//!
//! You can also specify multiple ranges in a single query or combine ranges with explicitly provided numbers. For example, `serial 1 3 to 6 10 12 - 14 17` will expand to `serial 1 3 4 5 6 10 12 13 14 17`.
//!
//! Open-ended ranges can be specified using the `<`, `>`, `<=`, and `>=` operators. For example, instead of writing `serial 1 to 180`, you can use `serial <= 180`. This will select all atoms with atom numbers lower than or equal to 180. Similarly, `resid > 33` will select atoms of all residues with residue number 34 or higher.
//!
//! Open-ended ranges can be combined with from-to ranges and explicitly provided numbers. For instance, `serial 1 3-6 >=20` will select all atoms with atom numbers 1, 3, 4, 5, 6, or 20 and higher.
//!
//! **Negations**
//!
//! Using keyword `not` or `!` in front of a query will negate it. For example, the query `not name CA` or `! name CA` will select all atoms which name does **not** correspond to CA. Similarly, `not resname POPE POPG` will select all atoms that correspond to residues with names other than POPE or POPG. `!Protein` will then select all atoms that are not part of the group named `Protein`.
//!
//! **Binary operations**
//!
//! You can combine basic queries by using `and` (`&&`) and `or` (`||`) operators.
//!
//! Joining two queries by `and` will select only atoms that were selected by **both** of the queries. For example, `resname POPE and name P` will select all atoms that belong to residues named POPE and that have the name P. Similarly, `resid 17 18 && serial 256 to 271` will select only atoms corresponding to residue 17 or 18 and with atom numbers between 256 and 271 (including 271).
//!
//! Joining two queries by `or` will select atoms that were selected by **either** of the queries (at least one of them). For example, `resname POPE or name P` will select all atoms that belong to residues named POPE as well as all atoms with the name P. Similarly, `resid 17 18 || serial 256 to 271` will select all atoms corresponding to residue 17 or 18 as well as all atoms with atom numbers between 256 and 271.
//!
//! In case multiple `and` and/or `or` operators are used in a single query, they are evaluated from left to right. For example, `resname POPE or name CA and not Protein` will select all atoms belonging to residues named POPE or having the atom name CA but all these atoms can not belong to the group named `Protein`.
//!
//! Autodetection macros can also be combined with other sub-queries using operators, i.e. `@membrane or group 'My Lipids'` will select all autodetected membrane lipids and all atoms of the group "My Lipids".
//!
//! **Parentheses**
//!
//! You can change the order in which the individual sub-queries and operations are evaluated by using parentheses `(` and `)`. Expressions enclosed in parentheses are evaluated first (think math). For example, `resname POPE or (name CA and not resid 18 to 21)` will select all atoms belonging to residues named POPE along with all atoms that
//! - have the atom name P **and**
//! - do not correspond to residues numbered 18 to 21.
//!
//! Meanwhile `(resname POPE or name CA) and not resid 18 to 21` is equivalent to `resname POPE or name CA and not resid 18 to 21`. This will select all atoms belonging to residues named POPE or having the atom name CA but all of these atoms can not belong to residues 18 to 21.
//!
//! You can place parenthetical expressions into other parenthetical expressions. For example `serial 1 to 6 or (name CA and resname POPE || (resid 1 to 7 or serial 123 to 128)) and Protein` is a valid query, albeit possibly way too convoluted.
//!
//! You can also place `not` (`!`) operator in front of a parenthetical expression. For example, `!(serial 1 to 6 && name P)` will select all atoms that do **not** have atom number between 1 and 6 while also having the atom name P.
//!
//! **Regular expressions**
//!
//! Atom, residue, and group names can be specified using regular expressions. For instance all hydrogens atoms in the system can be selected using `name r'^[1-9]?H.*'` (or using the `@hydrogen` macro which expands to the same regular expression). Similarly, all phosphatidylcholine lipids can be selected using `resname r'^.*PC'`. Query `group r'^P'` or just `r'^P'` will then select atoms of all groups which name starts with 'P' (case sensitive).
//!
//! Note that the regular expression must be enclosed inside "regular expression block" starting with `r'` and ending with `'`. You can specify multiple regular expressions in a single query and use them alongside normal strings.
//!
//! Regular expressions are evaluated using the `regex crate`. For more information about the supported regular expressions, visit <https://docs.rs/regex/latest/regex/>.
//!
//! **Note about whitespace**
//!
//! Operators and parentheses do not have to be separated by whitespace from the rest of the query, unless the meaning of the query would become unclear. For instance, `not(name CA)or(serial 1to45||Protein)` is a valid query, while `not(name CA)or(serial 1to45orProtein)` is **not** valid, as `orProtein` becomes uninterpretable. However, enclosing the `Protein` in parentheses, i.e. `not(name CA)or(serial 1to45or(Protein))`, turns the query into a valid one again.
//!
//! ## Error handling
//! Proper error handling and propagation is at heart of the `groan_rs` library.
//! The individual error types provided by the `groan_rs` are however not exported into the `prelude` module.
//!
//! If you want to use specific error type from the `groan_rs` library, you will have to include
//! it explicitly from the `errors` module. For instance, if you want to directly work with errors
//! that can occur when writing a pdb file, use:
//! ```
//! use groan_rs::errors::WritePdbError;
//! ```
//!
//! Note that `groan_rs` will still work correctly even if you do not explicitly include the error types.
//!
//! ## Features
//! - [x] reading and writing gro files
//! - [x] reading and writing xtc files
//! - [x] reading and writing ndx files
//! - [x] VMD-like selection language
//! - [x] basic geometry selection
//! - [x] center of geometry calculations
//! - [x] distance calculation respecting PBC
//! - [x] simulation frame centering
//! - [x] reading and writing pdb files
//! - [x] reading and writing trr files
//! - [ ] reading tpr files
//! - [ ] center of mass calculations
//! - [ ] support for non-orthogonal boxes
//!
//! ## Limitations
//! Currently, most of the `groan_rs` library only supports simulation boxes that are orthogonal!
//! If you intend to analyze simulations with non-orthogonal simulation boxes, look elsewhere.
//!
//! ## License
//! This library is released under the MIT License.

/// Current version of the `groan_rs` library.
pub const GROAN_VERSION: &str = env!("CARGO_PKG_VERSION");

pub mod errors;
pub mod files;
pub mod io {
    pub mod gro_io;
    mod ndx_io;
    pub mod pdb_io;
    pub mod trr_io;
    pub mod xdrfile;
    pub mod xtc_io;
}
pub mod iterators;
mod selections {
    mod name;
    mod numbers;
    pub mod select;
}
pub mod structures {
    pub mod atom;
    pub mod dimension;
    pub mod group;
    pub mod shape;
    pub mod simbox;
    pub mod vector3d;
}
pub mod system {
    mod analysis;
    pub mod general;
    mod groups;
    mod modifying;
    mod utility;
}

/// Reexported basic `groan_rs` structures and traits.
pub mod prelude {
    pub use crate::io::trr_io::{TrrGroupWriter, TrrRangeReader, TrrReader, TrrWriter};
    pub use crate::io::xdrfile::{XdrGroupWriter, XdrRangeReader, XdrReader, XdrWriter};
    pub use crate::io::xtc_io::{XtcGroupWriter, XtcRangeReader, XtcReader, XtcWriter};
    pub use crate::structures::atom::Atom;
    pub use crate::structures::dimension::Dimension;
    pub use crate::structures::shape::{Cylinder, Rectangular, Shape, Sphere};
    pub use crate::structures::simbox::SimBox;
    pub use crate::structures::vector3d::Vector3D;
    pub use crate::system::general::System;
}
