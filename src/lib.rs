// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! # groan_rs: Gromacs Analysis Library for Rust
//!
//! Rust library for analyzing Gromacs simulations.
//!
//! While the library is no longer in an _early_ stage of development, it is still unstable. Breaking changes can appear in any new version.
//!
//! ## What it can do
//! - Read and write [gro](`crate::io::gro_io::read_gro`), [pdb](`crate::io::pdb_io::read_pdb`), [pqr](`crate::io::pqr_io::read_pqr`), [ndx](`crate::system::System::read_ndx`), [xtc](`crate::system::System::xtc_iter`) and [trr](`crate::system::System::trr_iter`) files.
//! - Read topology and structure of the system from [tpr](`crate::io::tpr_io::read_tpr`) files (using the [`minitpr`](`minitpr`) crate).
//! - [Iterate over atoms](`crate::system::System::atoms_iter`) and access their [properties](`crate::structures::atom::Atom`), including connectivity (bonds).
//! - Select atoms using a [selection language](#groan-selection-language) similar to VMD.
//! - Calculate RMSD and perform RMSD fit for [individual](`crate::system::System::calc_rmsd`) [structures](`crate::system::System::calc_rmsd_and_fit`) and for [entire](`crate::prelude::RMSDTrajRead::calc_rmsd`) [trajectories](`crate::prelude::RMSDTrajRead::calc_rmsd_and_fit`).
//! - [Calculate distances between atoms](`crate::system::System::atoms_distance`) respecting periodic boundary conditions.
//! - [Select atoms based on geometric conditions.](`crate::system::System::group_create_from_geometry`)
//! - [Assign elements](`crate::system::System::guess_elements`) to atoms and [guess connectivity](`crate::system::System::guess_bonds`) between the atoms.
//! - Calculate [center of geometry](`crate::system::System::group_get_center`) or [mass](`crate::system::System::group_get_com`) for *any* group of atoms.
//! - [Center a group of atoms](`crate::system::System::atoms_center`) in a simulation box.
//! - Help with specific analyses by providing utility data structures (e.g., [`GridMaps`](`crate::structures::gridmap::GridMap`)).
//! - Some other minor stuff...
//!
//! ## What it CAN'T do (at the moment)
//! - Work with non-orthogonal periodic boundary conditions.
//! - Perform advanced analyses of structure and dynamics out of the box.
//!   (But `groan_rs` library tries to make it simple to implement your own!)
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
//! #### Analyzing structure files
//!
//! You can read structure files in various formats (GRO, PDB, PQR, TPR),
//! add groups from NDX files, and perform analyses.
//!
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error + Send + Sync>> {
//!     // read a structure file (in GRO, PDB, PQR, or TPR format)
//!     let mut system = System::from_file("system.gro")?;
//!
//!     // read an NDX file
//!     // NDX groups will be automatically created inside `system`
//!     system.read_ndx("index.ndx")?;
//!
//!     // calculate the center of geometry of a protein
//!     // 'Protein' is the name of a group from the NDX file
//!     let center = system.group_get_center("Protein")?;
//!
//!     // print the result
//!     println!("{:?}", center);
//!
//!     Ok(())
//! }
//! ```
//!
//! #### Selecting atoms
//!
//! You can select atoms using a query language similar to VMD (see 'groan selection language' below).
//!
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error + Send + Sync>> {
//!     let system = System::from_file("system.gro")?;
//!
//!     // select some atoms using the groan selection language and iterate through them
//!     for atom in system.selection_iter("serial 1-23 or (resname POPC and name P)")? {
//!         // perform some operation with the atom
//!     }
//!
//!     // you can temporarily store the selection iterator in a variable
//!     // see the `@protein` query? => groan can autodetect common structures
//!     // like atoms of membrane lipids or proteins
//!     let iterator = system.selection_iter("@protein")?;
//!
//!     Ok(())
//! }
//! ```
//!
//! #### Creating groups
//!
//! You can select atoms and save them into a group. This group is stored inside the system
//! and can be quickly accessed multiple times.
//!
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error + Send + Sync>> {
//!     let mut system = System::from_file("system.gro")?;
//!
//!     // select atoms and store them into a group 'Selected'
//!     system.group_create("Selected", "serial 1-23 or (resname POPC and name P)")?;
//!
//!     // iterate through the selected atoms
//!     for atom in system.group_iter("Selected")? {
//!         // perform some operation with the atom
//!     }
//!
//!     // use the previously created group to construct another group
//!     system.group_create("Some Atoms", "Selected || resid 87 to 124")?;
//!
//!     // modify the atoms in the group 'Some Atoms'
//!     for atom in system.group_iter_mut("Some Atoms")? {
//!         atom.set_residue_name("RES");
//!     }
//!
//!     // each system always contains two groups called 'all' and 'All'
//!     // these groups contain all the atoms in the system
//!     assert!(system.group_exists("all"));
//!     assert!(system.group_exists("All"));
//!
//!     Ok(())
//! }
//! ```
//!
//! #### Writing output structure files
//!
//! You can write structure files in GRO, PDB, or PQR format.
//!
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error + Send + Sync>> {
//!     let mut system = System::from_file("system.gro")?;
//!
//!     // write the entire system as a PDB file
//!     system.write_pdb("system.pdb", false)?;
//!
//!     // write only DNA atoms as a GRO file
//!     system.group_create("DNA", "@dna")?;
//!     system.group_write_gro("DNA", "dna.gro", false)?;
//!
//!     Ok(())
//! }
//! ```
//!
//! #### Geometry filtering of atoms
//!
//! You can select atoms based on geometric conditions.
//!
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error + Send + Sync>> {
//!     // when reading a TPR file, you get information about
//!     // the masses of atoms and connectivity between atoms
//!     let mut system = System::from_file("system.tpr")?;
//!
//!     // load groups into the system
//!     system.read_ndx("index.ndx")?;
//!
//!     // get the center of mass of the group 'Protein'
//!     let protein_com = system.group_get_com("Protein")?;
//!
//!     // construct a cylinder with its base in the protein center of mass,
//!     // a radius of 2 nm, height of 4 nm, and oriented along the z-axis
//!     let cylinder = Cylinder::new(protein_com, 2.0, 4.0, Dimension::Z);
//!
//!     // iterate over atoms of water molecules located inside the cylinder
//!     for atom in system.group_iter("Water")?.filter_geometry(cylinder) {
//!         // perform some operation with the atom
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! #### Analyzing trajectory files
//!
//! You can read trajectory files in XTC, TRR, or GRO format.
//!
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error + Send + Sync>> {
//!     let mut system = System::from_file("system.gro")?;
//!
//!     // load groups into the system
//!     system.read_ndx("index.ndx")?;
//!
//!     for frame in system.xtc_iter("trajectory.xtc")? {
//!         // check that the frame has been read correctly
//!         let frame = frame?;
//!
//!         // calculate distance between two groups of atoms
//!         let distance = frame.group_distance("Protein1", "Protein2", Dimension::XYZ)?;
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! You can also read only part of the trajectory file and/or skip some trajectory frames.
//! You can concatenate trajectory files and print the progress of any trajectory iteration.
//!
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error + Send + Sync>> {
//!     let mut system = System::from_file("system.gro")?;
//!
//!     // read multiple trajectory files with additional options
//!     for frame in system
//!         .xtc_cat_iter(&["md1.xtc", "md2.xtc", "md3.xtc"])?
//!         .with_range(50_000.0, 225_000.0)? // start at 50 ns, end at 225 ns
//!         .with_step(3)?                    // read every 3rd step
//!         .print_progress(ProgressPrinter::default()) {
//!
//!         // check that the frame has been read correctly
//!         let frame = frame?;
//!
//!         // continue working with the frame
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! #### Calculating RMSD
//!
//! You can calculate RMSD for two systems or for entire trajectories.
//!
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error + Send + Sync>> {
//!     let mut system = System::from_file("system.gro")?;
//!     let system2 = System::from_file("system2.gro")?;
//!
//!     // calculate RMSD between all atoms of `system` and `system2`
//!     let rmsd = system.calc_rmsd(&system2, "all")?;
//!
//!     // calculate RMSD for an entire trajectory
//!     let reference = system.clone();
//!     for frame_data in system.xtc_iter("trajectory.xtc")?.calc_rmsd(&reference, "all")? {
//!         let (frame, rmsd) = frame_data?;
//!
//!         println!("{}: {}", frame.get_simulation_time(), rmsd);
//!     }
//!
//!     Ok(())
//! }
//! ```
//!
//! #### Writing trajectory files
//!
//! You can not only read but also write trajectory files in XTC, TRR, or GRO format.
//!
//! ```no_run
//! use groan_rs::prelude::*;
//! use std::error::Error;
//!
//! fn main() -> Result<(), Box<dyn Error + Send + Sync>> {
//!     let mut system = System::from_file("system.gro")?;
//!
//!     // open an XTC file as output trajectory file
//!     system.xtc_writer_init("trajectory.xtc")?;
//!
//!     // read frames from a TRR trajectory and write them in an XTC format
//!     for frame in system.trr_iter("trajectory.trr")? {
//!         let frame = frame?;
//!
//!         // write the frame into all open trajectory files associated with the system
//!         frame.traj_write_frame()?;
//!     }
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
//! - Their **gro atom numbers** using `atomnum XYZ`. For instance, `atomnum 124` will select all atoms which have the atom number 124 in the `gro` or `pdb` file the system originates from. This can be multiple atoms.
//! - Their **chain identifiers** using `chain X`. For instance `chain A` will select all atoms belonging to the chain 'A'. The information about chains is usually present in pdb files. Note that if no chain information is present, using this keyword will select no atoms.
//! - Their **element names** using `element name XYZ` or `elname XYZ`. For instance, `element name carbon` will select all carbon atoms. Note that atoms must be assigned elements to be selected, read below.
//! - Their **element symbols** using `element symbol X` or `elsymbol X`. For instance, `element symbol C` will select all carbon atoms. Note that atoms must be assigned elements to be selected, read below.
//! - Their **labels** using `label XYZ`, see below.
//!
//! **Multiple identifiers**
//!
//! You can specify multiple identifiers in the query. For instance, by using `resname POPE POPG`, you will select all atoms of the system corresponding to residues named POPE as well as atoms corresponding to residues named POPG. See examples of similar queries below:
//!
//! - `resid 13 15 16 17` will select all atoms corresponding to residues with numbers 13, 15, 16, or 17.
//! - `name P CA HA` will select all atoms with atom names P, CA, or HA.
//! - `serial 245 267 269 271` will select atoms numbered 245, 267, 269, or 271.
//! - `chain A B C` will select atoms belonging to the chains 'A', 'B', or 'C'.
//! - `elname carbon hydrogen` will select all carbon and hydrogen atoms.
//! - `elsymbol C H` will select all carbon and hydrogen atoms.
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
//! You can select atoms using internally defined macros. Currently, `groan_rs` library provides six of such macros:
//! - `@protein` will select all atoms corresponding to amino acids (supports some 140 different amino acids).
//! - `@water` will select all atoms of water.
//! - `@ion` will select all atoms of ions.
//! - `@membrane` will select all atoms corresponding to membrane lipids (supports over 200 membrane lipid types).
//! - `@dna` will select all atoms belonging to a DNA molecule.
//! - `@rna` will select all atoms belonging to an RNA molecule.
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
//! **Selecting molecules**
//!
//! You can select all atoms that are part of a specific molecule using the `molecule with` (or `mol with`) operator. For instance, `molecule with serial 15` will select all atoms that are part of the same molecule as the atom with atom number 15 (including the atom 15).
//! Similarly, `molecule with resid 4 17 29` will select all atoms of all molecules which contain _any_ of the atoms of the residues 4, 17, or 29.
//! `molecule with name P` will then select all atoms of all molecules that contain an atom with name `P`. Note that a) by molecule we mean a collection of atoms connected by bonds and b) molecule is not a residue.
//!
//! The `molecule with` operator relates only to the first sub-query to its right. For instance, `molecule with serial 15 or name BB` will select all atoms that are either part of the same molecule as atom 15 or have name `BB`.
//! Meanwhile, `molecule with (serial 15 or name BB)` will select all molecules containing either the atom 15 or any atom with name `BB` (i.e., if a molecule contains an atom named BB, all its atoms will be selected).
//!
//! Note that to be able to select molecules, the system must contain topology information. Otherwise, no atoms are selected.
//!
//! **Labeling and selecting atoms**
//!
//! The `groan_rs` library allows you to label specific atoms by strings (see [`System::label_atom`](`crate::system::System::label_atom`)). Such labels can be used as identifiers in the Groan Selection Language. Each label is associated with a single atom which means that it can be used to select this atom. In other words, labels are similar to groups but are guaranteed to contain a single atom.
//!
//! For example, labeling a specific atom as `MyAtom` allows you to later select it using the query `label MyAtom`. You can also select multiple labeled atoms by chaining the labels together (e.g., `label MyAtom AnotherAtom OneMoreAtom`). Like groups, labels can be composed of multiple words. In such cases, the label must be enclosed by quotes (' or "), e.g., `label 'Very interesting atom'`.
//!
//! **Regular expressions**
//!
//! Atom, residue, and group names, as well as element names, element symbols, and labels can be specified using regular expressions. For instance all hydrogens atoms in the system can be selected using `name r'^[1-9]?H.*'` (or using `elname hydrogen`). Similarly, all phosphatidylcholine lipids can be selected using `resname r'^.*PC'`. Query `group r'^P'` or just `r'^P'` will then select atoms of all groups which name starts with 'P' (case sensitive).
//!
//! Note that the regular expression must be enclosed inside "regular expression block" starting with `r'` and ending with `'`. You can specify multiple regular expressions in a single query and use them alongside normal strings.
//!
//! Regular expressions are evaluated using the `regex crate`. For more information about the supported regular expressions, visit <https://docs.rs/regex/latest/regex/>.
//!
//! **Note on selecting elements**
//!
//! Note that the atoms of the system are not implicitly assigned elements by the `groan_rs` library. When using `element name` or `element symbol` keywords, the atoms will only be selected if they have been assigned an element.
//! Inside the `groan_rs` library, you can achieve this either by creating the `System` structure from a tpr file or by calling the `System::guess_elements` function.
//! If you are a user of a program employing the groan selection language, make sure that the program assigns elements to the atoms before relying on the `element` keyword.
//!
//! **Note on whitespace**
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
//! - **Serialization Support (`serde`)**: Enables the serialization and deserialization of `groan_rs` data structures through integration with the `serde` framework.
//! - **Concurrency Enhancement (`parallel`)**: Expands the `groan_rs` library with functions designed for multi-threaded execution.
//!
//! Install the `groan_rs` crate with a specific feature using `cargo add groan_rs --features [FEATURE]`.
//!
//! ## Limitations
//! - Currently, `groan_rs` library is not able to properly work with periodic simulation boxes that are **not orthogonal**.
//!   While it can read structures and trajectories with non-orthogonal boxes, calculated distances and similar properties may be incorrect!
//!   Tread very carefully!
//!
//! - While `groan_rs` can read double-precision trr and tpr files, it uses single-precision floating point numbers everywhere in its code.
//!   If you require double-precision for your analyses, look elsewhere.
//!
//! ## License
//! This library is released under the MIT License.

/// Current version of the `groan_rs` library.
pub const GROAN_VERSION: &str = env!("CARGO_PKG_VERSION");

mod auxiliary;
pub mod errors;
pub mod files;
pub mod io;
pub mod progress;
pub mod select;
pub mod structures;
pub mod system;
mod test_utilities;

/// Reexported basic `groan_rs` structures and traits.
pub mod prelude {
    pub use crate::files::FileType;
    pub use crate::io::gro_io::{GroReader, GroWriter};
    pub use crate::io::traj_read::{
        FrameData, FrameDataTime, TrajFile, TrajMasterRead, TrajRangeRead, TrajRangeReader,
        TrajRangeStepReader, TrajRead, TrajReadOpen, TrajReader, TrajStepRead, TrajStepReader,
        TrajStepTimeRead,
    };
    pub use crate::io::traj_write::TrajWrite;
    pub use crate::io::trr_io::{TrrReader, TrrWriter};
    pub use crate::io::xtc_io::{XtcReader, XtcWriter};
    pub use crate::progress::ProgressPrinter;
    pub use crate::structures::atom::Atom;
    pub use crate::structures::dimension::Dimension;
    pub use crate::structures::element::Elements;
    pub use crate::structures::gridmap::GridMap;
    pub use crate::structures::iterators::{
        AtomIterable, AtomIterator, AtomIteratorWithBox, AtomPairIterator, FilterAtomIterator,
        IntersectionAtomIterator, MoleculeIterator, MutAtomIterator, MutAtomIteratorWithBox,
        MutAtomPairIterator, MutFilterAtomIterator, MutMoleculeIterator, OrderedAtomIterator,
        OwnedAtomIterator, OwnedMutAtomIterator, UnionAtomIterator,
    };
    pub use crate::structures::shape::{Cylinder, Rectangular, Shape, Sphere, TriangularPrism};
    pub use crate::structures::simbox::SimBox;
    pub use crate::structures::vector3d::Vector3D;
    pub use crate::system::rmsd::RMSDTrajRead;
    #[cfg(any(feature = "parallel", doc))]
    pub use crate::system::ParallelTrajData;
    pub use crate::system::System;
}
