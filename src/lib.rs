// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! # groan_rs: Gromacs Analysis Library Written in Rust
//!
//! Rust library for analyzing Gromacs simulations.
//! Currently in very early stage of development:
//! anything can break, change or stop working at any time.
//!
//! ## Usage
//!
//! Add the following dependency to your `Cargo.toml` file:
//!
//! ```toml
//! [dependencies]
//! groan_rs = "0.1.0"
//! ```
//!
//! Import the crate in your Rust code:
//! ```
//! use groan_rs::*;
//! ```
//!
//! ## Examples
//!
//! Read a gro file and an ndx file and calculate the center of geometry of a protein.
//! ```no_run
//! use groan_rs::*;
//!
//! fn main() {
//!     // read a gro file
//!     let mut system = System::from_file("system.gro").unwrap();
//!     // read an ndx file
//!     system.read_ndx("index.ndx").unwrap();
//!     // calculate center of geometry of a protein
//!     let center = system.group_get_center("Protein").unwrap();
//!
//!     // print the result
//!     println!("{:?}", center);
//! }
//!
//! ```
//!
//! Read an xtc file and calculate distance between two groups of atoms as a function of time.
//! ```no_run
//! use groan_rs::*;
//!
//! fn main() {
//!     // read a gro file
//!     let mut system = System::from_file("system.gro").unwrap();
//!
//!     // create the groups we are interested in
//!     // groan_rs uses VMD-like selection language for specifying groups of atoms
//!     system.group_create("group 1", "serial 1 to 5").unwrap();
//!     system.group_create("group 2", "resid 45").unwrap();
//!
//!     // read the trajectory calculating distance between the groups for each frame
//!     // and collect the results into a vector
//!     let distances: Vec<f32> = system.xtc_iter("trajectory.xtc").unwrap()
//!         .map(|frame| {
//!             frame.unwrap().group_distance("group 1", "group 2", Dimension::XYZ).unwrap()
//!         })
//!         .collect();
//!
//!     // print the calculated distances
//!     println!("{:?}", distances);
//! }
//! ```
//!
//! (Note that in a real program, you would want proper error handling instead of using unwrap.)
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
//! - [ ] reading and writing pdb files
//! - [ ] reading and writing trr files
//! - [ ] reading tpr files
//! - [ ] center of mass calculations
//! - [ ] support for non-orthogonal boxes
//!
//! ## Warning
//! Currently, most of the `groan_rs` library only supports simulation boxes that are orthogonal!
//! If you intend to analyze simulations with non-orthogonal simulation boxes, look elsewhere.
//!
//! ## License
//! This library is released under the MIT License.

mod analysis;
//pub mod analyzer;
pub mod atom;
mod c_xdrfile;
pub mod dimension;
pub mod errors;
mod gro_io;
mod group;
pub mod iterators;
mod ndx_io;
mod pdb_io;
mod select;
pub mod shape;
pub mod simbox;
pub mod system;
mod utility;
pub mod vector3d;
pub mod xtc_io;

pub use atom::Atom;
pub use dimension::Dimension;
pub use simbox::SimBox;
pub use system::System;
pub use vector3d::Vector3D;
pub use xtc_io::XtcWriter;
