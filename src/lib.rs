// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

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
//! - Select atoms using a [selection language](https://ladme.github.io/gsl-guide/) similar to VMD.
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
//! You can select atoms using a [query language](https://ladme.github.io/gsl-guide/) similar to VMD.
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
//! # #[cfg(any(feature = "molly", not(feature = "no-xdrfile")))] {
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
//! # }
//! ```
//!
//! You can also read only part of the trajectory file and/or skip some trajectory frames.
//! You can concatenate trajectory files and print the progress of any trajectory iteration.
//!
//! ```no_run
//! # #[cfg(any(feature = "molly", not(feature = "no-xdrfile")))] {
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
//! # }
//! ```
//!
//! #### Calculating RMSD
//!
//! You can calculate RMSD for two systems or for entire trajectories.
//!
//! ```no_run
//! # #[cfg(any(feature = "molly", not(feature = "no-xdrfile")))] {
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
//! # }
//! ```
//!
//! #### Writing trajectory files
//!
//! You can not only read but also write trajectory files in XTC, TRR, or GRO format.
//!
//! ```no_run
//! # #[cfg(not(feature = "no-xdrfile"))] {
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
//! # }
//! ```
//!
//! ## Selecting atoms using GSL
//!
//! Groan selection language (GSL) is a query language for selecting atoms in the `groan_rs` library. In general, GSL is similar to the selection language used by VMD.
//! To understand the capabilities (and limitations) of the language, read [this guide](https://ladme.github.io/gsl-guide/).
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
//!
//! ### Default
//! - `molly`: **Blazingly Fast Reading of XTC Files**
//!   - Enables the use of the [`molly`](https://crates.io/crates/molly) crate for very fast reading of xtc files.
//!   - Enables the use of `GroupXtcReader` allowing partial read of XTC frames.
//!   - *This feature is enabled by default.*
//!   - If disabled, `xdrfile` library will be used instead to read XTC files and partial reading of XTC frames will not be supported.
//!
//! ### Additional
//! - `serde`: **Serialization Support**
//!   - Enables the serialization and deserialization of `groan_rs` data structures through integration with the `serde` framework.
//!
//! - `parallel`: **Concurrency**
//!   - Expands the `groan_rs` library with functions designed for multi-threaded execution.
//!
//! - `no-xdrfile`: **Pure Rust experience with no external C libraries**
//!   - **Removes** the compilation and use of `xdrfile` library.
//!   - If enabled, you will lose the ability to write XTC files and both read and write TRR files.
//!   - Only use this feature if the `xdrfile` library compilation fails.
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
    #[cfg(any(feature = "molly", not(feature = "no-xdrfile")))]
    pub use crate::io::xtc_io::XtcReader;

    #[cfg(not(feature = "no-xdrfile"))]
    pub use crate::io::trr_io::{TrrReader, TrrWriter};
    #[cfg(not(feature = "no-xdrfile"))]
    pub use crate::io::xtc_io::XtcWriter;

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
