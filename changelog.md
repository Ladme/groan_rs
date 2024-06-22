
## Changelog for the `groan_rs` library

### Version 0.8.0

#### Reading structures from TPR files
- In addition to system topology and simulation box, the positions, velocities, and forces of atoms and intermolecular bonds can be now also read from TPR files.

#### Reworked and expanded Atom iterators
- All custom iterators over atoms now implement one of two traits: `MasterAtomIterator` and `MasterMutAtomIterator`, depending on whether they provide references or *mutable* references to atoms.
- All mutable atom iterators can now call `MasterMutAtomIterator::translate` and `MasterMutAtomIterator::wrap` methods which translate and wrap, respectively, the atoms of the iterator.
- All immutable atom iterators can now call `MasterAtomIterator::get_center` and `MasterAtomIterator::get_com` methods to calculate center of geometry and center of mass, respectively, of the atoms of the iterator. This is useful for the calculation of local center of mass in the system.
- Iterators over atoms of the system can be now constructed directly from selection queries using `System::selection_iter` (for iteration over immutable atoms) and using `System::selection_iter_mut` (for iteration over mutable atoms). This allows selecting atoms without adding groups into the system.
- All custom iterators over atoms now implement `Debug` and `Clone`.

#### Labeled atoms
- Individual atoms can be now labeled with strings (`System::label_atom` and `System::select_and_label`).
- Labeled atoms can be accessed using their labels (`System::get_labeled_atom_as_ref`, `System::get_labeled_atom_as_mut`, `System::get_labeled_atom_copy`).
- You can iterate through labeled atoms using `System::labeled_atoms_iter` and `System::labeled_atoms_iter_mut`.
- Labeled atoms can be selected using the Groan Selection Language (`label 'atom X'`).

#### GridMap structure
- Implemented `GridMap` structure for easier analysis of properties of planar surfaces such as membranes.
- `GridMap` is a generic 2D array that can be navigated using coordinates of the molecular system.

#### Serde feature
- Serde support was added for `FileType`, `Dimension`, `Sphere`, `Rectangular`, `Cylinder`, and `TriangularPrism`.
- All (de)serializable `groan_rs` structures/enums now deny unknown fields.

#### Other changes
- **Breaking change:** `Group::name_is_valid` function has been moved to the `aux` module and is no longer public.
- **Breaking change:** Renamed several functions for consistency with the commonly used terminology:
  - `System::get_atoms_as_ref_mut` -> `System::get_atoms_as_mut`,
  - `System::get_groups_as_ref_mut` -> `System::get_groups_as_mut`,
  - `System::get_box_as_ref_mut` -> `System::get_box_as_mut`,
  - `System::get_atom_as_ref_mut` -> `System::get_atom_as_mut`.
- **Breaking Change:** `System::get_atoms_as_mut`, `System::get_groups_as_mut`, `System::group_remove`, and `System::group_rename` are no longer public (and also no longer marked as unsafe).Users should not be able to add/remove atoms or groups or change the properties and names of groups. Changing properties of atoms is obviously still allowed but that is better done using `System::atoms_iter_mut`.
- Implemented `System::traj_iter_map_reduce` for simple embarrassingly parallel iteration through simulation trajectories.
- Implemented `System::group_intersection` allowing to directly create groups that are intersections of other groups.
- Implemented `System::from_file_with_format` allowing to directly specify the format of the input file.
- When writing 'gro' and 'pdb' files, large coordinates can no longer 'overflow' the space allocated for them in the file format and instead an error is returned.
- Internal changes in the `ProgressPrinter` which potentially allow it to be used in multithreaded environment (untested!).
- Added unchecked (and unsafe!) variants of `System::get_atom_*` methods.

***

### Version 0.7.1
- Bug fix: GSL queries containing non-ASCII characters should no longer cause panic.
- `tpr_io::read_tpr` is now a public function.
- `TrajReadOpen` trait is now included in `groan_rs::prelude`.
- Replaced the `crossbeam_utils` crate with features from standard library.

### Version 0.7.0
#### Reading TPR files
- Implemented basic reading of TPR files using the `minitpr` crate. The `System::from_file` function can now recognize and parse TPR files.
- TPR files can also be explicitly read using the `tpr_io::read_tpr` function (import with `use groan_rs::io::tpr_io`).
- Currently, only the system name, simulation box dimensions, and the topology of the system (atoms and bonds) are being parsed.
- **Atom positions, velocities, and forces are not read from the TPR file by `groan_rs`.** All atoms in a system created by parsing a TPR file will have undefined positions and velocities.
- Reading intermolecular bonds and groups from TPR files is not currently supported.

#### Reading and writing PQR files
- The `System::from_file` function now supports reading whitespace-delimited PQR files.
- PQR files can also be explicitly read using the `pqr_io::read_pqr` function (import with `use groan_rs::io::pqr_io`).
- PQR files can be written with user-specified precision using the `System::write_pqr` and `System::group_write_pqr` functions.

#### `TrajConcatenator` structure
- All trajectory readers implementing certain traits (refer to `trajectory_readers.md` for more information) can now be concatenated using the `System::traj_cat_iter`.
- During concatenation, duplicate frames from trajectory boundaries are removed, allowing for simple iteration through simulation trajectories that directly follow each other.
- The `System::xtc_cat_iter` and `System::trr_cat_iter` functions are provided for direct concatenation of XTC and TRR files, respectively.
- Basic methods associated with trajectory readers, such as `with_range`, `with_step`, or `print_progress`, can be applied to the concatenation iterator.
- The `TrajConcatenator` structure ensures that the trajectory range is correctly identified even if it spans multiple trajectory files and that step iteration properly handles trajectory boundaries.

#### Replaced custom `Vector3D` implementation with faster `Vector3` from `nalgebra` crate
- **Breaking change:** When `Vector3D::to_unit` is applied to a vector with a length (magnitude) of zero, its fields are now set to `NaN` instead of the previously used `0.0`.

#### Restructuring of the module system
- **Breaking change:** The `Select` structure is now located directly in the `groan_rs::select` module, not in `groan_rs::selections::select`.
- **Breaking change:** The `System` structure is now located directly in the `groan_rs::system` module, not in `groan_rs::system::general`. This change does not affec users who utilize `groan_rs::prelude`.

#### `serde` feature
- Added a new feature, `serde`, which provides serialization and deserialization for several `groan_rs` structures, namely Vector3D, SimBox, Atom, AtomContainer, AtomBlock, Group, System, and Select.
- Select can be deserialized from a valid Groan Selection Language (GSL) query and serialized into a valid GSL query.
- This feature is not enabled by default and must be activated by using `cargo add groan_rs --features serde` or by adding `groan_rs = { version = "0.7.0", features = ["serde"] }` to your `Cargo.toml` file.
- Note that automated tests for serialization and deserialization currently only use the YAML file format, but other formats supported by the `serde` crate should also work.

#### `parallel` feature
- Introduced a new feature, `parallel`, which contains functions that can run on multiple threads.
- Activate with `cargo add groan_rs --features parallel` or by adding `groan_rs = { version = "0.7.0", features = ["parallel"] }`.
- Currently, only a parallelized version of `System::guess_bonds`, called `System::guess_bonds_parallel`, has been introduced.

#### Other changes
- **Breaking change:** The `parse_query` function from `crate::selections::select` has been transformed into an associated function of the `Select` structure and should now be called as `Select::parse_query`.
- Added a basic implementation of the `TriangularPrism` shape, which implements the `Shape` trait.
- Added a method, `System::group_create_from_select`, allowing groups to be created directly from `Select` structures.
- Clarified changes in the order of groups when removing or renaming a group from the System.

***

### Version 0.6.2
- Bug fix: `System::group_split_by_resid` and `System::group_split_by_resname` did not work correctly due to incorrect treatment of atom indices. This has now been fixed.

### Version 0.6.1
- Bug fix: When setting `ProgressPrinter` to write into a file, the printer did not start a new line once iteration was finished and instead printed new line into standard output. This has been now fixed.
- `ProgressPrinter` can be now set to NOT start a new line upon finishing the iteration through the trajectory (see `ProgressPrinter::with_newline_at_end`).

### Version 0.6.0
#### Simulation box is now an optional property of the system
- **Breaking change:** The simulation box is now an optional property of the system.
  - This change affects the return type of various functions that access the simulation box.
  - High-level methods of the `System` structure now return errors if they require a simulation box and it is not present.

#### Undefined atom positions no longer cause panic
- **Breaking change:** Higher-level `System` and `Atom` methods no longer panic when encountering an atom with an undefined position; instead, they return a recoverable error.

#### Assigning and guessing elements and bonds
- Added the `Elements` structure, which defines properties of supported elements.
- The `Elements` structure can be constructed from a YAML file containing element definitions by calling `Elements::from_file`.
- Default properties of the elements supported by `groan_rs` are stored in `src/config/elements.yaml`. A default `Elements` structure, containing elements and their properties from this YAML file, can be constructed using `Elements::default`.
- The `Elements` structure can be updated using `Elements::update`. In this case, another `Elements` structure must be provided, specifying the changes to be made to the original `Elements` structure.
- Elements can be guessed for the atoms of the system using `System::guess_elements`.
- Bonds can now be guessed for the system using `System::guess_bonds`. This requires van der Waals radii to be assigned to individual atoms. Van der Waals radii can be guessed based on the elements.

#### Changes to the `Atom` structure
- New optional fields added: `mass`, `element_name`, `element_symbol`, `vdw`, `expected_max_bonds`, and `expected_min_bonds`:
  - `mass` specifies the mass of the atom.
  - `element_name` specifies the name of the element assigned to this atom.
  - `element_symbol` specifies the symbol of the element.
  - `vdw` specifies the van der Waals radius used to guess bonds between atoms.
  - `expected_max_bonds` specifies the expected maximum number of bonds this atom may form.
  - `expected_min_bonds` specifies the expected minimum number of bonds this atom forms.
- Methods have been added to utilize and access these new fields.

#### Changes to the Groan Selection Language
- **Breaking change:** The macro `@hydrogen` is no longer supported. Use `element name hydrogen` instead.
- Bug fix: Invalid queries containing a command in parentheses followed by anything other than a binary operator (e.g., `(name CA CB) resname LYS`) no longer cause a panic but instead return a proper error.
- Atoms can now be selected based on their element name or element symbol (`element name carbon` or `element symbol C`). This requires atoms to have assigned elements.
- Atoms that are part of the same molecule can now be selected using the `molecule with` operator. For example, `molecule with serial 13` will select all atoms that are part of the same molecule as the atom with number 13. This requires the system to have topology information.

#### Center of mass calculations
- Introduced `System::group_get_com` for calculating the center of mass of a group of atoms using the Bai and Breen algorithm. Note that all atoms in the target group must have mass information.
- Introduced `System::atoms_center_mass` as an alternative to `System::atoms_center` but using center of mass instead of center of geometry.

#### Other changes
- **Breaking change:** The `System::set_mol_references` method is no longer public.
- Added more tests for the basic `Atom` `get_*`, `set_*`, `reset_*`, and `with_*` functions.
- Documentation examples no longer repeat that you have to use `groan_rs::prelude::*`.
- Added `Atom::distance_naive` for calculating the distance between atoms without considering periodic boundary conditions.
- Bug fix: `Vector3D::wrap_coordinate` and `Vector3D::min_image` now panic if `box_len` is 0 instead of looping indefinitely.
- The `System::group_get_center` function has been refactored and may return a result that is not binary identical to previous versions.
- Introduced `System::group_isempty` function for checking whether a target group of atoms is empty.

***

### Version 0.5.0
#### Important BREAKING CHANGES affecting the entire `groan_rs` library

- **Fields `position`, `velocity`, and `force` in the `Atom` structure are now of type `Option<Vector3D>`.**
  - `Atom::new` initializes these values to `None`. Use `Atom::with_position`, `Atom::with_velocity`, or `Atom::with_force` to set them.
  - Accessor methods (`Atom::get_position`, etc.) now return `Option<&Vector3D>`.
  - Distance functions require all atoms to have position information; absence leads to panic.
  - When writing files, undefined positions (and velocities) are output as 0.0 in all dimensions.
  - Reading xtc trajectories sets `velocity` and `force` of all atoms to `None`.
  - Atoms without `position` are excluded from geometric shape calculations in `(Mut)FilterAtomIterator` and `System::group_create_from_geometry`.
  - `System::has_forces`, `System::has_velocities`, `System::has_positions` return `true` only if *all* atoms have the respective property.
  - A position with x = 0.0, y = 0.0, z = 0.0 is now considered valid. However, atoms without positions are still written with 0.0 coordinates.

#### Other changes to the `Atom` structure

- New optional fields `charge` and `bonded` added:
  - `charge`: Specifies the atom's charge.
  - `bonded`: Contains indices of atoms bonded to the atom.
  - Methods added to use and access `charge` and `bonded`.
  - `Atom::wrap` method introduced for wrapping the atom into the simulation box.

#### `AtomContainer` structure

- **Breaking Change**: Introduction of `AtomContainer`:
  - Describes a collection of atoms in the `System` structure.
  - `Group` structure and all atom iterators reworked to use `AtomContainer`.
  - Many methods associated with `Group` removed/renamed/rewritten.

#### Reading and writing PDB files

- **Breaking Change**: `System::write_pdb` and `System::group_write_pdb` now require an argument for connectivity information.
- **Potentially Breaking Change**: Reading of PDB files now stops at `ENDMDL` or `END`. `END` keyword is written correctly at the end of files.
- Connectivity section can now be read with `System::add_bonds_from_pdb` and written with `System::write_pdb`/`System::group_write_pdb`.

#### Reading systems with non-orthogonal simulation boxes

- Support added for reading/writing systems with non-orthogonal simulation boxes from/to PDB, XTC, and TRR files.

#### New `System` connectivity methods

- New methods introduced:
  1. `System::add_bond`: Adds a bond between atoms.
  2. `System::bonded_atoms_iter` and `System::bonded_atoms_iter_mut`: Iterates through bonded atoms.
  3. `System::has_bonds`: Checks for available connectivity information.
  4. `System::molecule_iter` and `System::molecule_iter_mut`: Iterates through atoms in the same molecule.
  5. `System::make_molecules_whole`: Fixes molecules broken at periodic boundaries.

#### Quality-of-Life improvements

- `ProgressPrinter` introduced for trajectory reading progress. Activated via `TrajMasterRead::print_progress`.

#### Other changes

- **Breaking Change**: Atom access methods (`System::get_atom_as_ref`, etc.) now use index instead of `gmx_atom_number`.
- **Breaking Change**: `System::from_file` returns `Result<Self, Box<dyn Error + Send + Sync>>` for threading compatibility.
- `System::has_duplicate_atom_numbers`: Checks for duplicate atom numbers.
- `System::residues_renumber`: Renumbers residues in the system.
- `System::atoms_wrap` and `System::group_wrap`: Wraps atoms into the simulation box.
- `System::atoms_distance`: Calculates distance between two atoms.
- Panic groan errors now specify the originating function.

***

### Version 0.4.2
- Introduced `System::group_create_from_geometries` for group creation with multiple geometry constraints.
- Excluded test files from the `crates.io` crate distribution.
- Removed unnecessary fields from several error variants.
- Added new test for `System::group_create_from_geometry`.

### Version 0.4.1
- Added new `TrajMasterRead` trait that is implemented by `TrajReader`, `TrajRangeReader`, `TrajStepReader`, and `TrajRangeStepReader` allowing easier usage of these structures in generic functions.
- All new traits and structures for trajectory reading are now properly included in the `prelude`.

### Version 0.4.0
#### Changes to the Groan selection language
- **Breaking Change:** Group names can no longer include the characters '<', '>', or '='.
- Atom, residue, and group names can now be specified using regular expressions in the groan selection language.
- Operators (<, >, <=, and =>) have been introduced for open-ended ranges in the groan selection language.
- Added a new `@hydrogen` macro for automatic hydrogen atom detection.
- Implemented a new tokenizer for atom and residue numbers in the groan selection language.
#### Changes to reading xtc and trr files
- **Breaking Change:** Traits for reading xtc and trr files have been completely redone. Most notably, `Xdr`* methods, traits, and errors have been renamed to `Traj`*.
- All trajectory readers must now implement `TrajRead` trait and **do not** have to be iterable. Trajectory iterator is constructed by wrapping the `TrajRead`-implementing structure into `TrajReader`.
- Introduced several structures for efficient partial reading of trajectory files:
  - `TrajRangeReader` reads frames in a specified time range. Can be constructed for any structure implementing `TrajRangeRead` using `with_range` method.
  - `TrajStepReader` reads every `step`th frame of the trajectory. Can be constructed for any structure implementing `TrajStepRead` using `with_step` method.
  - `TrajRangeStepReader` reads every `step`th frame of the trajectory in a specified time range. Can be constructed for any structure implementing both `TrajRangeRead` and `TrajStepRead` using `with_range` and `with_step` methods.
#### Other changes
- Added `Atom::has_position` and `System::has_positions` methods that check whether the atom or all atoms in the system, respectively, have non-zero position(s).
- Added tests for reading double-precision trr files.
- Enhanced documentation for error variants within the `errors` module.

***

### Version 0.3.3
- `@ion` macro should no longer identify any part of proteins as ions.

### Version 0.3.2
- `@water` macro should no longer identify N-terminal (or any other) protein hydrogens as water hydrogens.

### Version 0.3.1
- Fixed failing test which occured due to the change in `GroupError::AlreadyExistsWarning` error message.

### Version 0.3.0
- **Breaking Change:** Renamed functions `System::group_by_resid` and `System::group_by_resname` to `System::atoms_split_by_resid` and `System::atoms_split_by_resname`, respectively.
- **Breaking Change:** Revised the return type of `System::group_create` from `Result<(), Box<dyn Error>>` to `Result<(), GroupError>`. This change simplifies pattern matching when handling the result of the function. Additionally, the error type `SelectError`, previously returned by the function, is now encapsulated within the new error variant `GroupError::InvalidQuery`.
- Updated selection macros according to Gromacs definitions, and introduced new macros for identifying DNA (`@dna`) and RNA (`@rna`).
- Introduced new functions: `System::group_split_by_resid` for splitting groups of atoms by their residue ID and `System::group_split_by_resname` for splitting by residue name.
- Added functions `System::group_remove` for removing groups and `System::group_rename` for renaming them.
- Introduced `System::group_create_from_geometry` which allows to construct groups of atoms that are inside some geometric shape in the current simulation frame.
- Restructured codebase:
  - Split the `system.rs` file into multiple smaller files.
  - Grouped files related to IO operations into a dedicated directory.
  - Organized files containing `System`-related methods into a separate directory.
  - Placed files defining fundamental data structures into their own directory.
- Made minor documentation tweaks.

***

### Version 0.2.0
- Added functionality to write atom groups into xtc and trr files using the `XdrGroupWriter` trait along with the `XtcGroupWriter` and `TrrGroupWriter` structures.
- Introduced an optional property `chain` to the `Atom` structure, allowing specification of the chain ID to which the atom belongs. This property is automatically loaded from pdb files containing chain identifiers and it can also be set manually using the `Atom::with_chain` function.
- When writing pdb files, if the `chain` property is set for the `Atom`, it is now written into the pdb file.
- A new keyword `chain` has been added to the Groan Selection Language (GSL) that allows to select atoms based on the chain they belong to.
- Made the constructors for `AtomIterator` and `MutAtomIterator` public.
- Implemented operations for ndx-writable and ndx-nonwritable groups, including `System::group_make_writable`, `System::group_make_nonwritable`, and `System::group_names_writable`.
- Performed small code refactors.
- Made minor documentation fixes.

***