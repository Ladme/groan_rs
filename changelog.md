
## Changelog for the `groan_rs` library

### Version 0.6.0

#### Simulation box is now optional property of the system
- **Breaking change:** Simulation box is now an optional property of the system.
  - This causes changes in the return type of various functions accessing the simulation box.
  - High-level methods of the `System` structure return errors when they require a simulation box but it is not present.

#### Undefined atom positions no longer cause panic
- **Breaking change:** Higher-level `System` and `Atom` methods now do not panic when encountering an atom with undefined position and instead return a recoverable error.

#### Assigning and guessing elements and bonds
- Added `Elements` structure which defines properties of the supported elements.
- `Elements` structure can be constructed from a YAML file containing the element definitions by calling `Elements::from_file`.
- Default properties of the elements supported by `groan_rs` are stored in `src/config/elements.yaml`. Default `Elements` structure, containing elements and their properties from this YAML file, can be constructed using `Elements::default`.
- `Elements` structure can be updated using `Elements::update`. In such case, another `Elements` structure must be provided which specifies changes to be made in the original `Elements` structure.
- Bonds can be now guessed for the system using `System::guess_bonds`. This requires van der Waals radii to be assigned to the individual atoms. Van der Waals radii can be guessed from the elements.

#### Changes to the `Atom` structure
- New optional fields `mass`, `element_name`, `element_symbol`, `vdw`, and `expected_max_bonds` added:
  - `mass` specifies mass of the atom.
  - `element_name` specifies name of the element assigned to this atom.
  - `element_symbol` specifies the symbol of the element.
  - `vdw` specifies the van der Waals radius used to guess bonds between atoms.
  - `expected_max_bonds` specifies the expected maximal number of bonds this atom may form.
  - `expected_min_bonds` specifies the expected minimal number of bonds this atom may form.
- Methods were added to use and access the new fields.

#### Groan Selection Language
- **Breaking change:** Macro `@hydrogen` is no longer supported. Instead, use `element name hydrogen`.
- Bug fix: Invalid queries containing a command in parentheses followed by anything other than a binary operator (e.g. `(name CA CB) resname LYS`) no longer cause a panic, but return a proper error.
- Atoms can be now selected based on their element name or element symbol (`element name carbon` or `element symbol C`). This requires atoms to be assigned elements.
- Atoms that are part of the same molecule can be now selected using the `molecule with` operator. For instance, `molecule with serial 13` will select all atoms that are part of the same molecule as atom with number 13. This requires the system to have topology information.

#### Other changes
- **Breaking change:** `System::set_mol_references` method is no longer public.
- Added more tests for basic `Atom` `get_*`, `set_*`, `reset_*`, and `with_*` functions.
- Documentation examples no longer repeat that you have to use `groan_rs::prelude::*`.
- Added `Atom::distance_naive` for calculating distance between atoms without taking periodic boundary conditions into consideration.
- Bug fix: `Vector3D::wrap_coordinate` and `Vector3D::min_image` now panic if `box_len` is 0 instead of looping indefinitely.

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