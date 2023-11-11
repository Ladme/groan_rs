
## Changelog for the `groan_rs` library

### Version 0.5.0
#### Changes to the Atom structure
- Two optional fields have been aded to the `Atom` structure: `charge` and `bonded`. `charge` specifies the charge of the atom. `bonded` is a list of atom indices of atoms that are bonded to the given atom.
- Several methods have been introduced that allow using and accessing the new `charge` and `bonded` properties.
#### AtomContainer structure
- **Breaking Change**: Introduced `AtomContainer`: a general structure describing a collection of atoms in the `System` structure. `Group` structure and all atom iterators have been reworked to employ `AtomContainer`. Many methods associated with `Group` have consequently been removed/renamed/rewritten.
#### Reading and writing PDB files
- **Breaking Change**: `System::write_pdb` and `System::group_write_pdb` now require additional argument specifying whether connectivity information should be written into the pdb file.
- **Potentially Breaking Change**: Reading of PDB files now properly ends once `ENDMDL` keyword is reached.
- The connectivity section of PDB files can be now read using `System::add_bonds_from_pdb` and written using `System::write_pdb`/`System::group_write_pdb`.
#### Quality-of-life improvements
- Introduced `ProgressPrinter` for printing the progress of trajectory reading. Progress printing can be turned on for any trajectory iteration by using the `TrajMasterRead::print_progress` method.
#### Other changes
- **Breaking Change**: `System::get_atom_as_ref`, `System::get_atom_as_ref_mut`, and `System::get_atom_copy` now take `index` of the atom in the `System` instead of `gmx_atom_number`. Indexing is now consistent with similar functions as it start from 0.
- Introduced `System::has_duplicate_atom_numbers` method which checks whether there are any atoms in the `System` structure sharing atom number.
- Reworked all panic groan errors to specify function from which they have been called.

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