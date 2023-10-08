
## Changelog for the `groan_rs` library

### Version 0.4.0
#### Changes to the Groan selection language
- **Breaking Change:** Group names can no longer include the characters '<', '>', or '='.
- Atom, residue, and group names can now be specified using regular expressions in the groan selection language.
- Operators (<, >, <=, and =>) have been introduced for open-ended ranges in the groan selection language.
- Added a new `@hydrogen` macro for automatic hydrogen atom detection.
- Implemented a new tokenizer for atom and residue numbers in the groan selection language.
#### Changes to reading xtc and trr files
- Introduced `XdrRangeReader` trait and `XtcRangeReader` and `TrrRangeReader` structures for efficient partial reading of xtc and trr files based on time ranges. `XdrRangeReader` and `TrrRangeReader` skip frames with times below the specified start time (atom properties from these frames are not read at all) and stop reading when the end time is reached.
#### Other changes
- Enhanced documentation for error variants within the `errors` module.

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


### Version 0.2.0
- Added functionality to write atom groups into xtc and trr files using the `XdrGroupWriter` trait along with the `XtcGroupWriter` and `TrrGroupWriter` structures.
- Introduced an optional property `chain` to the `Atom` structure, allowing specification of the chain ID to which the atom belongs. This property is automatically loaded from pdb files containing chain identifiers and it can also be set manually using the `Atom::with_chain` function.
- When writing pdb files, if the `chain` property is set for the `Atom`, it is now written into the pdb file.
- A new keyword `chain` has been added to the Groan Selection Language (GSL) that allows to select atoms based on the chain they belong to.
- Made the constructors for `AtomIterator` and `MutAtomIterator` public.
- Implemented operations for ndx-writable and ndx-nonwritable groups, including `System::group_make_writable`, `System::group_make_nonwritable`, and `System::group_names_writable`.
- Performed small code refactors.
- Made minor documentation fixes.
