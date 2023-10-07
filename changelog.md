
## Changelog for the `groan_rs` library

### Version 0.4.0
- **Breaking Change:** Characters '<', '>', and '=' are no longer allowed in group names.
- Atom, residue, and group names can be now specified using regular expressions.
- Introduced operators for open-ended ranges (<, >, <=, and =>) to the groan selection language.
- Added `@hydrogen` macro which can autodetect hydrogen atoms.
- Introduced new tokenizer for groan selection language atom and residue numbers.
- Added better documentation for the individual error variants in the `errors` module.

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
