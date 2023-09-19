
## Changelog for the `groan_rs` library

### Version 0.2.0
- Added functionality to write atom groups into xtc and trr files using the `XdrGroupWriter` trait along with the `XtcGroupWriter` and `TrrGroupWriter` structures.
- Introduced an optional property `chain` to the `Atom` structure, allowing specification of the chain ID to which the atom belongs. This property is automatically loaded from pdb files containing chain identifiers and it can also be set manually using the `Atom::with_chain` function.
- When writing pdb files, if the `chain` property is set for the `Atom`, it is now written into the pdb file.
- A new keyword `chain` has been added to the Groan Selection Language (GSL) that allows to select atoms based on the chain they belong to.
- Made the constructors for `AtomIterator` and `MutAtomIterator` public.
- Implemented operations for ndx-writable and ndx-nonwritable groups, including `System::group_make_writable`, `System::group_make_nonwritable`, and `System::group_names_writable`.
- Performed small code refactors.
- Made minor documentation fixes.

### Version 0.3.0
- Updated selection macros according to Gromacs definitions, and introduced new macros for identifying DNA (`@dna`) and RNA (`@rna`).
- Introduced the `System::group_create_ignore_warnings` function, which replicates the behavior of `System::group_create` while ignoring any WARNINGS that may occur during its execution. This new function addresses the issue of easily filtering out warnings returned by `System::group_create`, which uses a generic `Box<dyn Error>` return type. All other errors are propagated as usual since they indicate that the function failed.
- Renamed functions `System::group_by_resid` and `System::group_by_resname` to `System::atoms_split_by_resid` and `System::atoms_split_by_resname`, respectively. **This is a BREAKING CHANGE.**
- Introduced functions `System::group_split_by_resid` and `System::group_split_by_resname` for splitting groups of atoms by their residue ID and residue name, respectively.