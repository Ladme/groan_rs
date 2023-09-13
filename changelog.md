
## Changelog for the `groan_rs` library

### Version 0.2.0
- Writing atom groups into xtc and trr files using the `XdrGroupWriter` trait and `XtcGroupWriter` and `TrrGroupWriter` structures
- `AtomIterator` and `MutAtomIterator` constructors made public
- Operations for ndx-writable and ndx-nonwritable groups (`System::group_make_writable`, `System::group_make_nonwritable`, `System::group_names_writable`)
- Small refactors of code
- Small documentation fixes