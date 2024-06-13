// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Simple, auxiliary functions used through the `groan_rs` library.

/// Check whether the name for the group is a valid group name or a label name.
/// Characters '"&|!@()<>= are not allowed. Names containing whitespace only are also not allowed.
pub(crate) fn name_is_valid(string: &str) -> bool {
    if string.trim().is_empty() {
        return false;
    }

    let forbidden_chars = "'\"&|!@()<>=";

    for c in string.chars() {
        if forbidden_chars.contains(c) {
            return false;
        }
    }

    true
}
