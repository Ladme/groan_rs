// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Simple, auxiliary functions and constants used through the `groan_rs` library.

use std::{f32::consts, path::Path};

use crate::structures::{simbox::SimBox, vector3d::Vector3D};

/******************************/
/*         CONSTANTS          */
/******************************/

/// PI times 2.
pub(crate) const PI_X2: f32 = consts::PI * 2.0f32;

/// Smallest coordinate supported by PDB. The actual minimal supported coordinate is
/// -99.9999 nm but due to floating point shenanigans, we are slightly more restrictive to be safe.
pub(crate) const PDB_MIN_COORDINATE: f32 = -99.0;
/// Largest coordinate supported by PDB. The actual maximal supported coordinate is
/// 999.9999 nm but due to floating point shenanigans, we are slightly more restrictive to be safe.
pub(crate) const PDB_MAX_COORDINATE: f32 = 999.0;

/// Smallest coordinate supported by GRO. The actual minimal supported coordinate is
/// -999.999 nm but due to floating point shenanigans, we are slightly more restrictive to be safe.
pub(crate) const GRO_MIN_COORDINATE: f32 = -999.0;
/// Largest coordinate supported by GRO. The actual maximal supported coordinate is
/// 9999.999 nm but due to floating point shenanigans, we are slightly more restrictive to be safe.
pub(crate) const GRO_MAX_COORDINATE: f32 = 9999.0;

/******************************/
/*      GROUP/LABEL NAMES     */
/******************************/

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

/******************************/
/*   CENTER OF MASS/GEOMETRY  */
/******************************/

/// Calculate contribution of an atom to the center of mass.
#[inline(always)]
pub(crate) fn center_atom_contribution(
    mut position: Vector3D,
    scaling: &Vector3D,
    simbox: &SimBox,
    mass: f32,
    sum_xi: &mut Vector3D,
    sum_zeta: &mut Vector3D,
) {
    // wrap position into the box
    position.wrap(simbox);

    let theta = Vector3D::new(
        position.x * scaling.x,
        position.y * scaling.y,
        position.z * scaling.z,
    );

    sum_xi.x += mass * theta.x.cos();
    sum_xi.y += mass * theta.y.cos();
    sum_xi.z += mass * theta.z.cos();

    sum_zeta.x += mass * theta.x.sin();
    sum_zeta.y += mass * theta.y.sin();
    sum_zeta.z += mass * theta.z.sin();
}

/// Convert coordinates from their representation on a circle to their representation on a line.
#[inline(always)]
pub(crate) fn from_circle_to_line(zeta: Vector3D, xi: Vector3D, scaling: &Vector3D) -> Vector3D {
    let theta = Vector3D::new(
        (-zeta.x).atan2(-xi.x) + consts::PI,
        (-zeta.y).atan2(-xi.y) + consts::PI,
        (-zeta.z).atan2(-xi.z) + consts::PI,
    );

    Vector3D::new(
        theta.x / scaling.x,
        theta.y / scaling.y,
        theta.z / scaling.z,
    )
}

/******************************/
/*           OTHER            */
/******************************/

/// Convert `impl AsRef<Path>` to `String` panicking with an error message in case the conversion fails.
#[inline(always)]
pub(crate) fn path2string(path: impl AsRef<Path>) -> String {
    path.as_ref()
        .to_str()
        .expect("FATAL GROAN ERROR | auxiliary::path2string | Could not convert Path to &str.")
        .to_owned()
}
