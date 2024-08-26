// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of functions for reading and writing various structure, trajectory, and auxiliary files.

use crate::structures::atom::Atom;

pub mod gro_io;
pub mod ndx_io;
pub mod pdb_io;
pub mod pqr_io;
pub mod tpr_io;
pub mod traj_cat;
pub mod traj_io;
pub mod trr_io;
mod xdrfile;
pub mod xtc_io;

/// Check that no coordinate of any atom is smaller than `min_limit` or higher than `max_limit`.
/// This has to be done before writing the atoms into an output, otherwise a partial output could be writen.
/// Returns `true` if all coordinates are in the supported range (or undefined). Otherwise returns `false`.
fn check_coordinate_sizes<'a>(
    mut atoms: impl Iterator<Item = &'a Atom>,
    min_limit: f32,
    max_limit: f32,
) -> bool {
    atoms.all(|atom| {
        if let Some(pos) = atom.get_position() {
            [pos.x, pos.y, pos.z]
                .into_iter()
                .all(|coor| coor >= min_limit && coor <= max_limit)
        } else {
            true
        }
    })
}
