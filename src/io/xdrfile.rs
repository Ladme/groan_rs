// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Rust bindings for the `xdrfile` library.

use std::os::raw::{c_char, c_float, c_int};
use std::{
    ffi::{CString, NulError},
    path::Path,
};

use crate::errors::{ReadTrajError, TrajError};
use crate::prelude::TrajFile;
use crate::structures::simbox::SimBox;

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct CXdrFile {
    _empty: [u8; 0],
}

/// Module containing wrappers of C functions from the `xdrfile` library.
#[cfg(not(feature = "molly"))]
pub(crate) mod xtc {
    use super::*;

    extern "C" {
        pub fn read_xtc(
            xd: *mut CXdrFile,
            natoms: c_int,
            step: *mut c_int,
            time: *mut c_float,
            box_vec: *mut [[c_float; 3usize]; 3usize],
            x: *mut [c_float; 3usize],
            prec: *mut c_float,
        ) -> c_int;

        pub fn read_xtc_natoms(filename: *const c_char, natoms: *mut c_int) -> c_int;

        /// Jump to the frame which time is higher than or equal to `target_time`.
        pub fn xtc_jump_to_start(xd: *mut CXdrFile, target_time: c_float) -> c_int;

        /// Skip to the next frame in the xtc file without reading it.
        pub fn xtc_skip_frame(xd: *mut CXdrFile) -> c_int;

        /// Skip to the next frame in the xtc file without reading it. Read only information about the time of the frame.
        pub fn xtc_skip_frame_with_time(xd: *mut CXdrFile, time: *mut c_float) -> c_int;
    }
}

extern "C" {
    pub fn xdrfile_open(path: *const c_char, mode: *const c_char) -> *mut CXdrFile;

    pub fn xdrfile_close(xfp: *mut CXdrFile) -> c_int;

    pub fn write_xtc(
        xd: *mut CXdrFile,
        natoms: c_int,
        step: c_int,
        time: c_float,
        box_vec: *mut [[c_float; 3usize]; 3usize],
        x: *mut [c_float; 3usize],
        prec: c_float,
    ) -> c_int;

    /// Jump to the frame which time is higher than or equal to `target_time`
    pub fn trr_jump_to_start(xd: *mut CXdrFile, target_time: c_float) -> c_int;

    /// Skip to the next frame in the trr file without reading it.
    pub fn trr_skip_frame(xd: *mut CXdrFile) -> c_int;

    /// Skip to the next frame in the trr file without reading it. Read only information about the time of the frame.
    pub fn trr_skip_frame_with_time(xd: *mut CXdrFile, time: *mut c_float) -> c_int;

    pub fn read_trr(
        xd: *mut CXdrFile,
        natoms: c_int,
        step: *mut c_int,
        time: *mut c_float,
        lambda: *mut c_float,
        box_vec: *mut [[c_float; 3usize]; 3usize],
        x: *mut [c_float; 3usize],
        v: *mut [c_float; 3usize],
        f: *mut [c_float; 3usize],
    ) -> c_int;

    pub fn read_trr_natoms(filename: *const c_char, natoms: *mut c_int) -> c_int;

    pub fn write_trr(
        xd: *mut CXdrFile,
        natoms: c_int,
        step: c_int,
        time: c_float,
        lambda: c_float,
        box_vec: *mut [[c_float; 3usize]; 3usize],
        x: *mut [c_float; 3usize],
        v: *mut [c_float; 3usize],
        f: *mut [c_float; 3usize],
    ) -> c_int;
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum OpenMode {
    Read,
    Write,
}

#[derive(Debug)]
pub struct XdrFile {
    pub handle: *mut CXdrFile,
}

impl Drop for XdrFile {
    /// Close the file as XdrFile gets dropped.
    fn drop(&mut self) {
        unsafe {
            xdrfile_close(self.handle);
        }
    }
}

impl XdrFile {
    /// Open an xdr file returning a handle to the file.
    pub(super) fn open_xdr(filename: impl AsRef<Path>, mode: OpenMode) -> Result<Self, TrajError> {
        unsafe {
            let c_path = match path2cstring(filename.as_ref()) {
                Ok(x) => x,
                Err(_) => return Err(TrajError::InvalidPath(Box::from(filename.as_ref()))),
            };

            let handle = xdrfile_open(c_path.as_ptr(), mode2cstring(mode).as_ptr());

            if !handle.is_null() {
                Ok(XdrFile { handle })
            } else {
                Err(TrajError::FileNotFound(Box::from(filename.as_ref())))
            }
        }
    }
}

impl TrajFile for XdrFile {}

/// Convert Rust path to null-terminated C string.
pub fn path2cstring(path: impl AsRef<Path>) -> Result<CString, NulError> {
    CString::new(
        path.as_ref()
            .to_str()
            .expect(
                "FATAL GROAN ERROR | xdrfile::path2cstring | Could not convert path to CString.",
            )
            .as_bytes(),
    )
}

/// Convert XdrFile OpenMode to C string.
pub fn mode2cstring(mode: OpenMode) -> CString {
    match mode {
        OpenMode::Read => std::ffi::CStr::from_bytes_with_nul(b"r\0")
            .unwrap()
            .to_owned(),
        OpenMode::Write => std::ffi::CStr::from_bytes_with_nul(b"w\0")
            .unwrap()
            .to_owned(),
    }
}

/// Convert C box matrix from an xtc file to SimBox.
///
/// ## Returns
/// `SimBox` or `ReadTrajError::InvalidSimBox` if the simulation box is not valid.
pub fn matrix2simbox(matrix: [[c_float; 3usize]; 3usize]) -> Result<SimBox, ReadTrajError> {
    if matrix[0][1] != 0.0 || matrix[0][2] != 0.0 || matrix[1][2] != 0.0 {
        return Err(ReadTrajError::InvalidSimBox);
    }

    Ok([
        matrix[0][0],
        matrix[1][1],
        matrix[2][2],
        matrix[0][1],
        matrix[0][2],
        matrix[1][0],
        matrix[1][2],
        matrix[2][0],
        matrix[2][1],
    ]
    .into())
}

/// Convert SimBox to C box matrix for an xtc file.
/// If `simbox` is `None`, returns zeroed array.
pub fn simbox2matrix(simbox: Option<&SimBox>) -> [[c_float; 3usize]; 3usize] {
    match simbox {
        Some(simbox) => [
            [simbox.v1x, simbox.v1y, simbox.v1z],
            [simbox.v2x, simbox.v2y, simbox.v2z],
            [simbox.v3x, simbox.v3y, simbox.v3z],
        ],
        None => [[0.0; 3]; 3],
    }
}
