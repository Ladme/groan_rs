// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Rust bindings for the `xdrfile` library.

use std::os::raw::{c_char, c_float, c_int};
use std::{
    ffi::{CString, NulError},
    path::Path,
};

use crate::errors::TrajError;
use crate::structures::simbox::SimBox;

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct CXdrFile {
    _empty: [u8; 0],
}

extern "C" {
    pub fn xdrfile_open(path: *const c_char, mode: *const c_char) -> *mut CXdrFile;
}

extern "C" {
    pub fn xdrfile_close(xfp: *mut CXdrFile) -> c_int;
}

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
}

extern "C" {
    pub fn read_xtc_natoms(filename: *const c_char, natoms: *mut c_int) -> c_int;
}

extern "C" {
    pub fn write_xtc(
        xd: *mut CXdrFile,
        natoms: c_int,
        step: c_int,
        time: c_float,
        box_vec: *mut [[c_float; 3usize]; 3usize],
        x: *mut [c_float; 3usize],
        prec: c_float,
    ) -> c_int;
}

extern "C" {
    /// Jump to the frame which time is higher than or equal to `target_time`.
    pub fn xtc_jump_to_start(xd: *mut CXdrFile, target_time: c_float) -> c_int;
}

extern "C" {
    /// Jump to the frame which time is higher than or equal to `target_time`
    pub fn trr_jump_to_start(xd: *mut CXdrFile, target_time: c_float) -> c_int;
}

extern "C" {
    /// Skip to the next frame in the xtc file without reading it.
    pub fn xtc_skip_frame(xd: *mut CXdrFile) -> c_int;
}

extern "C" {
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
}

extern "C" {
    pub fn read_trr_natoms(filename: *const c_char, natoms: *mut c_int) -> c_int;
}

extern "C" {
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
    pub fn open_xdr(filename: impl AsRef<Path>, mode: OpenMode) -> Result<Self, TrajError> {
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

/// Convert Rust path to null-terminated C string.
pub fn path2cstring(path: impl AsRef<Path>) -> Result<CString, NulError> {
    CString::new(
        path.as_ref()
            .to_str()
            .expect("Groan error. Could not convert path to string.")
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
/// ## Warning
/// Currently only works with orthogonal simulation boxes.
pub fn matrix2simbox(matrix: [[c_float; 3usize]; 3usize]) -> SimBox {
    [matrix[0][0], matrix[1][1], matrix[2][2]].into()
}

/// Convert SimBox to C box matrix for an xtc file.
///
/// ## Warning
/// Currently only works with orthogonal simulation boxes.
pub fn simbox2matrix(simbox: &SimBox) -> [[c_float; 3usize]; 3usize] {
    [
        [simbox.x, 0.0, 0.0],
        [0.0, simbox.y, 0.0],
        [0.0, 0.0, simbox.z],
    ]
}
