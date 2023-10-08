// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Rust bindings for the `xdrfile` library.

/******************************/
/*   C bindings for XdrFile   */
/******************************/

use std::os::raw::{c_char, c_float, c_int};
use std::{
    ffi::{CString, NulError},
    path::Path,
};

use crate::errors::{ReadXdrError, WriteXdrError, XdrError};
use crate::structures::simbox::SimBox;
use crate::system::general::System;

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

/******************************/
/*        Rust XdrFile        */
/******************************/

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
    pub fn open_xdr(filename: impl AsRef<Path>, mode: OpenMode) -> Result<Self, XdrError> {
        unsafe {
            let c_path = match path2cstring(filename.as_ref()) {
                Ok(x) => x,
                Err(_) => return Err(XdrError::InvalidPath(Box::from(filename.as_ref()))),
            };

            let handle = xdrfile_open(c_path.as_ptr(), mode2cstring(mode).as_ptr());

            if !handle.is_null() {
                Ok(XdrFile { handle })
            } else {
                Err(XdrError::FileNotFound(Box::from(filename.as_ref())))
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

/// Check that the specified times are valid. Returns Ok if valid, else ReadXdrError.
pub fn sanity_check_timerange(start_time: f32, end_time: f32) -> Result<(), ReadXdrError> {
    if start_time < 0.0 {
        return Err(ReadXdrError::TimeRangeNegative(start_time.to_string()));
    }

    if end_time < 0.0 {
        return Err(ReadXdrError::TimeRangeNegative(end_time.to_string()));
    }

    if start_time > end_time {
        return Err(ReadXdrError::InvalidTimeRange(
            start_time.to_string(),
            end_time.to_string(),
        ));
    }

    Ok(())
}

/***********************************/
/*  Traits for reading and writing */
/***********************************/

/// Any structure implementing `XdrReader` can be used to read an xdr file.
pub trait XdrReader<'a>: Iterator<Item = Result<&'a mut System, ReadXdrError>> {
    /// Open an xdr (= xtc / trr) file creating an iterator over it.
    ///
    /// ## Example
    /// Using `XdrReader::new` in a generic function.
    /// ```no_run
    /// use groan_rs::prelude::*;
    /// use std::path::Path;
    ///
    /// // this function can read any xtc or trr file
    /// fn read_xdr_file<'a, Reader>(system: &'a mut System, file: impl AsRef<Path>)
    ///     where Reader: XdrReader<'a>
    /// {
    ///     // open the xtc/trr file for reading
    ///     let iterator = Reader::new(system, file).unwrap();
    ///
    ///     // read the xtc/trr file
    ///     for raw_frame in iterator {
    ///         let frame = raw_frame.unwrap();
    ///
    ///         // perform some operation with frame
    ///     }
    /// }
    ///
    /// // `read_xdr_file` can be then called to read either an xtc file or a trr file
    /// let mut system = System::from_file("system.gro").unwrap();
    /// read_xdr_file::<XtcReader>(&mut system, "trajectory.xtc");
    /// read_xdr_file::<TrrReader>(&mut system, "trajectory.trr");
    /// ```
    fn new(system: &'a mut System, filename: impl AsRef<Path>) -> Result<Self, ReadXdrError>
    where
        Self: Sized;
}

/// Any structure that implements the XdrRangeReader interface can be utilized for performing partial
/// reads of an xdr file, utilizing the provided time range.
pub trait XdrRangeReader<'a>: Iterator<Item = Result<&'a mut System, ReadXdrError>> {
    /// Open an xdr (xtc/trr) file and create an iterator for a specified range of frames.
    /// The `start_time` and `end_time` parameters should be provided in picoseconds.
    fn new(
        system: &'a mut System,
        filename: impl AsRef<Path>,
        start_time: f32,
        end_time: f32,
    ) -> Result<Self, ReadXdrError>
    where
        Self: Sized;
}

/// Any structure implementing `XdrWriter` can be used to write an xdr file.
pub trait XdrWriter {
    /// Open a new xdr (= xtc or trr) file for writing.
    ///
    /// ## Example
    /// Using `XdrWriter::new` in a generic function.
    /// ```no_run
    /// use groan_rs::prelude::*;
    /// use std::path::Path;
    ///
    /// // this function can write an xtc or a trr file
    /// fn write_xdr_file<Writer>(system: &System, file: impl AsRef<Path>)
    ///     where Writer: XdrWriter
    /// {
    ///     // open the xtc/trr file for writing
    ///     let mut writer = Writer::new(system, file).unwrap();
    ///
    ///     // write frame into the xtc/trr file
    ///     writer.write_frame().unwrap();
    /// }
    ///
    /// // `write_xdr_file` can be then called to write either an xtc file or a trr file
    /// let system = System::from_file("system.gro").unwrap();
    /// write_xdr_file::<XtcWriter>(&system, "trajectory.xtc");
    /// write_xdr_file::<TrrWriter>(&system, "trajectory.trr");
    /// ```
    fn new(system: &System, filename: impl AsRef<Path>) -> Result<Self, WriteXdrError>
    where
        Self: Sized;

    /// Write the current state of the system into an open xdr file (xtc or trr).
    fn write_frame(&mut self) -> Result<(), WriteXdrError>;
}

/// Any structure implementing `XdrGroupWriter` can be used to write an xdr file.
pub trait XdrGroupWriter {
    /// Open a new xdr file for writing and associate a specific group from a specific system with it.
    ///
    /// ## Example
    /// Using `XdrGroupWriter::new` in a generic function.
    /// ```no_run
    /// use groan_rs::prelude::*;
    /// use std::path::Path;
    ///
    /// // this function can write an xtc or a trr file
    /// fn write_xdr_file<Writer>(system: &System, file: impl AsRef<Path>)
    ///     where Writer: XdrGroupWriter
    /// {
    ///     // open the xtc/trr file for writing
    ///     let mut writer = Writer::new(system, "Protein", file).unwrap();
    ///
    ///     // write frame into the xtc/trr file
    ///     writer.write_frame().unwrap();
    /// }
    ///
    /// // `write_xdr_file` can be then called to write the group `Protein` either into an xtc file or a trr file
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("Protein", "@protein").unwrap();
    /// write_xdr_file::<XtcGroupWriter>(&system, "trajectory.xtc");
    /// write_xdr_file::<TrrGroupWriter>(&system, "trajectory.trr");
    /// ```
    fn new(system: &System, group: &str, filename: impl AsRef<Path>) -> Result<Self, WriteXdrError>
    where
        Self: Sized;

    /// Write the current state of the specified group into an open xdr file (xtc or trr).
    fn write_frame(&mut self) -> Result<(), WriteXdrError>;
}
