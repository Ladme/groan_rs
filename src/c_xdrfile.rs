// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Rust bindings for the xdrfile library.

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct XDRFILE {
    _empty: [u8; 0],
}

extern "C" {
    pub fn xdrfile_open(
        path: *const ::std::os::raw::c_char,
        mode: *const ::std::os::raw::c_char,
    ) -> *mut XDRFILE;
}

extern "C" {
    pub fn xdrfile_close(xfp: *mut XDRFILE) -> ::std::os::raw::c_int;
}

extern "C" {
    pub fn read_xtc(
        xd: *mut XDRFILE,
        natoms: ::std::os::raw::c_int,
        step: *mut ::std::os::raw::c_int,
        time: *mut ::std::os::raw::c_float,
        box_vec: *mut [[::std::os::raw::c_float; 3usize]; 3usize],
        x: *mut [::std::os::raw::c_float; 3usize],
        prec: *mut ::std::os::raw::c_float,
    ) -> ::std::os::raw::c_int;
}

extern "C" {
    pub fn read_xtc_natoms(
        filename: *const ::std::os::raw::c_char,
        natoms: *mut ::std::os::raw::c_int,
    ) -> ::std::os::raw::c_int;
}

extern "C" {
    pub fn write_xtc(
        xd: *mut XDRFILE,
        natoms: ::std::os::raw::c_int,
        step: ::std::os::raw::c_int,
        time: ::std::os::raw::c_float,
        box_vec: *mut [[::std::os::raw::c_float; 3usize]; 3usize],
        x: *mut [::std::os::raw::c_float; 3usize],
        prec: ::std::os::raw::c_float,
    ) -> ::std::os::raw::c_int;
}
