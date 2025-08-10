// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Implementation of functions for reading and writing xtc files.

use std::marker::PhantomData;
use std::os::raw::{c_float, c_int};
use std::path::Path;

use crate::errors::{ReadTrajError, TrajError};
use crate::io::traj_read::{
    FrameData, FrameDataTime, TrajFullReadOpen, TrajRangeRead, TrajRead, TrajReadOpen,
    TrajStepRead, TrajStepTimeRead,
};
use crate::io::xdrfile::{self, OpenMode, XdrFile};
use crate::structures::{simbox::SimBox, vector3d::Vector3D};
use crate::system::System;

/**************************/
/*       READING XTC      */
/**************************/

/// Structure containing data read from an xtc file.
#[derive(Debug)]
pub struct XtcFrameData {
    step: c_int,
    time: c_float,
    simbox: SimBox,
    precision: c_float,
    coordinates: Vec<[f32; 3]>,
}

impl FrameData for XtcFrameData {
    type TrajFile = XdrFile;

    /// Read `XtcFrameData` from an xtc frame.
    ///
    /// ## Returns
    /// - `None` if the file has been read completely
    /// - `Some(XtcFrameData)` if the frame has been successfully read.
    /// - `Some(ReadTrajError)` if the frame could not be read.
    fn from_frame(
        xdrfile: &mut Self::TrajFile,
        system: &System,
    ) -> Option<Result<Self, ReadTrajError>> {
        let mut step: c_int = 0;
        let mut time: c_float = 0.0;
        let mut boxvector: [[c_float; 3usize]; 3usize] = [[0.0; 3]; 3];
        let mut precision: c_float = 0.0;
        let mut coordinates = vec![[0.0, 0.0, 0.0]; system.get_n_atoms()];

        unsafe {
            let return_code = xdrfile::xtc::read_xtc(
                xdrfile.handle,
                system.get_n_atoms() as c_int,
                &mut step,
                &mut time,
                &mut boxvector,
                coordinates.as_mut_ptr(),
                &mut precision,
            );

            match return_code {
                // reading successful and there is more to read
                0 => {
                    let simbox = match xdrfile::matrix2simbox(boxvector) {
                        Ok(x) => x,
                        Err(e) => return Some(Err(e)),
                    };

                    Some(Ok(XtcFrameData {
                        step,
                        time,
                        simbox,
                        precision,
                        coordinates,
                    }))
                }
                // file is completely read
                11 => None,
                // error occured
                _ => Some(Err(ReadTrajError::FrameNotFound)),
            }
        }
    }

    /// Update the `System` structure based on data from `XtcFrameData`.
    fn update_system(self, system: &mut System) {
        unsafe {
            for (i, atom) in system.get_atoms_mut().iter_mut().enumerate() {
                atom.set_position(Vector3D::from(*self.coordinates.get_unchecked(i)));
                atom.reset_velocity();
                atom.reset_force();
            }

            // update the system
            system.set_simulation_step(self.step as u64);
            system.set_simulation_time(self.time);
            system.set_box(self.simbox);
            system.set_precision(self.precision as u64);
        }
    }
}

impl FrameDataTime for XtcFrameData {
    fn get_time(&self) -> f32 {
        self.time
    }
}

/// Iterator over an xtc file.
pub struct XtcReader<'a> {
    system: *mut System,
    xtc: XdrFile,
    phantom: PhantomData<&'a mut System>,
}

impl<'a> TrajRead<'a> for XtcReader<'a> {
    type FrameData = XtcFrameData;

    fn get_system(&mut self) -> *mut System {
        self.system
    }

    fn get_file_handle(&mut self) -> &mut XdrFile {
        &mut self.xtc
    }
}

impl<'a> TrajReadOpen<'a> for XtcReader<'a> {
    /// Create an iterator over an xtc file.
    ///
    /// ## Panic
    /// Panics if the `group` is **not** None.
    ///
    /// ## Note
    /// Prefer using [`XtcReader::new`] which does not panic.
    fn initialize(
        system: &'a mut System,
        filename: impl AsRef<Path>,
        group: Option<&str>,
    ) -> Result<Self, ReadTrajError>
    where
        Self: Sized,
    {
        match group {
            None => XtcReader::new(system, filename),
            Some(_) => panic!("FATAL GROAN ERROR | XtcReader::initialize | XtcReader does not support partial-frame reading. You may want to enable `molly` feature and use `GroupXtcReader`."),
        }
    }
}

impl<'a> TrajFullReadOpen<'a> for XtcReader<'a> {
    /// Create an iterator over an xtc file.
    fn new(
        system: &'a mut System,
        filename: impl AsRef<Path>,
    ) -> Result<XtcReader<'a>, ReadTrajError> {
        let n_atoms = system.get_n_atoms();

        // sanity check the number of atoms
        match XdrFile::check_xtc(filename.as_ref(), n_atoms) {
            Err(e) => return Err(e),
            Ok(false) => {
                return Err(ReadTrajError::AtomsNumberMismatch(Box::from(
                    filename.as_ref(),
                )))
            }
            Ok(true) => (),
        };

        // open the xtc file and save the handle to it
        let xtc = match XdrFile::open_xdr(filename.as_ref(), OpenMode::Read) {
            Ok(x) => x,
            Err(TrajError::FileNotFound(x)) => return Err(ReadTrajError::FileNotFound(x)),
            Err(TrajError::InvalidPath(x)) => return Err(ReadTrajError::InvalidPath(x)),
        };

        let xtc_reader = XtcReader {
            system,
            xtc,
            phantom: PhantomData,
        };

        Ok(xtc_reader)
    }
}

impl<'a> TrajRangeRead<'a> for XtcReader<'a> {
    fn jump_to_start(&mut self, start_time: f32) -> Result<(), ReadTrajError> {
        unsafe {
            if xdrfile::xtc::xtc_jump_to_start(self.get_file_handle().handle, start_time as c_float)
                != 0
            {
                Err(ReadTrajError::StartNotFound(start_time.to_string()))
            } else {
                Ok(())
            }
        }
    }
}

impl<'a> TrajStepRead<'a> for XtcReader<'a> {
    fn skip_frame(&mut self) -> Result<bool, ReadTrajError> {
        unsafe {
            match xdrfile::xtc::xtc_skip_frame(self.get_file_handle().handle) {
                0 => Ok(true),
                1 => Err(ReadTrajError::SkipFailed),
                2 => Ok(false),
                number => panic!(
                    "FATAL GROAN ERROR | XtcReader::skip_frame | `xdrfile::xtc_skip_frame` returned an unsupported number '{}'",
                    number
                ),
            }
        }
    }
}

impl<'a> TrajStepTimeRead<'a> for XtcReader<'a> {
    fn skip_frame_time(&mut self) -> Result<Option<f32>, ReadTrajError> {
        unsafe {
            let mut time: c_float = 0.0;

            match xdrfile::xtc::xtc_skip_frame_with_time(self.get_file_handle().handle, &mut time as *mut c_float) {
                0 => Ok(Some(time)),
                1 => Err(ReadTrajError::SkipFailed),
                2 => Ok(None),
                number => panic!(
                    "FATAL GROAN ERROR | XtcReader::skip_frame_time | `xdrfile::xtc_skip_frame_with_time` returned an unsupported number '{}'",
                    number
                ),
            }
        }
    }
}

/******************************/
/*     PRIVATE FUNCTIONS      */
/******************************/

impl XdrFile {
    /// Check that the number of atoms in an unopened xtc file matches the expected number.
    fn check_xtc(filename: impl AsRef<Path>, n_atoms: usize) -> Result<bool, ReadTrajError> {
        unsafe {
            let c_path = match xdrfile::path2cstring(filename.as_ref()) {
                Ok(x) => x,
                Err(_) => return Err(ReadTrajError::InvalidPath(Box::from(filename.as_ref()))),
            };

            let mut xtc_atoms: c_int = 0;

            if xdrfile::xtc::read_xtc_natoms(c_path.as_ptr(), &mut xtc_atoms) != 0 {
                // reading the file failed
                return Err(ReadTrajError::FileNotFound(Box::from(filename.as_ref())));
            }

            // if reading was successful
            if xtc_atoms == n_atoms as c_int {
                Ok(true)
            } else {
                Ok(false)
            }
        }
    }
}
