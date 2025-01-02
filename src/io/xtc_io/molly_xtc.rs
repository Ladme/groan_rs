// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Bindings to the `molly` crate.

use molly::{selection::AtomSelection, Frame};

use crate::{
    errors::ReadTrajError,
    prelude::{
        FrameData, FrameDataTime, TrajFile, TrajRangeRead, TrajRead, TrajReadOpen, TrajStepRead,
        TrajStepTimeRead, Vector3D,
    },
    system::System,
};
use std::{
    fs::File,
    io::{ErrorKind, Read, Seek, SeekFrom},
    marker::PhantomData,
    path::Path,
};

/**************************/
/*       READING XTC      */
/**************************/

/// Used when jumping to the start of iteration.
const TIME_PRECISION: f32 = 0.001;

/// Wrapper around `molly`'s XTCReader.
pub struct MollyXtc {
    reader: molly::XTCReader<File>,
    _atom_selection: AtomSelection,
}

impl MollyXtc {
    /// Open a new xtc file for reading.
    /// Checks that the number of atoms in the XTC file matches the number of atoms in the System.
    /// (This should compare the total number of atoms in the system.)
    fn new(filename: impl AsRef<Path>, n_atoms: usize) -> Result<MollyXtc, ReadTrajError> {
        let file = std::fs::File::open(&filename)
            .map_err(|_| ReadTrajError::FileNotFound(Box::from(filename.as_ref())))?;

        let mut xtc = MollyXtc {
            reader: molly::XTCReader { file, step: 0 },
            _atom_selection: AtomSelection::All,
        };

        // check magic number
        match xtc.read_i32() {
            Ok(x) if x == 1995 || x == 2023 => x,
            Ok(_) => return Err(ReadTrajError::NotXtc(Box::from(filename.as_ref()))),
            Err(_) => return Err(ReadTrajError::CouldNotReadMagic),
        };

        // check that the numbers of atoms match
        match xtc.read_i32() {
            Ok(x) if x == n_atoms as i32 => (),
            Ok(_) => {
                return Err(ReadTrajError::AtomsNumberMismatch(Box::from(
                    filename.as_ref(),
                )))
            }
            Err(_) => return Err(ReadTrajError::FrameNotFound),
        };

        // return back to the start
        xtc.reader
            .home()
            .expect("FATAL GROAN ERROR | MollyXtc::new | Could not jump to the start of the file.");

        Ok(xtc)
    }

    /// Skip one trajectory frame of an xtc file.
    fn skip_frame(&mut self) -> Result<bool, ReadTrajError> {
        // read magic number
        let magic_number = match self.read_i32() {
            Ok(x) => x,
            // this should only fail if we have reached the end of the file
            Err(_) => return Ok(false),
        };

        // check the magic number is correct (validity of the frame)
        if magic_number != 1995 && magic_number != 2023 {
            return Err(ReadTrajError::FrameNotFound);
        }

        // get the number of bytes to the next frame
        self.xdr_jump(84).map_err(|_| ReadTrajError::SkipFailed)?;
        let mut size = self.read_i32().map_err(|_| ReadTrajError::SkipFailed)? as i64;

        // add padding
        if size % 4 != 0 {
            size += 4 - (size % 4);
        }

        // this should only fail if we have reached the end of the file
        if self.xdr_jump(size).is_err() {
            return Ok(false);
        }

        Ok(true)
    }

    /// Skip one trajectory frame of an xtc file and report the simulation time of this frame.
    fn skip_frame_with_time(&mut self) -> Result<Option<f32>, ReadTrajError> {
        // read magic number
        let magic_number = match self.read_i32() {
            Ok(x) => x,
            // this should only fail if we have reached the end of the file
            Err(_) => return Ok(None),
        };

        // check the magic number is correct
        if magic_number != 1995 && magic_number != 2023 {
            return Err(ReadTrajError::FrameNotFound);
        }

        // jump to time information and read it
        self.xdr_jump(8).map_err(|_| ReadTrajError::SkipFailed)?;
        let time = self.read_f32().map_err(|_| ReadTrajError::SkipFailed)?;

        // get the number of bytes to the next frame
        self.xdr_jump(72).map_err(|_| ReadTrajError::SkipFailed)?;
        let mut size = self.read_i32().map_err(|_| ReadTrajError::SkipFailed)? as i64;

        // add padding
        if size % 4 != 0 {
            size += 4 - (size % 4);
        }

        // this should only fail if we have reached the end of the file
        if self.xdr_jump(size).is_err() {
            return Ok(None);
        }

        Ok(Some(time))
    }

    /// Jump to the start of the iteration.
    fn jump_to_start(&mut self, target_time: f32) -> Result<(), ReadTrajError> {
        // repeatedly jump over frames until the target frame is reached
        loop {
            match self.get_jump_info(target_time, TIME_PRECISION)? {
                0 => break, // target time reached
                jump => self.xdr_jump(jump).map_err(|_| ReadTrajError::SkipFailed)?,
            }
        }

        // jump back to the start of the frame
        self.xdr_jump(-16).map_err(|_| ReadTrajError::SkipFailed)?;

        Ok(())
    }

    /// Get the number of bytes required to jump to the next frame in an xtc file.
    /// Use in [`MollyXtc::jump_to_start`].
    fn get_jump_info(
        &mut self,
        target_time: f32,
        time_precision: f32,
    ) -> Result<i64, ReadTrajError> {
        // read magic number
        let magic_number = self
            .read_i32()
            .map_err(|_| ReadTrajError::StartNotFound(target_time.to_string()))?;

        // check that the magic number is correct
        if magic_number != 1995 && magic_number != 2023 {
            return Err(ReadTrajError::FrameNotFound);
        }

        // jump to time information and read it
        self.xdr_jump(8)
            .map_err(|_| ReadTrajError::StartNotFound(target_time.to_string()))?;
        let time = self
            .read_f32()
            .map_err(|_| ReadTrajError::StartNotFound(target_time.to_string()))?;

        // check whether the time is higher than or equal to the target time
        if time >= target_time - time_precision {
            return Ok(0);
        }

        // get the number of bytes to the next frame
        self.xdr_jump(72)
            .map_err(|_| ReadTrajError::StartNotFound(target_time.to_string()))?;
        let mut size = self
            .read_i32()
            .map_err(|_| ReadTrajError::StartNotFound(target_time.to_string()))?
            as i64;

        // Add padding
        if size % 4 != 0 {
            size += 4 - (size % 4);
        }

        Ok(size)
    }

    /// Read 32-bit integer from the xtc file.
    /// (Taken from `molly`.)
    #[inline(always)]
    fn read_i32(&mut self) -> std::io::Result<i32> {
        let mut buf = [0u8; 4];
        self.reader.file.read_exact(&mut buf)?;
        Ok(i32::from_be_bytes(buf))
    }

    /// Read 32-bit float from the xtc file.
    /// (Taken from `molly`.)
    #[inline(always)]
    fn read_f32(&mut self) -> std::io::Result<f32> {
        let mut buf = [0u8; 4];
        self.reader.file.read_exact(&mut buf)?;
        Ok(f32::from_be_bytes(buf))
    }

    /// Jump forward in the open xtc file.
    #[inline(always)]
    fn xdr_jump(&mut self, offset: i64) -> std::io::Result<()> {
        self.reader.file.seek(SeekFrom::Current(offset))?;
        Ok(())
    }
}

/// Iterator over an xtc file.
pub struct XtcReader<'a> {
    system: *mut System,
    xtc: MollyXtc,
    phantom: PhantomData<&'a mut System>,
}

impl TrajFile for MollyXtc {}

impl FrameData for molly::Frame {
    type TrajFile = MollyXtc;

    #[inline(always)]
    fn from_frame(traj_file: &mut MollyXtc, _system: &System) -> Option<Result<Self, ReadTrajError>>
    where
        Self: Sized,
    {
        let mut frame = Frame::default();

        match traj_file.reader.read_frame(&mut frame) {
            Ok(_) => Some(Ok(frame)),
            // expecting that this is not an error but the end of file was just reached
            Err(e) if e.kind() == ErrorKind::UnexpectedEof => None,
            Err(e) => Some(Err(ReadTrajError::MollyXtcError(e.to_string()))),
        }
    }

    #[inline]
    fn update_system(self, system: &mut System) {
        let positions = self.coords();
        for (pos, atom) in positions.zip(system.get_atoms_mut().iter_mut()) {
            atom.set_position(Vector3D::new(pos.x, pos.y, pos.z));
            atom.reset_velocity();
            atom.reset_force();
        }

        // update the system
        system.set_simulation_step(self.step as u64);
        system.set_simulation_time(self.time);
        let b = self.boxvec;
        system.set_box(
            [
                b.col(0).x,
                b.col(1).y,
                b.col(2).z,
                b.col(0).y,
                b.col(0).z,
                b.col(1).x,
                b.col(1).z,
                b.col(2).x,
                b.col(2).y,
            ]
            .into(),
        );
        system.set_precision(self.precision as u64);
    }
}

impl FrameDataTime for molly::Frame {
    fn get_time(&self) -> f32 {
        self.time
    }
}

impl<'a> TrajRead<'a> for XtcReader<'a> {
    type FrameData = molly::Frame;

    fn get_system(&mut self) -> *mut System {
        self.system
    }

    fn get_file_handle(&mut self) -> &mut MollyXtc {
        &mut self.xtc
    }
}

impl<'a> TrajReadOpen<'a> for XtcReader<'a> {
    /// Create an iterator over an xtc file.
    fn new(system: &'a mut System, filename: impl AsRef<Path>) -> Result<XtcReader, ReadTrajError> {
        let xtc = MollyXtc::new(&filename, system.get_n_atoms())?;

        Ok(XtcReader {
            system,
            xtc,
            phantom: PhantomData,
        })
    }
}

impl<'a> TrajRangeRead<'a> for XtcReader<'a> {
    fn jump_to_start(&mut self, start_time: f32) -> Result<(), ReadTrajError> {
        self.xtc.jump_to_start(start_time)
    }
}

impl<'a> TrajStepRead<'a> for XtcReader<'a> {
    #[inline(always)]
    fn skip_frame(&mut self) -> Result<bool, ReadTrajError> {
        self.xtc.skip_frame()
    }
}

impl<'a> TrajStepTimeRead<'a> for XtcReader<'a> {
    #[inline(always)]
    fn skip_frame_time(&mut self) -> Result<Option<f32>, ReadTrajError> {
        self.xtc.skip_frame_with_time()
    }
}

// TODO:
// [x] Add to molly: jump_to_start function, skip_frame function, skip_frame_time function
// [x] Remove offsets from MollyXtc and implement changes that follow. Remove Frame selection from MollyXtc.
// Check MAGIC NUMBER 2023: number of atoms and number of bytes
// Implement atom number checks.
// Implement partial trajectory reader using atom selections.
// Implement xtc file writing using molly.
// Create a feature `xdrfile-trr` and separate trr file reading and writing to it.
// Make `molly` feature default. If `molly` is not specified, `xdrfile` is used.
// Make build.rs dependent on xdrfile feature.
// Reflect changes in the documentation (e.g. support for new magic number).
