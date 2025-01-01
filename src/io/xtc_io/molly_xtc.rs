// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Bindings to the `molly` crate.

use molly::{
    selection::{AtomSelection, FrameSelection},
    Frame,
};

use crate::{
    errors::ReadTrajError,
    prelude::{
        FrameData, FrameDataTime, TrajFile, TrajRangeRead, TrajRead, TrajReadOpen, TrajStepRead,
        TrajStepTimeRead, Vector3D,
    },
    system::System,
};
use std::{fs::File, io::ErrorKind, marker::PhantomData, path::Path};

/**************************/
/*       READING XTC      */
/**************************/

/// Used when jumping to the start of iteration.
const TIME_PRECISION: f32 = 0.001;

/// Structure for storing the current offset to use.
#[derive(Debug, Clone)]
struct OffsetIterator {
    counter: usize,
    offsets: Box<[u64]>,
}

impl Iterator for OffsetIterator {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        self.offsets.get(self.counter).map(|&x| {
            self.counter += 1;
            x
        })
    }
}

/// Wrapper around `molly`'s XTCReader.
pub struct MollyXtc {
    reader: molly::XTCReader<File>,
    /// Only set if needed, i.e. if we are doing frame or atom selection.
    offsets: Option<OffsetIterator>,
    frame_selection: FrameSelection,
    atom_selection: AtomSelection,
}

impl MollyXtc {
    /// Open a new xtc file for reading.
    fn new(filename: impl AsRef<Path>) -> Result<MollyXtc, ReadTrajError> {
        Ok(MollyXtc {
            reader: molly::XTCReader::open(&filename)
                .map_err(|_| ReadTrajError::FileNotFound(Box::from(filename.as_ref())))?,
            offsets: None,
            frame_selection: FrameSelection::All,
            atom_selection: AtomSelection::All,
        })
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

    fn from_frame(traj_file: &mut MollyXtc, _system: &System) -> Option<Result<Self, ReadTrajError>>
    where
        Self: Sized,
    {
        if let Some(offsets) = traj_file.offsets.as_mut() {
            // read the next frame using offset
            let offset = offsets.next()?;
            let mut frame = Frame::default();

            match traj_file.reader.read_frame_at_offset::<true>(
                &mut frame,
                offset,
                &traj_file.atom_selection,
            ) {
                Ok(_) => Some(Ok(frame)),
                Err(e) => Some(Err(ReadTrajError::MollyXtcError(e.to_string()))),
            }
        } else {
            // read the next frame naively
            let mut frame = Frame::default();

            match traj_file.reader.read_frame(&mut frame) {
                Ok(_) => Some(Ok(frame)),
                // expecting that this is not an error but end of file was just reached
                Err(e) if e.kind() == ErrorKind::UnexpectedEof => None,
                Err(e) => Some(Err(ReadTrajError::MollyXtcError(e.to_string()))),
            }
        }
    }

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
        let xtc = MollyXtc::new(&filename)?;

        Ok(XtcReader {
            system,
            xtc,
            phantom: PhantomData,
        })
    }
}

impl<'a> TrajRangeRead<'a> for XtcReader<'a> {
    fn jump_to_start(&mut self, start_time: f32) -> Result<(), ReadTrajError> {
        // TODO: optimize this; currently this reads all frames
        let system = self.get_system();
        loop {
            // funnily, we do not even need the system wrapped inside unsafe, but the trait requires it...
            let frame = match Frame::from_frame(&mut self.xtc, unsafe { &*system }) {
                None => return Err(ReadTrajError::StartNotFound(start_time.to_string())),
                Some(Err(e)) => return Err(e),
                Some(Ok(x)) => x,
            };

            if frame.get_time() >= start_time - TIME_PRECISION {
                return Ok(());
            }
        }
    }
}

impl<'a> TrajStepRead<'a> for XtcReader<'a> {
    fn skip_frame(&mut self) -> Result<bool, ReadTrajError> {
        // TODO: optimize this; currently this reads all frames
        let system = self.get_system();
        match Frame::from_frame(&mut self.xtc, unsafe { &*system }) {
            None => Ok(false),
            Some(Err(e)) => Err(e),
            Some(Ok(_)) => Ok(true),
        }
    }
}

impl<'a> TrajStepTimeRead<'a> for XtcReader<'a> {
    fn skip_frame_time(&mut self) -> Result<Option<f32>, ReadTrajError> {
        // TODO: optimize this; currently, this reads all frames
        let system = self.get_system();
        match Frame::from_frame(&mut self.xtc, unsafe { &*system }) {
            None => Ok(None),
            Some(Err(e)) => Err(e),
            Some(Ok(f)) => Ok(Some(f.get_time())),
        }
    }
}

// TODO:
// Add to molly: jump_to_start function, skip_frame function, skip_frame_time function
// Remove offsets from MollyXtc and implement changes that follow. Remove Frame selection from MollyXtc.
// Implement atom number checks.
// Implement partial trajectory reader using atom selections.
// Implement xtc file writing using molly.
// Create a feature `xdrfile-trr` and separate trr file reading and writing to it.
// Make `molly` feature default. If `molly` is not specified, `xdrfile` is used.
// Make build.rs dependent on xdrfile feature.
// Reflect changes in the documentation (e.g. support for new magic number).
