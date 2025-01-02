// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Bindings to the `molly` crate.

use molly::{selection::AtomSelection, Frame};

use crate::{
    errors::ReadTrajError,
    io::traj_read::TrajGroupReadOpen,
    prelude::{
        FrameData, FrameDataTime, TrajFile, TrajFullReadOpen, TrajRangeRead, TrajRead,
        TrajReadOpen, TrajReader, TrajStepRead, TrajStepTimeRead, Vector3D,
    },
    structures::{container::AtomContainer, group::Group},
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
    group: AtomContainer,
    atom_selection: AtomSelection,
}

impl MollyXtc {
    /// Open a new xtc file for reading.
    /// Checks that the number of atoms in the XTC file matches the number of atoms in the System.
    /// (This should compare the total number of atoms in the system.)
    fn new(
        filename: impl AsRef<Path>,
        n_atoms: usize,
        group: Option<&Group>,
    ) -> Result<MollyXtc, ReadTrajError> {
        let file = std::fs::File::open(&filename)
            .map_err(|_| ReadTrajError::FileNotFound(Box::from(filename.as_ref())))?;

        let mut xtc = match group {
            None => MollyXtc {
                reader: molly::XTCReader { file, step: 0 },
                group: AtomContainer::from_ranges(vec![(0, n_atoms)], n_atoms),
                atom_selection: AtomSelection::All,
            },
            Some(group) => {
                let number = group.get_atoms().last().map(|x| x + 1).unwrap_or(0);
                MollyXtc {
                    reader: molly::XTCReader { file, step: 0 },
                    group: group.get_atoms().clone(),
                    atom_selection: AtomSelection::Until(number as u32),
                }
            }
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

    /// Read the number of bytes allocated for the coordinates in the current frame + PADDING.
    /// Reads 4 bytes for magic number 1995 and 8 bytes for magic number 2023.
    #[inline]
    fn read_n_bytes(&mut self, magic: i32) -> Result<i64, ReadTrajError> {
        let bytes = match magic {
            1995 => self.read_i32().map_err(|_| ReadTrajError::SkipFailed)? as i64,
            2023 => self.read_i64().map_err(|_| ReadTrajError::SkipFailed)?,
            _ => panic!("FATAL GROAN ERROR | MollyXtc::read_n_bytes | Unexpected magic number `{}` provided.", magic),
        };

        // add padding
        Ok(bytes + (4 - (bytes % 4)) % 4)
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
        let size = self.read_n_bytes(magic_number)?;

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
        let size = self.read_n_bytes(magic_number)?;

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
        let size = self.read_n_bytes(magic_number)?;

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

    /// Read 64-bit integer from the xtc file.
    pub(crate) fn read_i64(&mut self) -> std::io::Result<i64> {
        let mut buf: [u8; 8] = Default::default();
        self.reader.file.read_exact(&mut buf)?;
        Ok(i64::from_be_bytes(buf))
    }

    /// Jump forward in the open xtc file.
    #[inline(always)]
    fn xdr_jump(&mut self, offset: i64) -> std::io::Result<()> {
        self.reader.file.seek(SeekFrom::Current(offset))?;
        Ok(())
    }
}

/// Iterator over an xtc file.
/// Constructed using [`System::xtc_iter`].
pub struct XtcReader<'a> {
    system: *mut System,
    xtc: MollyXtc,
    phantom: PhantomData<&'a mut System>,
}

impl TrajFile for MollyXtc {}

/// Wrapper around `molly`'s Frame.
pub struct XtcFrameData {
    frame: molly::Frame,
    /// Pointer to a copy of an atom container stored in `MollyXtc`.
    /// Safety: XtcFrameData only exists for a short time during iteration of a trajectory. It does not leaks outside.
    group: *const AtomContainer,
}

impl FrameData for XtcFrameData {
    type TrajFile = MollyXtc;

    #[inline(always)]
    fn from_frame(traj_file: &mut MollyXtc, _system: &System) -> Option<Result<Self, ReadTrajError>>
    where
        Self: Sized,
    {
        let mut frame = Frame::default();

        match traj_file
            .reader
            .read_frame_with_selection(&mut frame, &traj_file.atom_selection)
        {
            Ok(_) => {
                let frame_data = XtcFrameData {
                    frame,
                    group: &traj_file.group as *const AtomContainer,
                };

                Some(Ok(frame_data))
            }
            // expecting that this is not an error but the end of file was just reached
            Err(e) if e.kind() == ErrorKind::UnexpectedEof => None,
            Err(e) => Some(Err(ReadTrajError::MollyXtcError(e.to_string()))),
        }
    }

    #[inline]
    fn update_system(self, system: &mut System) {
        let positions = self.frame.coords();
        // safety: `XtcFrameData` only exists during iteration
        // group is obtained from `MollyXtc` which is located inside `XtcReader` which must by definition
        // live through the entire iteration
        let atoms = unsafe { &*self.group };

        for (index, (pos, atom)) in positions.zip(system.get_atoms_mut().iter_mut()).enumerate() {
            if atoms.isin(index) {
                atom.set_position(Vector3D::new(pos.x, pos.y, pos.z));
                atom.reset_velocity();
                atom.reset_force();
            }
            // atoms that are not part of the group keep their original properties
        }

        // update the system
        system.set_simulation_step(self.frame.step as u64);
        system.set_simulation_time(self.frame.time);
        let b = self.frame.boxvec;
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
        system.set_precision(self.frame.precision as u64);
    }
}

impl FrameDataTime for XtcFrameData {
    #[inline(always)]
    fn get_time(&self) -> f32 {
        self.frame.time
    }
}

impl<'a> TrajRead<'a> for XtcReader<'a> {
    type FrameData = XtcFrameData;

    #[inline(always)]
    fn get_system(&mut self) -> *mut System {
        self.system
    }

    #[inline(always)]
    fn get_file_handle(&mut self) -> &mut MollyXtc {
        &mut self.xtc
    }
}

impl<'a> TrajFullReadOpen<'a> for XtcReader<'a> {
    /// Create an iterator over an xtc file.
    #[inline(always)]
    fn new(system: &'a mut System, filename: impl AsRef<Path>) -> Result<XtcReader, ReadTrajError> {
        let xtc = MollyXtc::new(&filename, system.get_n_atoms(), None)?;

        Ok(XtcReader {
            system,
            xtc,
            phantom: PhantomData,
        })
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
    #[inline(always)]
    fn initialize(
        system: &'a mut System,
        filename: impl AsRef<Path>,
        group: Option<&str>,
    ) -> Result<Self, ReadTrajError>
    where
        Self: Sized,
    {
        match group {
            None => TrajFullReadOpen::new(system, filename),
            Some(_) => panic!("FATAL GROAN ERROR | XtcReader::initialize | XtcReader does not support partial-frame reading. Use `GroupXtcReader` instead."),
        }
    }
}

impl<'a> TrajRangeRead<'a> for XtcReader<'a> {
    #[inline(always)]
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

/// Partial iterator over an xtc file. Reads only specified atoms from each frame.
/// Constructed using [`System::group_xtc_iter`].
pub struct GroupXtcReader<'a> {
    system: *mut System,
    xtc: MollyXtc,
    phantom: PhantomData<&'a mut System>,
}

impl<'a> TrajRead<'a> for GroupXtcReader<'a> {
    type FrameData = XtcFrameData;

    #[inline(always)]
    fn get_system(&mut self) -> *mut System {
        self.system
    }

    #[inline(always)]
    fn get_file_handle(&mut self) -> &mut MollyXtc {
        &mut self.xtc
    }
}

impl<'a> TrajReadOpen<'a> for GroupXtcReader<'a> {
    /// Create an iterator over an xtc file that only reads information about atoms from a specific group.
    ///
    /// ## Panic
    /// Panics if the `group` is None.
    ///
    /// ## Note
    /// Prefer using [`GroupXtcReader::new`] which does not panic.
    #[inline(always)]
    fn initialize(
        system: &'a mut System,
        filename: impl AsRef<Path>,
        group: Option<&str>,
    ) -> Result<Self, ReadTrajError>
    where
        Self: Sized,
    {
        match group {
            None => panic!("FATAL GROAN ERROR | GroupXtcReader::initialize | GroupXtcReader requires a group to be initialized. Use `XtcReader` instead."),
            Some(x) => TrajGroupReadOpen::new(system, filename, x),
        }
    }
}

impl<'a> TrajGroupReadOpen<'a> for GroupXtcReader<'a> {
    /// Create an iterator over an xtc file that only reads information about atoms from a specific group.
    #[inline(always)]
    fn new(
        system: &'a mut System,
        filename: impl AsRef<Path>,
        group: &str,
    ) -> Result<Self, ReadTrajError>
    where
        Self: Sized,
    {
        let group = match system.get_groups().get(group) {
            Some(x) => x,
            None => return Err(ReadTrajError::GroupNotFound(group.to_owned())),
        };

        let xtc = MollyXtc::new(&filename, system.get_n_atoms(), Some(group))?;
        Ok(GroupXtcReader {
            system,
            xtc,
            phantom: PhantomData,
        })
    }
}

impl<'a> TrajRangeRead<'a> for GroupXtcReader<'a> {
    #[inline(always)]
    fn jump_to_start(&mut self, start_time: f32) -> Result<(), ReadTrajError> {
        self.xtc.jump_to_start(start_time)
    }
}

impl<'a> TrajStepRead<'a> for GroupXtcReader<'a> {
    #[inline(always)]
    fn skip_frame(&mut self) -> Result<bool, ReadTrajError> {
        self.xtc.skip_frame()
    }
}

impl<'a> TrajStepTimeRead<'a> for GroupXtcReader<'a> {
    #[inline(always)]
    fn skip_frame_time(&mut self) -> Result<Option<f32>, ReadTrajError> {
        self.xtc.skip_frame_with_time()
    }
}

impl System {
    /// Create a `GroupXtcReader` structure which is an iterator over an xtc
    /// file that **efficiently** reads only properties of atoms of the specified group.
    /// If you haven't already, read the documentation of [`System::xtc_iter`] first.
    ///
    /// ## Returns
    /// `TrajReader<GroupXtcReader>` if the xtc file exists and matches the structure file.
    /// Else returns `ReadTrajError`.
    ///
    /// ## Examples
    /// Iterating through a trajectory while reading only the properties of protein atoms.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use groan_rs::errors::ReadTrajError;
    /// #
    /// // load system from file
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// // create the protein group
    /// system.group_create("Protein", "@protein").unwrap();
    ///
    /// // loop through the trajectory
    /// for raw_frame in system.group_xtc_iter("trajectory.xtc", "Protein").unwrap() {
    ///     let frame = raw_frame.unwrap();
    ///     println!("{:?}", frame.group_get_center("Protein"));
    /// }
    /// ```
    ///
    /// [`GroupXtcReader`] also supports range iteration using [`TrajReader::with_range`],
    /// stepping using [`TrajReader::with_step`], and progress printing using
    /// [`TrajMasterRead::print_progress`](crate::prelude::TrajMasterRead::print_progress).
    ///
    /// ## Notes
    /// - If the specified group does not exist, an error is returned.
    /// - The function checks whether the number of atoms in the **system** (not the group) corresponds to the number of atoms in the xtc file.
    /// - The `System` structure is modified while iterating through the xtc file. **Only properties of atoms
    ///   from the specified group are changed, all other atoms are left unchanged!**
    /// - The `velocity` and `force` information is set to `None` for all atoms from the group (but not for other atoms).
    /// - Group trajectory iteration can be extremely efficient thanks to the [`molly`] crate especially if your group is located close
    ///   to the beginning of the structure. If your group is located close to the end of the structure, the performance
    ///   gain compared to simple [`System::xtc_iter`] might be small (but still probably worth it).
    /// - When reading a group of zero size, properties of no atoms are changed, but properties of the system are updated.
    /// - Supports reading xtc files in the 2023 format (generally used for giant systems).
    ///   **Reading of the 2023 format is not well tested. Be careful!**
    pub fn group_xtc_iter(
        &mut self,
        filename: impl AsRef<Path>,
        group: &str,
    ) -> Result<TrajReader<GroupXtcReader>, ReadTrajError> {
        Ok(TrajReader::wrap_traj(GroupXtcReader::new(
            self, filename, group,
        )?))
    }
}

#[cfg(test)]
mod tests_group {
    use crate::{
        errors::ReadTrajError,
        system::System,
        test_utilities::utilities::{compare_atoms, compare_box},
    };

    #[test]
    fn group_xtc_iter_all() {
        let mut group_system = System::from_file("test_files/example.gro").unwrap();
        let mut full_system = System::from_file("test_files/example.gro").unwrap();

        for (group_frame, full_frame) in group_system
            .group_xtc_iter("test_files/short_trajectory.xtc", "all")
            .unwrap()
            .zip(
                full_system
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap(),
            )
        {
            let g = group_frame.unwrap();
            let f = full_frame.unwrap();

            compare_box(g.get_box().unwrap(), f.get_box().unwrap());

            for (a1, a2) in g.atoms_iter().zip(f.atoms_iter()) {
                compare_atoms(a1, a2);
            }
        }
    }

    #[test]
    fn group_xtc_iter_popc() {
        let mut group_system = System::from_file("test_files/example.gro").unwrap();
        group_system
            .group_create("Membrane", "resname POPC")
            .unwrap();
        let mut full_system = System::from_file("test_files/example.gro").unwrap();
        let original_system = full_system.clone();

        for (group_frame, full_frame) in group_system
            .group_xtc_iter("test_files/short_trajectory.xtc", "Membrane")
            .unwrap()
            .zip(
                full_system
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap(),
            )
        {
            let g = group_frame.unwrap();
            let f = full_frame.unwrap();

            compare_box(g.get_box().unwrap(), f.get_box().unwrap());

            for ((a1, a2), a3) in g
                .atoms_iter()
                .zip(f.atoms_iter())
                .zip(original_system.atoms_iter())
            {
                // only atoms of `Membrane` should be modified
                if g.group_isin("Membrane", a1.get_index()).unwrap() {
                    compare_atoms(a1, a2);
                } else {
                    compare_atoms(a1, a3);
                }
            }
        }
    }

    #[test]
    fn group_xtc_iter_protein() {
        let mut group_system = System::from_file("test_files/example.gro").unwrap();
        group_system.group_create("Protein", "@protein").unwrap();
        let mut full_system = System::from_file("test_files/example.gro").unwrap();
        let original_system = full_system.clone();

        for (group_frame, full_frame) in group_system
            .group_xtc_iter("test_files/short_trajectory.xtc", "Protein")
            .unwrap()
            .zip(
                full_system
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap(),
            )
        {
            let g = group_frame.unwrap();
            let f = full_frame.unwrap();

            compare_box(g.get_box().unwrap(), f.get_box().unwrap());

            for ((a1, a2), a3) in g
                .atoms_iter()
                .zip(f.atoms_iter())
                .zip(original_system.atoms_iter())
            {
                // only atoms of `Protein` should be modified
                if g.group_isin("Protein", a1.get_index()).unwrap() {
                    compare_atoms(a1, a2);
                } else {
                    compare_atoms(a1, a3);
                }
            }
        }
    }

    #[test]
    fn group_xtc_iter_single_particle() {
        let mut group_system = System::from_file("test_files/example.gro").unwrap();
        group_system.group_create("Particle", "serial 1").unwrap();
        let mut full_system = System::from_file("test_files/example.gro").unwrap();
        let original_system = full_system.clone();

        for (group_frame, full_frame) in group_system
            .group_xtc_iter("test_files/short_trajectory.xtc", "Particle")
            .unwrap()
            .zip(
                full_system
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap(),
            )
        {
            let g = group_frame.unwrap();
            let f = full_frame.unwrap();

            compare_box(g.get_box().unwrap(), f.get_box().unwrap());

            for ((a1, a2), a3) in g
                .atoms_iter()
                .zip(f.atoms_iter())
                .zip(original_system.atoms_iter())
            {
                // only atoms of `Particle` should be modified
                if g.group_isin("Particle", a1.get_index()).unwrap() {
                    compare_atoms(a1, a2);
                } else {
                    compare_atoms(a1, a3);
                }
            }
        }
    }

    #[test]
    fn group_xtc_iter_no_particle() {
        let mut group_system = System::from_file("test_files/example.gro").unwrap();
        group_system.group_create("Nothing", "not all").unwrap();
        let mut full_system = System::from_file("test_files/example.gro").unwrap();
        let original_system = full_system.clone();

        for (group_frame, full_frame) in group_system
            .group_xtc_iter("test_files/short_trajectory.xtc", "Nothing")
            .unwrap()
            .zip(
                full_system
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap(),
            )
        {
            let g = group_frame.unwrap();
            let f = full_frame.unwrap();

            compare_box(g.get_box().unwrap(), f.get_box().unwrap());

            // no atom should be modified
            for (a1, a3) in g.atoms_iter().zip(original_system.atoms_iter()) {
                compare_atoms(a1, a3);
            }
        }
    }

    #[test]
    fn group_xtc_iter_range() {
        let mut group_system = System::from_file("test_files/example.gro").unwrap();
        group_system
            .group_create("Membrane", "resname POPC")
            .unwrap();
        let mut full_system = System::from_file("test_files/example.gro").unwrap();
        let original_system = full_system.clone();

        for (group_frame, full_frame) in group_system
            .group_xtc_iter("test_files/short_trajectory.xtc", "Membrane")
            .unwrap()
            .with_range(300.0, 800.0)
            .unwrap()
            .zip(
                full_system
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(300.0, 800.0)
                    .unwrap(),
            )
        {
            let g = group_frame.unwrap();
            let f = full_frame.unwrap();

            compare_box(g.get_box().unwrap(), f.get_box().unwrap());

            for ((a1, a2), a3) in g
                .atoms_iter()
                .zip(f.atoms_iter())
                .zip(original_system.atoms_iter())
            {
                // only atoms of `Membrane` should be modified
                if g.group_isin("Membrane", a1.get_index()).unwrap() {
                    compare_atoms(a1, a2);
                } else {
                    compare_atoms(a1, a3);
                }
            }
        }
    }

    #[test]
    fn group_xtc_range_start_not_found() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.group_create("Membrane", "resname POPC").unwrap();

        match system
            .group_xtc_iter("test_files/short_trajectory.xtc", "Membrane")
            .unwrap()
            .with_range(12000.0, 20000.0)
        {
            Ok(_) => panic!("Iterator should not have been constructed."),
            Err(ReadTrajError::StartNotFound(_)) => (),
            Err(e) => panic!("Incorrect error type {} returned", e),
        }
    }

    #[test]
    fn group_xtc_step_3() {
        let mut group_system = System::from_file("test_files/example.gro").unwrap();
        group_system
            .group_create("Membrane", "resname POPC")
            .unwrap();
        let mut full_system = System::from_file("test_files/example.gro").unwrap();
        let original_system = full_system.clone();

        for (group_frame, full_frame) in group_system
            .group_xtc_iter("test_files/short_trajectory.xtc", "Membrane")
            .unwrap()
            .with_step(3)
            .unwrap()
            .zip(
                full_system
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_step(3)
                    .unwrap(),
            )
        {
            let g = group_frame.unwrap();
            let f = full_frame.unwrap();

            compare_box(g.get_box().unwrap(), f.get_box().unwrap());

            for ((a1, a2), a3) in g
                .atoms_iter()
                .zip(f.atoms_iter())
                .zip(original_system.atoms_iter())
            {
                // only atoms of `Membrane` should be modified
                if g.group_isin("Membrane", a1.get_index()).unwrap() {
                    compare_atoms(a1, a2);
                } else {
                    compare_atoms(a1, a3);
                }
            }
        }
    }

    #[test]
    fn group_xtc_iter_range_and_step() {
        let mut group_system = System::from_file("test_files/example.gro").unwrap();
        group_system
            .group_create("Membrane", "resname POPC")
            .unwrap();
        let mut full_system = System::from_file("test_files/example.gro").unwrap();
        let original_system = full_system.clone();

        for (group_frame, full_frame) in group_system
            .group_xtc_iter("test_files/short_trajectory.xtc", "Membrane")
            .unwrap()
            .with_range(300.0, 800.0)
            .unwrap()
            .with_step(2)
            .unwrap()
            .zip(
                full_system
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(300.0, 800.0)
                    .unwrap()
                    .with_step(2)
                    .unwrap(),
            )
        {
            let g = group_frame.unwrap();
            let f = full_frame.unwrap();

            compare_box(g.get_box().unwrap(), f.get_box().unwrap());

            for ((a1, a2), a3) in g
                .atoms_iter()
                .zip(f.atoms_iter())
                .zip(original_system.atoms_iter())
            {
                // only atoms of `Membrane` should be modified
                if g.group_isin("Membrane", a1.get_index()).unwrap() {
                    compare_atoms(a1, a2);
                } else {
                    compare_atoms(a1, a3);
                }
            }
        }
    }

    #[test]
    fn group_xtc_nonexistent_file() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.group_create("Membrane", "resname POPC").unwrap();

        match system.group_xtc_iter("test_files/nonexistent.xtc", "Membrane") {
            Err(ReadTrajError::FileNotFound(_)) => (),
            _ => panic!("XTC file should not exist."),
        }
    }

    #[test]
    fn group_xtc_nonexistent_group() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system.group_xtc_iter("test_files/short_trajectory.xtc", "Membrane") {
            Err(ReadTrajError::GroupNotFound(x)) => assert_eq!(x, "Membrane"),
            _ => panic!("XTC file should not exist."),
        }
    }
}
