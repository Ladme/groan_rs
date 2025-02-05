// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Integration with the trajectory readers provided by `chemfiles`.
//!
//! See [`ChemfilesReader`] for more information.

use std::{marker::PhantomData, path::Path};

use crate::{
    errors::ReadTrajError,
    files::FileType,
    prelude::{
        FrameData, FrameDataTime, SimBox, TrajFile, TrajFullReadOpen, TrajRangeRead, TrajRead,
        TrajReadOpen, TrajStepRead, TrajStepTimeRead, Vector3D,
    },
    system::System,
};

/// Used when jumping to the start of iteration.
const TIME_PRECISION: f32 = 0.001;
/// Used when rounding box dimensions.
const SPATIAL_PRECISION: f32 = 1e-8;

/// Iterator over a trajectory file using `chemfiles`. Constructed using [`System::traj_iter<ChemfilesReader>`].
///
/// All trajectory formats [supported by `chemfiles`](https://chemfiles.org/chemfiles/latest/formats.html) can be read.
/// However, only XTC, TRR, TNG, GRO, PDB, DCD, Amber NetCDF (NC), and LAMMPSTRJ have been tested inside the `groan_rs` library.
///
/// - For reading XTC and TRR files, it is better to use [`XtcReader`](crate::prelude::XtcReader) and [`TrrReader`](crate::prelude::TrrReader),
///   respectively, as these are more efficient, especially when using the [`with_range`](crate::prelude::TrajReader::with_range)
///   or [`with_step`](crate::prelude::TrajReader::with_step) methods.
/// - For reading GRO trajectories, it is better to use [`GroReader`](crate::prelude::GroReader), since it is also much faster and
///   it can read information about the simulation time (if available).
///
/// ## Limitations
/// - **DCD**: Simulation step information is not available. Instead, the frame number is returned. Simulation time is **assumed** to be in ps.
/// - **Amber NetCDF**: Simulation step information is not available. Instead, the frame number is returned. Simulation time is currently not read.
/// - **GRO**: Simulation step and simulation time information is not available. Unlike `GroReader`, velocities are set to 0, not `None`, when not present.
/// - **PDB**: Simulation step and simulation time information is not available.
/// - **TNG**: Due to the bugginess of the `chemfiles` library, TNG trajectories can only be read when their first frame corresponds to the simulation step of 0.
#[derive(Debug)]
pub struct ChemfilesReader<'a> {
    system: *mut System,
    trajectory: ChemfilesTrajectory,
    phantom: PhantomData<&'a mut System>,
}

/// Wrapper around a chemfiles trajectory.
#[derive(Debug)]
pub struct ChemfilesTrajectory {
    filename: Box<Path>,
    filetype: FileType,
    traj: chemfiles::Trajectory,
    current_frame: usize,
    max_frames: usize,
}

/// Wrapper around a chemfiles frame.
#[derive(Debug)]
pub struct ChemfilesFrame {
    frame: chemfiles::Frame,
    // set either from the current frame or from the previous frame
    current_time: f32,
    // set either from the current frame or from the previous frame
    precision: u64,
}

impl ChemfilesFrame {
    /// Create a new chemfiles frame.
    #[inline(always)]
    fn new(n_atoms: usize, time: f32, precision: u64) -> Self {
        let mut frame = chemfiles::Frame::new();
        frame.resize(n_atoms);

        Self {
            frame,
            current_time: time,
            precision,
        }
    }

    /// Set new simulation time based on information in the frame.
    /// Keep the previous time, if this information is not present.
    #[inline(always)]
    fn set_time(&mut self, file_type: FileType) {
        match self.frame.get("time") {
            Some(chemfiles::Property::Double(time)) => {
                if file_type == FileType::LAMMPSTRJ {
                    let units = match self.frame.get("lammps_units").expect("FATAL GROAN ERROR | ChemfilesFrame::set_time | `lammps_units` property not found.") {
                        chemfiles::Property::String(x) => x,
                        _ => panic!("FATAL GROAN ERROR | ChemfilesFrames::set_time | `lammps_units` should be a string."),
                    };

                    self.current_time = lammps_time_convert(time, &units);
                } else {
                    self.current_time = time as f32;
                }
            }

            _ => (),
        }
    }

    /// Set new precision based on information in the frame.
    /// Keep the previous precision, if this information is not present.
    #[inline(always)]
    fn set_precision(&mut self) {
        match self.frame.get("xtc_precision") {
            Some(chemfiles::Property::Double(prec)) => self.precision = prec as u64,
            _ => (),
        }
    }
}

impl ChemfilesTrajectory {
    #[inline]
    fn new(filename: impl AsRef<Path>, system: &System) -> Result<Self, ReadTrajError> {
        // disable intrusive warnings by the chemfiles library
        chemfiles::set_warning_callback(|_| {});

        let mut traj =
            chemfiles::Trajectory::open(&filename, 'r').map_err(ReadTrajError::ChemfilesError)?;

        // read the first frame to check whether this is a valid trajectory file
        let mut frame = chemfiles::Frame::new();
        traj.read(&mut frame)
            .map_err(ReadTrajError::ChemfilesError)?;

        // check that the number of atoms matches the number of atoms in the System
        if frame.positions().len() != system.get_n_atoms() {
            return Err(ReadTrajError::AtomsNumberMismatch(Box::from(
                filename.as_ref(),
            )));
        }

        // open the trajectory again to return to the start
        let mut traj =
            chemfiles::Trajectory::open(&filename, 'r').map_err(ReadTrajError::ChemfilesError)?;

        // we have to get the number of frames in the trajectory to be able to stop at the appropriate frame
        let max_frames = traj.nsteps();
        Ok(Self {
            filename: Box::from(filename.as_ref()),
            filetype: FileType::from_name(
            filename
                    .as_ref()
                    .to_str()
                    .expect("FATAL GROAN ERROR | ChemfilesTrajectory::new | Could not convert filename to string.")),
            traj,
            max_frames,
            current_frame: 0,
        })
    }
}

impl TrajFile for ChemfilesTrajectory {}

impl FrameData for ChemfilesFrame {
    type TrajFile = ChemfilesTrajectory;

    #[inline]
    fn from_frame(
        traj_file: &mut Self::TrajFile,
        system: &System,
    ) -> Option<Result<Self, ReadTrajError>>
    where
        Self: Sized,
    {
        let mut frame = ChemfilesFrame::new(
            system.get_n_atoms(),
            system.get_simulation_time(),
            system.get_precision(),
        );

        if traj_file.current_frame >= traj_file.max_frames {
            return None;
        }

        match traj_file.traj.read(&mut frame.frame) {
            Ok(_) => {
                traj_file.current_frame += 1;
                frame.set_time(traj_file.filetype);
                frame.set_precision();
                Some(Ok(frame))
            }
            Err(e) => Some(Err(ReadTrajError::ChemfilesError(e))),
        }
    }

    fn update_system(self, system: &mut System) {
        // trr file frames do not have to contain positions
        let positions = match self.frame.get("has_positions") {
            None | Some(chemfiles::Property::Bool(true)) => Some(self.frame.positions()),
            Some(chemfiles::Property::Bool(false)) => None,
            _ => panic!("FATAL GROAN ERROR | ChemfilesFrame::update_system | `has_positions` must be boolean."),
        };

        let velocities = self.frame.velocities();
        let atoms = system.atoms_iter_mut();

        match (positions, velocities) {
            (Some(positions), Some(velocities)) => {
                for (atom, (pos, vel)) in atoms.zip(positions.iter().zip(velocities.iter())) {
                    unsafe {
                        atom.set_position(vector_from_slice(pos));
                        atom.set_velocity(vector_from_slice(vel));
                    }
                    atom.reset_force();
                }
            }

            (Some(positions), None) => {
                for (atom, pos) in atoms.zip(positions.iter()) {
                    unsafe {
                        atom.set_position(vector_from_slice(pos));
                    }
                    atom.reset_velocity();
                    atom.reset_force();
                }
            }

            (None, Some(velocities)) => {
                for (atom, vel) in atoms.zip(velocities.iter()) {
                    unsafe {
                        atom.set_velocity(vector_from_slice(vel));
                    }
                    atom.reset_position();
                    atom.reset_force();
                }
            }

            (None, None) => {
                for atom in atoms {
                    atom.reset_position();
                    atom.reset_velocity();
                    atom.reset_force();
                }
            }
        }

        system.set_simulation_step(self.frame.step() as u64);
        system.set_simulation_time(self.current_time);
        system.set_precision(self.precision);
        let cell = self.frame.cell().matrix();
        let simbox = SimBox::from([
            round_box_dim((cell[0][0] / 10.0) as f32),
            round_box_dim((cell[1][1] / 10.0) as f32),
            round_box_dim((cell[2][2] / 10.0) as f32),
            round_box_dim((cell[1][0] / 10.0) as f32),
            round_box_dim((cell[2][0] / 10.0) as f32),
            round_box_dim((cell[0][1] / 10.0) as f32),
            round_box_dim((cell[2][1] / 10.0) as f32),
            round_box_dim((cell[0][2] / 10.0) as f32),
            round_box_dim((cell[1][2] / 10.0) as f32),
        ]);
        system.set_box(simbox);
    }
}

/// Round box dimension. This is necessary to avoid some floating point errors.
#[inline(always)]
fn round_box_dim(value: f32) -> f32 {
    if value.abs() < SPATIAL_PRECISION {
        0.0
    } else {
        value
    }
}

// Helper to convert a slice to a Vector3D.
#[inline(always)]
unsafe fn vector_from_slice(slice: &[f64; 3]) -> Vector3D {
    Vector3D::new(
        (*slice.get_unchecked(0) / 10.0) as f32,
        (*slice.get_unchecked(1) / 10.0) as f32,
        (*slice.get_unchecked(2) / 10.0) as f32,
    )
}

impl<'a> TrajRead<'a> for ChemfilesReader<'a> {
    type FrameData = ChemfilesFrame;

    #[inline(always)]
    fn get_system(&mut self) -> *mut System {
        self.system
    }

    #[inline(always)]
    fn get_file_handle(&mut self) -> &mut ChemfilesTrajectory {
        &mut self.trajectory
    }
}

impl<'a> TrajReadOpen<'a> for ChemfilesReader<'a> {
    #[inline(always)]
    fn initialize(
        system: &'a mut System,
        filename: impl AsRef<std::path::Path>,
        group: Option<&str>,
    ) -> Result<Self, ReadTrajError>
    where
        Self: Sized,
    {
        match group {
            None => Self::new(system, filename),
            Some(_) => panic!("FATAL GROAN ERROR | ChemfilesReader::initialize | Chemfiles do not support partial-frame reading."),
        }
    }
}

impl<'a> TrajFullReadOpen<'a> for ChemfilesReader<'a> {
    #[inline]
    fn new(system: &'a mut System, filename: impl AsRef<Path>) -> Result<Self, ReadTrajError>
    where
        Self: Sized,
    {
        Ok(ChemfilesReader {
            system: system as *mut System,
            trajectory: ChemfilesTrajectory::new(filename, system)?,
            phantom: PhantomData,
        })
    }
}

impl FrameDataTime for ChemfilesFrame {
    #[inline(always)]
    fn get_time(&self) -> f32 {
        self.current_time
    }
}

impl<'a> TrajRangeRead<'a> for ChemfilesReader<'a> {
    fn jump_to_start(&mut self, start_time: f32) -> Result<(), ReadTrajError> {
        // optimization is possible by assuming frames are saved at equal time intervals
        // we could then calculate the target frame index in O(1) using the time difference (dt) between the first two frames
        // then directly access the expected starting frame and only fall back to iteration if the estimate overshoots
        // however, if there are multiple frames with the same time, we may start from the wrong one
        // we must always start from the first frame with the corresponding time, not from any frame with that time
        // therefore, such heuristics could lead to errors (although only in rare cases) which is unacceptable
        // it might be worth it to add this optimization as an opt-in feature though
        let mut frame_index = 0;
        loop {
            let system = unsafe { &*self.get_system() };
            match ChemfilesFrame::from_frame(&mut self.trajectory, system) {
                None => return Err(ReadTrajError::StartNotFound(start_time.to_string())),
                Some(Err(e)) => return Err(e),
                Some(Ok(frame)) => {
                    if frame.current_time >= start_time - TIME_PRECISION {
                        if frame_index > 0 {
                            // read again the previous frame; this sets the file pointer to the correct position
                            let mut frame = chemfiles::Frame::new();
                            frame.resize(system.get_n_atoms());
                            self.trajectory.traj
                                .read_step(frame_index - 1, &mut frame)
                                .expect("FATAL GROAN ERROR | ChemfilesReader::jump_to_start | Could not re-read a trajectory frame.");

                            // chemfiles is buggy and does not properly propagate the file pointer when using `read_step` with LAMMPSTRJ and DCD,
                            // so we have to move it manually
                            if self.trajectory.filetype == FileType::LAMMPSTRJ
                                || self.trajectory.filetype == FileType::DCD
                            {
                                self.trajectory.traj
                                    .read(&mut frame)
                                    .expect("FATAL GROAN ERROR | ChemfilesReader::jump_to_start | Could not re-read trajectory frame (LAMMPSTRJ/DCD).");
                            }
                        } else {
                            // or if this is the first frame, open the trajectory again
                            self.trajectory.traj = chemfiles::Trajectory::open(&self.trajectory.filename, 'r')
                                .expect("FATAL GROAN ERROR | ChemfilesRead::jump_to_start | Could not re-open a trajectory file.");
                        }

                        // set the current frame index
                        self.trajectory.current_frame = frame_index;

                        return Ok(());
                    }
                }
            }

            frame_index += 1;
        }
    }
}

impl<'a> TrajStepRead<'a> for ChemfilesReader<'a> {
    #[inline(always)]
    fn skip_frame(&mut self) -> Result<bool, ReadTrajError> {
        // we could skip frames more efficiently using `read_step`
        // but that would require opting out of or redoing the trajectory reader trait system
        let system = unsafe { &*self.get_system() };
        match ChemfilesFrame::from_frame(&mut self.trajectory, system) {
            None => Ok(false),
            Some(Err(e)) => Err(e),
            Some(Ok(_)) => Ok(true),
        }
    }
}

impl<'a> TrajStepTimeRead<'a> for ChemfilesReader<'a> {
    #[inline(always)]
    fn skip_frame_time(&mut self) -> Result<Option<f32>, ReadTrajError> {
        let system = unsafe { &*self.get_system() };
        match ChemfilesFrame::from_frame(&mut self.trajectory, system) {
            None => Ok(None),
            Some(Err(e)) => Err(e),
            Some(Ok(frame)) => Ok(Some(frame.current_time)),
        }
    }
}

/// Converts time in the specified LAMMPS units to picoseconds.
fn lammps_time_convert(time: f64, units: &str) -> f32 {
    let converter = match units {
        "lj" => 1.0,        // dimension-less, assume ps
        "real" => 1e-3,     // fs
        "metal" => 1.0,     // ps
        "si" => 1e12,       // s
        "cgs" => 1e12,      // s
        "electron" => 1e-3, // fs
        "micro" => 1e6,     // μs
        "nano" => 1e3,      // ns,
        x => panic!(
            "FATAL GROAN ERROR | chemfiles::lammps_time_convert | Unknown LAMMPS units `{}`.",
            x
        ),
    };

    (time * converter) as f32
}

#[cfg(any(feature = "molly", not(feature = "no-xdrfile")))]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utilities::utilities::{compare_atoms_without_forces, compare_box};
    use float_cmp::assert_approx_eq;

    /// Compare two generic trajectory iterators.
    fn compare_iterators<'a>(
        mut iter1: impl Iterator<Item = Result<&'a mut System, ReadTrajError>>,
        mut iter2: impl Iterator<Item = Result<&'a mut System, ReadTrajError>>,
    ) {
        while let (Some(frame1), Some(frame2)) = (iter1.next(), iter2.next()) {
            let frame1 = frame1.unwrap();
            let frame2 = frame2.unwrap();

            /*println!(
                "{} {}",
                frame1.get_simulation_time(),
                frame2.get_simulation_time()
            );*/

            compare_box(frame1.get_box().unwrap(), frame2.get_box().unwrap());
            assert_eq!(frame1.get_simulation_step(), frame2.get_simulation_step());
            assert_approx_eq!(
                f32,
                frame1.get_simulation_time(),
                frame2.get_simulation_time()
            );

            for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                // chemfiles does not load forces even if they are available...
                compare_atoms_without_forces(atom1, atom2);
            }
        }

        // check that both operators are exhausted
        assert!(iter1.next().is_none() && iter2.next().is_none());
    }

    #[cfg(any(feature = "molly", not(feature = "no-xdrfile")))]
    mod tests_xtc {
        use crate::test_utilities::utilities::compare_atoms_without_forces;

        use super::*;

        /// Compare two xtc iterators.
        fn compare_xtc_iterators<'a>(
            mut iter1: impl Iterator<Item = Result<&'a mut System, ReadTrajError>>,
            mut iter2: impl Iterator<Item = Result<&'a mut System, ReadTrajError>>,
        ) {
            while let (Some(frame1), Some(frame2)) = (iter1.next(), iter2.next()) {
                let frame1 = frame1.unwrap();
                let frame2 = frame2.unwrap();

                /*println!(
                    "{} {}",
                    frame1.get_simulation_time(),
                    frame2.get_simulation_time()
                );*/

                compare_box(frame1.get_box().unwrap(), frame2.get_box().unwrap());
                assert_eq!(frame1.get_precision(), frame2.get_precision());
                assert_eq!(frame1.get_simulation_step(), frame2.get_simulation_step());
                assert_approx_eq!(
                    f32,
                    frame1.get_simulation_time(),
                    frame2.get_simulation_time()
                );

                for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                    compare_atoms_without_forces(atom1, atom2);
                }
            }

            // check that both operators are exhausted
            assert!(iter1.next().is_none() && iter2.next().is_none());
        }

        #[test]
        fn read_xtc_isolated() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            let n_frames = system
                .traj_iter::<ChemfilesReader>("test_files/short_trajectory.xtc")
                .unwrap()
                .count();

            assert_eq!(n_frames, 11);
        }

        #[test]
        fn read_xtc_pass() {
            for (gro_file, xtc_file) in [
                "test_files/example.gro",
                "test_files/triclinic.gro",
                "test_files/octahedron.gro",
                "test_files/dodecahedron.gro",
            ]
            .into_iter()
            .zip(
                [
                    "test_files/short_trajectory.xtc",
                    "test_files/triclinic_trajectory.xtc",
                    "test_files/octahedron_trajectory.xtc",
                    "test_files/dodecahedron_trajectory.xtc",
                ]
                .into_iter(),
            ) {
                let mut system_chem = System::from_file(gro_file).unwrap();
                let mut system_xtc = system_chem.clone();

                let xtc_iter = system_xtc.xtc_iter(xtc_file).unwrap();
                let chem_iter = system_chem.traj_iter::<ChemfilesReader>(xtc_file).unwrap();

                compare_xtc_iterators(xtc_iter, chem_iter);
            }
        }

        #[test]
        fn read_xtc_unmatching() {
            let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/short_trajectory.xtc") {
                Ok(_) => panic!("XTC file should not be valid."),
                Err(ReadTrajError::AtomsNumberMismatch(_)) => (),
                Err(e) => panic!("Unexpected error type `{}` returned.", e),
            }
        }

        #[test]
        fn read_xtc_nonexistent() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/nonexistent.xtc") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                _ => panic!("XTC file should not exist."),
            }
        }

        #[test]
        fn read_xtc_not_xtc() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/fake_xtc.xtc") {
                Ok(_) => panic!("XTC file should not be valid."),
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error type `{}` returned.", e),
            }
        }

        #[test]
        fn read_xtc_ranges() {
            for (start, end) in [
                (0.0, 100_000.0),
                (200.0, 600.0),
                (300.0, 500.0),
                (500.0, 500.0),
                (300.0, 100_000.0),
            ] {
                let mut system_chem = System::from_file("test_files/example.gro").unwrap();
                let mut system_xtc = system_chem.clone();

                let xtc_iter = system_xtc
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(start, end)
                    .unwrap();

                let chem_iter = system_chem
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(start, end)
                    .unwrap();

                compare_xtc_iterators(xtc_iter, chem_iter);
            }
        }

        #[test]
        fn read_xtc_steps() {
            for step in [1, 2, 3, 5, 23] {
                let mut system_chem = System::from_file("test_files/example.gro").unwrap();
                let mut system_xtc = system_chem.clone();

                let xtc_iter = system_xtc
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                let chem_iter = system_chem
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                compare_xtc_iterators(xtc_iter, chem_iter);
            }
        }

        #[test]
        fn read_xtc_ranges_steps() {
            for (start, end, step) in [(0.0, 100_000.0, 1), (300.0, 800.0, 2), (100.0, 900.0, 4)] {
                let mut system_chem = System::from_file("test_files/example.gro").unwrap();
                let mut system_xtc = system_chem.clone();

                let xtc_iter = system_xtc
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(start, end)
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                let chem_iter = system_chem
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_step(step)
                    .unwrap()
                    .with_range(start, end)
                    .unwrap();

                compare_xtc_iterators(xtc_iter, chem_iter);
            }
        }

        #[test]
        fn read_xtc_fail_start_not_found() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system
                .traj_iter::<ChemfilesReader>("test_files/short_trajectory.xtc")
                .unwrap()
                .with_range(1100.0, 2000.0)
            {
                Ok(_) => panic!("Function should have failed."),
                Err(ReadTrajError::StartNotFound(_)) => (),
                Err(e) => panic!("Unexpected error type `{}` returned.", e),
            }
        }
    }

    #[cfg(not(feature = "no-xdrfile"))]
    mod tests_trr {
        use super::*;

        #[test]
        fn read_trr_isolated() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            let n_frames = system
                .traj_iter::<ChemfilesReader>("test_files/short_trajectory.trr")
                .unwrap()
                .count();

            assert_eq!(n_frames, 9);
        }

        #[test]
        fn read_trr_pass() {
            for (gro_file, trr_file) in [
                "test_files/example.gro",
                "test_files/example.gro",
                "test_files/triclinic.gro",
                "test_files/octahedron.gro",
                "test_files/dodecahedron.gro",
            ]
            .into_iter()
            .zip(
                [
                    "test_files/short_trajectory.trr",
                    "test_files/short_trajectory_double.trr",
                    "test_files/triclinic_trajectory.trr",
                    "test_files/octahedron_trajectory.trr",
                    "test_files/dodecahedron_trajectory.trr",
                ]
                .into_iter(),
            ) {
                let mut system_chem = System::from_file(gro_file).unwrap();
                let mut system_trr = system_chem.clone();

                let trr_iter = system_trr.trr_iter(trr_file).unwrap();
                let chem_iter = system_chem.traj_iter::<ChemfilesReader>(trr_file).unwrap();

                compare_iterators(trr_iter, chem_iter);
            }
        }

        #[test]
        fn read_trr_unmatching() {
            let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/short_trajectory.trr") {
                Err(ReadTrajError::AtomsNumberMismatch(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("TRR file should not be valid."),
            }
        }

        #[test]
        fn read_trr_nonexistent() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/nonexistent.trr") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("TRR file should not exist."),
            }
        }

        #[test]
        fn read_trr_not_trr() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/fake_trr.trr") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("File should not be a trr file."),
            }
        }

        #[test]
        fn read_trr_ranges() {
            for (start, end) in [
                (0.0, 1000.0),
                (130.0, 400.0),
                (200.0, 600.0),
                (480.0, 480.0),
                (200.0, 1000.0),
            ] {
                let mut system_chem = System::from_file("test_files/example.gro").unwrap();
                let mut system_xtc = system_chem.clone();

                let trr_iter = system_xtc
                    .trr_iter("test_files/short_trajectory.trr")
                    .unwrap()
                    .with_range(start, end)
                    .unwrap();

                let chem_iter = system_chem
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.trr")
                    .unwrap()
                    .with_range(start, end)
                    .unwrap();

                compare_iterators(trr_iter, chem_iter);
            }
        }

        #[test]
        fn read_trr_steps() {
            for step in [1, 2, 3, 5, 23] {
                let mut system_chem = System::from_file("test_files/example.gro").unwrap();
                let mut system_xtc = system_chem.clone();

                let trr_iter = system_xtc
                    .trr_iter("test_files/short_trajectory.trr")
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                let chem_iter = system_chem
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.trr")
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                compare_iterators(trr_iter, chem_iter);
            }
        }

        #[test]
        fn read_trr_ranges_steps() {
            for (start, end, step) in [(0.0, 1000.0, 1), (200.0, 600.0, 2), (100.0, 500.0, 3)] {
                let mut system_chem = System::from_file("test_files/example.gro").unwrap();
                let mut system_xtc = system_chem.clone();

                let trr_iter = system_xtc
                    .trr_iter("test_files/short_trajectory.trr")
                    .unwrap()
                    .with_range(start, end)
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                let chem_iter = system_chem
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.trr")
                    .unwrap()
                    .with_step(step)
                    .unwrap()
                    .with_range(start, end)
                    .unwrap();

                compare_iterators(trr_iter, chem_iter);
            }
        }

        #[test]
        fn read_trr_fail_start_not_found() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system
                .traj_iter::<ChemfilesReader>("test_files/short_trajectory.trr")
                .unwrap()
                .with_range(720.0, 1000.0)
            {
                Ok(_) => panic!("Function should have failed."),
                Err(ReadTrajError::StartNotFound(_)) => (),
                Err(e) => panic!("Unexpected error type `{}` returned.", e),
            }
        }
    }

    #[cfg(any(feature = "molly", not(feature = "no-xdrfile")))]
    mod tests_tng {
        use super::*;

        #[test]
        fn read_tng_isolated() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            let n_frames = system
                .traj_iter::<ChemfilesReader>("test_files/short_trajectory.tng")
                .unwrap()
                .count();

            assert_eq!(n_frames, 11);
        }

        #[test]
        fn read_tng_pass() {
            for (gro_file, xtc_file, tng_file) in [
                (
                    "test_files/example.gro",
                    "test_files/short_trajectory.xtc",
                    "test_files/short_trajectory.tng",
                ),
                (
                    "test_files/octahedron.gro",
                    "test_files/octahedron_trajectory.xtc",
                    "test_files/octahedron_trajectory.tng",
                ),
            ] {
                let mut system_tng = System::from_file(gro_file).unwrap();
                let mut system_xtc = system_tng.clone();

                let xtc_iter = system_xtc.xtc_iter(xtc_file).unwrap();
                let tng_iter = system_tng.traj_iter::<ChemfilesReader>(tng_file).unwrap();

                compare_iterators(xtc_iter, tng_iter);
            }
        }

        #[test]
        fn read_tng_unmatching() {
            let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/short_trajectory.tng") {
                Err(ReadTrajError::AtomsNumberMismatch(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("TNG file should not be valid."),
            }
        }

        #[test]
        fn read_tng_nonexistent() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/nonexistent.tng") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("TNG file should not exist."),
            }
        }

        #[test]
        fn read_tng_not_tng() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/fake_tng.tng") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("File should not be a tng file."),
            }
        }

        #[test]
        fn read_tng_ranges() {
            for (start, end) in [
                (0.0, 100_000.0),
                (200.0, 600.0),
                (300.0, 500.0),
                (500.0, 500.0),
                (300.0, 100_000.0),
            ] {
                let mut system_tng = System::from_file("test_files/example.gro").unwrap();
                let mut system_xtc = system_tng.clone();

                let xtc_iter = system_xtc
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(start, end)
                    .unwrap();

                let tng_iter = system_tng
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.tng")
                    .unwrap()
                    .with_range(start, end)
                    .unwrap();

                compare_iterators(xtc_iter, tng_iter);
            }
        }

        #[test]
        fn read_tng_steps() {
            for step in [1, 2, 3, 5, 23] {
                let mut system_tng = System::from_file("test_files/example.gro").unwrap();
                let mut system_xtc = system_tng.clone();

                let xtc_iter = system_xtc
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                let tng_iter = system_tng
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.tng")
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                compare_iterators(xtc_iter, tng_iter);
            }
        }

        #[test]
        fn read_tng_ranges_steps() {
            for (start, end, step) in [(0.0, 100_000.0, 1), (300.0, 800.0, 2), (100.0, 900.0, 4)] {
                let mut system_tng = System::from_file("test_files/example.gro").unwrap();
                let mut system_xtc = system_tng.clone();

                let xtc_iter = system_xtc
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(start, end)
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                let tng_iter = system_tng
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.tng")
                    .unwrap()
                    .with_step(step)
                    .unwrap()
                    .with_range(start, end)
                    .unwrap();

                compare_iterators(xtc_iter, tng_iter);
            }
        }
    }

    #[cfg(any(feature = "molly", not(feature = "no-xdrfile")))]
    mod tests_dcd {
        use super::*;

        /// Compare two iterators, one of which is the DCD iterator.
        fn compare_dcd_iterators<'a>(
            mut iter1: impl Iterator<Item = Result<&'a mut System, ReadTrajError>>,
            mut iter2: impl Iterator<Item = Result<&'a mut System, ReadTrajError>>,
        ) {
            while let (Some(frame1), Some(frame2)) = (iter1.next(), iter2.next()) {
                let frame1 = frame1.unwrap();
                let frame2 = frame2.unwrap();

                /*println!(
                    "{} {}",
                    frame1.get_simulation_time(),
                    frame2.get_simulation_time()
                );*/

                compare_box(frame1.get_box().unwrap(), frame2.get_box().unwrap());
                assert_approx_eq!(
                    f32,
                    frame1.get_simulation_time(),
                    frame2.get_simulation_time()
                );

                for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                    // chemfiles does not load forces even if they are available...
                    compare_atoms_without_forces(atom1, atom2);
                }
            }

            // check that both operators are exhausted
            assert!(iter1.next().is_none() && iter2.next().is_none());
        }

        #[test]
        fn read_dcd_isolated() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            let n_frames = system
                .traj_iter::<ChemfilesReader>("test_files/short_trajectory.dcd")
                .unwrap()
                .count();

            assert_eq!(n_frames, 11);
        }

        #[test]
        fn read_dcd_pass() {
            for (gro_file, xtc_file, dcd_file) in [
                (
                    "test_files/example.gro",
                    "test_files/short_trajectory.xtc",
                    "test_files/short_trajectory.dcd",
                ),
                (
                    "test_files/octahedron.gro",
                    "test_files/octahedron_trajectory.xtc",
                    "test_files/octahedron_trajectory.dcd",
                ),
            ] {
                let mut system_dcd = System::from_file(gro_file).unwrap();
                let mut system_xtc = system_dcd.clone();

                let xtc_iter = system_xtc.xtc_iter(xtc_file).unwrap();
                let dcd_iter = system_dcd.traj_iter::<ChemfilesReader>(dcd_file).unwrap();

                compare_dcd_iterators(xtc_iter, dcd_iter);
            }
        }

        #[test]
        fn read_dcd_unmatching() {
            let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/short_trajectory.dcd") {
                Err(ReadTrajError::AtomsNumberMismatch(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("DCD file should not be valid."),
            }
        }

        #[test]
        fn read_dcd_nonexistent() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/nonexistent.dcd") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("DCD file should not exist."),
            }
        }

        #[test]
        fn read_dcd_not_dcd() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/fake_dcd.dcd") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("File should not be a DCD file."),
            }
        }

        #[test]
        fn read_dcd_ranges() {
            for (start, end) in [
                (0.0, 100_000.0),
                (200.0, 600.0),
                (300.0, 500.0),
                (500.0, 500.0),
                (300.0, 100_000.0),
            ] {
                let mut system_dcd = System::from_file("test_files/example.gro").unwrap();
                let mut system_xtc = system_dcd.clone();

                let xtc_iter = system_xtc
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(start, end)
                    .unwrap();

                let dcd_iter = system_dcd
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.dcd")
                    .unwrap()
                    .with_range(start, end)
                    .unwrap();

                compare_dcd_iterators(xtc_iter, dcd_iter);
            }
        }

        #[test]
        fn read_dcd_steps() {
            for step in [1, 2, 3, 5, 23] {
                let mut system_dcd = System::from_file("test_files/example.gro").unwrap();
                let mut system_xtc = system_dcd.clone();

                let xtc_iter = system_xtc
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                let dcd_iter = system_dcd
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.dcd")
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                compare_dcd_iterators(xtc_iter, dcd_iter);
            }
        }

        #[test]
        fn read_dcd_ranges_steps() {
            for (start, end, step) in [(0.0, 100_000.0, 1), (300.0, 800.0, 2), (100.0, 900.0, 4)] {
                let mut system_dcd = System::from_file("test_files/example.gro").unwrap();
                let mut system_xtc = system_dcd.clone();

                let xtc_iter = system_xtc
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(start, end)
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                let dcd_iter = system_dcd
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.dcd")
                    .unwrap()
                    .with_step(step)
                    .unwrap()
                    .with_range(start, end)
                    .unwrap();

                compare_dcd_iterators(xtc_iter, dcd_iter);
            }
        }
    }

    #[cfg(any(feature = "molly", not(feature = "no-xdrfile")))]
    mod tests_nc {
        use super::*;

        /// Compare two trajectory iterators, where one of them is Amber NetCDF iterator.
        fn compare_nc_iterators<'a>(
            mut iter1: impl Iterator<Item = Result<&'a mut System, ReadTrajError>>,
            mut iter2: impl Iterator<Item = Result<&'a mut System, ReadTrajError>>,
        ) {
            while let (Some(frame1), Some(frame2)) = (iter1.next(), iter2.next()) {
                let frame1 = frame1.unwrap();
                let frame2 = frame2.unwrap();

                /*println!(
                    "{} {}",
                    frame1.get_simulation_time(),
                    frame2.get_simulation_time()
                );*/

                compare_box(frame1.get_box().unwrap(), frame2.get_box().unwrap());

                for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                    // chemfiles does not load forces even if they are available...
                    compare_atoms_without_forces(atom1, atom2);
                }
            }

            // check that both operators are exhausted
            assert!(iter1.next().is_none() && iter2.next().is_none());
        }

        #[test]
        fn read_nc_isolated() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            let n_frames = system
                .traj_iter::<ChemfilesReader>("test_files/short_trajectory.nc")
                .unwrap()
                .count();

            assert_eq!(n_frames, 11);
        }

        #[test]
        fn read_nc_pass() {
            for (gro_file, xtc_file, nc_file) in [
                (
                    "test_files/example.gro",
                    "test_files/short_trajectory.xtc",
                    "test_files/short_trajectory.nc",
                ),
                (
                    "test_files/octahedron.gro",
                    "test_files/octahedron_trajectory.xtc",
                    "test_files/octahedron_trajectory.nc",
                ),
            ] {
                let mut system_nc = System::from_file(gro_file).unwrap();
                let mut system_xtc = system_nc.clone();

                let xtc_iter = system_xtc.xtc_iter(xtc_file).unwrap();
                let nc_iter = system_nc.traj_iter::<ChemfilesReader>(nc_file).unwrap();

                compare_nc_iterators(xtc_iter, nc_iter);
            }
        }

        #[test]
        fn read_nc_unmatching() {
            let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/short_trajectory.nc") {
                Err(ReadTrajError::AtomsNumberMismatch(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("NC file should not be valid."),
            }
        }

        #[test]
        fn read_nc_nonexistent() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/nonexistent.nc") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("NC file should not exist."),
            }
        }

        #[test]
        fn read_nc_not_nc() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/fake_nc.nc") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("File should not be an NC file."),
            }
        }
    }

    #[cfg(any(feature = "molly", not(feature = "no-xdrfile")))]
    mod tests_lammps {
        use super::super::lammps_time_convert;
        use super::*;

        #[test]
        fn tests_lammps_convert() {
            let time = 1500.0;
            let units = [
                "lj", "real", "metal", "si", "cgs", "electron", "micro", "nano",
            ];
            let converted = [1500.0, 1.5, 1500.0, 1.5e15, 1.5e15, 1.5, 1.5e9, 1.5e6];

            for (unit, conv) in units.into_iter().zip(converted.into_iter()) {
                assert_approx_eq!(f32, lammps_time_convert(time, unit), conv);
            }
        }

        #[test]
        fn read_lammps_isolated() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            let n_frames = system
                .traj_iter::<ChemfilesReader>("test_files/short_trajectory.lammpstrj")
                .unwrap()
                .count();

            assert_eq!(n_frames, 11);
        }

        #[test]
        fn read_lammps_pass() {
            for (gro_file, xtc_file, lammps_file) in [
                (
                    "test_files/example.gro",
                    "test_files/short_trajectory.xtc",
                    "test_files/short_trajectory.lammpstrj",
                ),
                (
                    "test_files/octahedron.gro",
                    "test_files/octahedron_trajectory.xtc",
                    "test_files/octahedron_trajectory.lammpstrj",
                ),
            ] {
                let mut system_lammps = System::from_file(gro_file).unwrap();
                let mut system_xtc = system_lammps.clone();

                let xtc_iter = system_xtc.xtc_iter(xtc_file).unwrap();
                let lammps_iter = system_lammps
                    .traj_iter::<ChemfilesReader>(lammps_file)
                    .unwrap();

                compare_iterators(xtc_iter, lammps_iter);
            }
        }

        #[test]
        fn read_lammps_unmatching() {
            let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/short_trajectory.lammpstrj") {
                Err(ReadTrajError::AtomsNumberMismatch(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("LAMMPSTRJ file should not be valid."),
            }
        }

        #[test]
        fn read_lammps_nonexistent() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/nonexistent.lammpstrj") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("LAMMPSTRJ file should not exist."),
            }
        }

        #[test]
        fn read_lammps_not_lammps() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/fake_lammps.lammpstrj") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("File should not be a LAMMPSTRJ file."),
            }
        }

        #[test]
        fn read_lammps_ranges() {
            for (start, end) in [
                (0.0, 100_000.0),
                (200.0, 600.0),
                (300.0, 500.0),
                (500.0, 500.0),
                (300.0, 100_000.0),
            ] {
                let mut system_lammps = System::from_file("test_files/example.gro").unwrap();
                let mut system_xtc = system_lammps.clone();

                let xtc_iter = system_xtc
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(start, end)
                    .unwrap();

                let lammps_iter = system_lammps
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.lammpstrj")
                    .unwrap()
                    .with_range(start, end)
                    .unwrap();

                compare_iterators(xtc_iter, lammps_iter);
            }
        }

        #[test]
        fn read_lammps_steps() {
            for step in [1, 2, 3, 5, 23] {
                let mut system_lammps = System::from_file("test_files/example.gro").unwrap();
                let mut system_xtc = system_lammps.clone();

                let xtc_iter = system_xtc
                    .xtc_iter("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                let lammps_iter = system_lammps
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.lammpstrj")
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                compare_iterators(xtc_iter, lammps_iter);
            }
        }
    }

    #[test]
    fn read_lammps_ranges_steps() {
        for (start, end, step) in [(0.0, 100_000.0, 1), (300.0, 800.0, 2), (100.0, 900.0, 4)] {
            let mut system_lammps = System::from_file("test_files/example.gro").unwrap();
            let mut system_xtc = system_lammps.clone();

            let xtc_iter = system_xtc
                .xtc_iter("test_files/short_trajectory.xtc")
                .unwrap()
                .with_range(start, end)
                .unwrap()
                .with_step(step)
                .unwrap();

            let lammps_iter = system_lammps
                .traj_iter::<ChemfilesReader>("test_files/short_trajectory.lammpstrj")
                .unwrap()
                .with_step(step)
                .unwrap()
                .with_range(start, end)
                .unwrap();

            compare_iterators(xtc_iter, lammps_iter);
        }
    }

    mod tests_gro {
        use crate::test_utilities::utilities::compare_atoms_without_forces_and_velocities;

        use super::*;

        /// Compare two iterators at least one of which is a gro iterator.
        fn compare_gro_iterators<'a>(
            mut iter1: impl Iterator<Item = Result<&'a mut System, ReadTrajError>>,
            mut iter2: impl Iterator<Item = Result<&'a mut System, ReadTrajError>>,
        ) {
            while let (Some(frame1), Some(frame2)) = (iter1.next(), iter2.next()) {
                let frame1 = frame1.unwrap();
                let frame2 = frame2.unwrap();

                compare_box(frame1.get_box().unwrap(), frame2.get_box().unwrap());

                for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                    // if velocities are missing from a gro file, chemfiles set them to 0 which is not consistent with
                    // either other chemfiles formats nor the groan_rs convention
                    compare_atoms_without_forces_and_velocities(atom1, atom2);
                }
            }

            // check that both operators are exhausted
            assert!(iter1.next().is_none() && iter2.next().is_none());
        }

        #[test]
        fn read_gro_isolated() {
            let mut system = System::from_file("test_files/protein_trajectory.gro").unwrap();

            let n_frames = system
                .traj_iter::<ChemfilesReader>("test_files/protein_trajectory.gro")
                .unwrap()
                .count();

            assert_eq!(n_frames, 11);
        }

        #[test]
        fn read_gro_pass() {
            let mut system = System::from_file("test_files/protein_trajectory.gro").unwrap();
            let mut system_chemfiles = system.clone();

            let gro_iter = system
                .gro_iter("test_files/protein_trajectory.gro")
                .unwrap();
            let chemfiles_iter = system_chemfiles
                .traj_iter::<ChemfilesReader>("test_files/protein_trajectory.gro")
                .unwrap();

            compare_gro_iterators(gro_iter, chemfiles_iter);
        }

        #[test]
        fn read_gro_unmatching() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/protein_trajectory.gro") {
                Err(ReadTrajError::AtomsNumberMismatch(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("GRO file should not be valid."),
            }
        }

        #[test]
        fn read_gro_nonexistent() {
            let mut system = System::from_file("test_files/protein_trajectory.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/nonexistent.gro") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("GRO file should not exist."),
            }
        }

        #[test]
        fn read_gro_not_gro() {
            let mut system = System::from_file("test_files/protein_trajectory.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/example_empty.gro") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("File should not be a GRO file."),
            }
        }

        #[test]
        fn read_gro_steps() {
            for step in [1, 2, 3, 5, 23] {
                let mut system = System::from_file("test_files/protein_trajectory.gro").unwrap();
                let mut system_chemfiles = system.clone();

                let gro_iter = system
                    .gro_iter("test_files/protein_trajectory.gro")
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                let chemfiles_iter = system_chemfiles
                    .traj_iter::<ChemfilesReader>("test_files/protein_trajectory.gro")
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                compare_gro_iterators(gro_iter, chemfiles_iter);
            }
        }
    }

    mod tests_pdb {
        use crate::test_utilities::utilities::compare_box_low_precision;

        use super::*;

        /// Compare two iterators at least one of which is a pdb iterator.
        fn compare_pdb_iterators<'a>(
            mut iter1: impl Iterator<Item = Result<&'a mut System, ReadTrajError>>,
            mut iter2: impl Iterator<Item = Result<&'a mut System, ReadTrajError>>,
        ) {
            while let (Some(frame1), Some(frame2)) = (iter1.next(), iter2.next()) {
                let frame1 = frame1.unwrap();
                let frame2 = frame2.unwrap();

                compare_box_low_precision(frame1.get_box().unwrap(), frame2.get_box().unwrap());

                for (atom1, atom2) in frame1.atoms_iter().zip(frame2.atoms_iter()) {
                    compare_atoms_without_forces(atom1, atom2);
                }
            }

            // check that both operators are exhausted
            assert!(iter1.next().is_none() && iter2.next().is_none());
        }

        #[test]
        fn read_pdb_isolated() {
            let mut system = System::from_file("test_files/protein_trajectory.pdb").unwrap();

            let n_frames = system
                .traj_iter::<ChemfilesReader>("test_files/protein_trajectory.pdb")
                .unwrap()
                .count();

            assert_eq!(n_frames, 11);
        }

        #[test]
        fn read_pdb_pass() {
            let mut system_gro = System::from_file("test_files/protein_trajectory.pdb").unwrap();
            let mut system_pdb = system_gro.clone();

            let gro_iter = system_gro
                .gro_iter("test_files/protein_trajectory.gro")
                .unwrap();
            let pdb_iter = system_pdb
                .traj_iter::<ChemfilesReader>("test_files/protein_trajectory.pdb")
                .unwrap();

            compare_pdb_iterators(gro_iter, pdb_iter);
        }

        #[test]
        fn read_pdb_unmatching() {
            let mut system = System::from_file("test_files/example.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/protein_trajectory.pdb") {
                Err(ReadTrajError::AtomsNumberMismatch(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("PDB file should not be valid."),
            }
        }

        #[test]
        fn read_pdb_nonexistent() {
            let mut system = System::from_file("test_files/protein_trajectory.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/nonexistent.pdb") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("PDB file should not exist."),
            }
        }

        // this test panics because `chemfiles` is buggy as hell
        /*#[test]
        fn read_pdb_not_pdb() {
            let mut system = System::from_file("test_files/protein_trajectory.gro").unwrap();

            match system.traj_iter::<ChemfilesReader>("test_files/fake_pdb.pdb") {
                Err(ReadTrajError::ChemfilesError(_)) => (),
                Err(e) => panic!("Unexpected error `{}` returned.", e),
                Ok(_) => panic!("File should not be a PDB file."),
            }
        }*/

        #[test]
        fn read_pdb_steps() {
            for step in [1, 2, 3, 5, 23] {
                let mut system_gro =
                    System::from_file("test_files/protein_trajectory.gro").unwrap();
                let mut system_pdb = system_gro.clone();

                let gro_iter = system_gro
                    .gro_iter("test_files/protein_trajectory.gro")
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                let pdb_iter = system_pdb
                    .traj_iter::<ChemfilesReader>("test_files/protein_trajectory.pdb")
                    .unwrap()
                    .with_step(step)
                    .unwrap();

                compare_pdb_iterators(gro_iter, pdb_iter);
            }
        }
    }
}
