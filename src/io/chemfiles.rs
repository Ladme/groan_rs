// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Integration with the trajectory readers provided by `chemfiles`.

use std::{marker::PhantomData, path::Path};

use crate::{
    errors::ReadTrajError,
    prelude::{
        FrameData, FrameDataTime, TrajFile, TrajFullReadOpen, TrajRangeRead, TrajRead,
        TrajReadOpen, TrajStepRead, TrajStepTimeRead, Vector3D,
    },
    system::System,
};

/// Used when jumping to the start of iteration.
const TIME_PRECISION: f32 = 0.001;

/// Iterator over a trajectory file using `chemfiles`.
#[derive(Debug)]
pub struct ChemfilesReader<'a> {
    system: *mut System,
    trajectory: ChemfilesTrajectory,
    phantom: PhantomData<&'a mut System>,
}

/// Wrapper around a chemfiles trajectory.
#[derive(Debug)]
pub struct ChemfilesTrajectory {
    path: Box<Path>,
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
    fn set_time(&mut self) {
        match self.frame.get("time") {
            Some(chemfiles::Property::Double(time)) => self.current_time = time as f32,
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
    fn new(filename: impl AsRef<Path>) -> Result<Self, ReadTrajError> {
        let mut traj =
            chemfiles::Trajectory::open(&filename, 'r').map_err(ReadTrajError::ChemfilesError)?;

        // we have to get the number of frames in the trajectory to be able to stop at the appropriate frame,
        // since chemfiles does not return sensible errors
        let max_frames = traj.nsteps();
        Ok(Self {
            path: Box::from(filename.as_ref()),
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
                // we check for mismatch at every frame since chemfiles does not support checking at the start of iteration
                if frame.frame.positions().len() != system.get_n_atoms() {
                    return Some(Err(ReadTrajError::AtomsNumberMismatch(
                        traj_file.path.clone(),
                    )));
                }
                traj_file.current_frame += 1;
                frame.set_time();
                frame.set_precision();
                Some(Ok(frame))
            }
            Err(e) => Some(Err(ReadTrajError::ChemfilesError(e))),
        }
    }

    fn update_system(self, system: &mut System) {
        let positions = self.frame.positions();
        let velocities = self.frame.velocities();
        let atoms = system.atoms_iter_mut();

        if let Some(velocities) = velocities {
            for (atom, (pos, vel)) in atoms.zip(positions.iter().zip(velocities.iter())) {
                unsafe {
                    atom.set_position(Vector3D::new(
                        (*pos.get_unchecked(0) / 10.0) as f32,
                        (*pos.get_unchecked(1) / 10.0) as f32,
                        (*pos.get_unchecked(2) / 10.0) as f32,
                    ));
                    atom.set_velocity(Vector3D::new(
                        (*vel.get_unchecked(0) / 10.0) as f32,
                        (*vel.get_unchecked(1) / 10.0) as f32,
                        (*vel.get_unchecked(2) / 10.0) as f32,
                    ));
                }
                atom.reset_force();
            }
        } else {
            for (atom, pos) in atoms.zip(positions.iter()) {
                unsafe {
                    atom.set_position(Vector3D::new(
                        (*pos.get_unchecked(0) / 10.0) as f32,
                        (*pos.get_unchecked(1) / 10.0) as f32,
                        (*pos.get_unchecked(2) / 10.0) as f32,
                    ));
                }
                atom.reset_velocity();
                atom.reset_force();
            }
        }

        system.set_simulation_step(self.frame.step() as u64);
        system.set_simulation_time(self.current_time);
        system.set_precision(self.precision);
        let cell = self.frame.cell().matrix();
        system.set_box(
            [
                (cell[0][0] / 10.0) as f32,
                (cell[1][1] / 10.0) as f32,
                (cell[2][2] / 10.0) as f32,
                (cell[1][0] / 10.0) as f32,
                (cell[2][0] / 10.0) as f32,
                (cell[0][1] / 10.0) as f32,
                (cell[2][1] / 10.0) as f32,
                (cell[0][2] / 10.0) as f32,
                (cell[1][2] / 10.0) as f32,
            ]
            .into(),
        );
    }
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
            trajectory: ChemfilesTrajectory::new(filename)?,
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
    #[inline]
    fn jump_to_start(&mut self, start_time: f32) -> Result<(), ReadTrajError> {
        // TODO: this reads one more frame than it should
        loop {
            let system = unsafe { &*self.get_system() };
            match ChemfilesFrame::from_frame(&mut self.trajectory, system) {
                None => return Err(ReadTrajError::StartNotFound(start_time.to_string())),
                Some(Err(e)) => return Err(e),
                Some(Ok(frame)) => {
                    if frame.current_time >= start_time - TIME_PRECISION {
                        return Ok(());
                    }
                }
            }
        }
    }
}

impl<'a> TrajStepRead<'a> for ChemfilesReader<'a> {
    #[inline(always)]
    fn skip_frame(&mut self) -> Result<bool, ReadTrajError> {
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

#[cfg(any(feature = "molly", not(feature = "no-xdrfile")))]
#[cfg(test)]
mod tests_read {
    use float_cmp::assert_approx_eq;

    use crate::test_utilities::utilities::{compare_atoms, compare_box};

    use super::*;

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

            for (frame_xtc, frame_chem) in system_xtc
                .xtc_iter(xtc_file)
                .unwrap()
                .zip(system_chem.traj_iter::<ChemfilesReader>(xtc_file).unwrap())
            {
                let frame_xtc = frame_xtc.unwrap();
                let frame_chem = frame_chem.unwrap();

                compare_box(frame_chem.get_box().unwrap(), frame_xtc.get_box().unwrap());
                assert_eq!(frame_chem.get_precision(), frame_xtc.get_precision());
                assert_eq!(
                    frame_chem.get_simulation_step(),
                    frame_xtc.get_simulation_step()
                );
                assert_approx_eq!(
                    f32,
                    frame_chem.get_simulation_time(),
                    frame_xtc.get_simulation_time()
                );

                for (atom_chem, atom_xtc) in frame_chem.atoms_iter().zip(frame_xtc.atoms_iter()) {
                    compare_atoms(atom_chem, atom_xtc);
                }
            }
        }
    }

    #[test]
    fn read_xtc_unmatching() {
        let mut system = System::from_file("test_files/example_novelocities.gro").unwrap();

        let mut iterator = system
            .traj_iter::<ChemfilesReader>("test_files/short_trajectory.xtc")
            .unwrap();
        match iterator.next() {
            None => panic!("Frame should be present."),
            Some(Ok(_)) => panic!("XTC file should not be valid."),
            Some(Err(ReadTrajError::AtomsNumberMismatch(_))) => (),
            Some(Err(e)) => panic!("Unexpected error type `{}` returned.", e),
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

        let mut iterator = system
            .traj_iter::<ChemfilesReader>("test_files/triclinic.gro")
            .unwrap();

        match iterator.next() {
            None => panic!("Frame should be present."),
            Some(Ok(_)) => panic!("XTC file should not be valid."),
            Some(Err(ReadTrajError::ChemfilesError(_))) => (),
            Some(Err(e)) => panic!("Unexpected error type `{}` returned.", e),
        }
    }

    #[test]
    fn read_xtc_range() {
        let mut system_chem = System::from_file("test_files/example.gro").unwrap();
        let mut system_xtc = system_chem.clone();

        for (frame_xtc, frame_chem) in system_xtc
            .xtc_iter("test_files/short_trajectory.xtc")
            .unwrap()
            .with_range(300.0, 800.0)
            .unwrap()
            .zip(
                system_chem
                    .traj_iter::<ChemfilesReader>("test_files/short_trajectory.xtc")
                    .unwrap()
                    .with_range(300.0, 800.0)
                    .unwrap(),
            )
        {
            let frame_xtc = frame_xtc.unwrap();
            let frame_chem = frame_chem.unwrap();

            println!("{:?}", frame_chem.get_simulation_time());
            println!("{:?}", frame_xtc.get_simulation_time());

            compare_box(frame_chem.get_box().unwrap(), frame_xtc.get_box().unwrap());
            assert_eq!(frame_chem.get_precision(), frame_xtc.get_precision());
            assert_eq!(
                frame_chem.get_simulation_step(),
                frame_xtc.get_simulation_step()
            );
            assert_approx_eq!(
                f32,
                frame_chem.get_simulation_time(),
                frame_xtc.get_simulation_time()
            );

            for (atom_chem, atom_xtc) in frame_chem.atoms_iter().zip(frame_xtc.atoms_iter()) {
                compare_atoms(atom_chem, atom_xtc);
            }
        }
    }
}
