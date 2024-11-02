// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of functions for reading and writing gro files as trajectories.

use std::io::BufRead;
use std::marker::PhantomData;
use std::path::Path;
use std::str::FromStr;
use std::{fs::File, io::BufReader};

use regex::Regex;

use crate::prelude::{TrajRead, TrajReadOpen, TrajReader, TrajStepRead, Vector3D};
use crate::{
    errors::ReadTrajError,
    prelude::{FrameData, SimBox, TrajFile},
    system::System,
};

#[derive(Debug)]
pub struct GroReader<'a> {
    system: *mut System,
    gro: GroFile,
    phantom: PhantomData<&'a mut System>,
}

#[derive(Debug)]
pub struct GroFile {
    buffer: BufReader<File>,
    filename: Box<Path>,
}

impl TrajFile for GroFile {}

#[derive(Debug)]
pub struct GroFrameData {
    time: Option<f32>,
    step: Option<u64>,
    simbox: SimBox,
    positions: Vec<[f32; 3]>,
    velocities: Vec<Option<[f32; 3]>>,
}

/// Extract time and step from string. Returns `None` if time and step could not be read.
fn extract_time_step(string: &str) -> Option<(f32, u64)> {
    let re = Regex::new(r"t=\s*([\d\.\-]+)\s+step=\s*(\d+)").expect(
        "FATAL GROAN ERROR | gro_io::extract_time_step | Could not construct regular expression.",
    );

    if let Some(caps) = re.captures(string) {
        let time_str = caps.get(1)?.as_str();
        let step_str = caps.get(2)?.as_str();

        let time = f32::from_str(time_str).ok()?;
        let step = u64::from_str(step_str).ok()?;

        Some((time, step))
    } else {
        None
    }
}

/// Attempt to extract frame time and step from the title of the frame.
/// Returns `None` if the time or step is not present.
fn read_time_step(reader: &mut GroFile) -> Result<Option<(f32, u64)>, ReadTrajError> {
    let title = super::get_title(&mut reader.buffer, reader.filename.clone())
        .map_err(|_| ReadTrajError::FrameNotFound)?;

    Ok(extract_time_step(&title))
}

/// Get the position (and velocity, if present) of an atom from a single line of a gro file.
fn read_position_velocity(
    reader: &mut GroFile,
) -> Result<([f32; 3], Option<[f32; 3]>), ReadTrajError> {
    let mut line = String::new();
    reader
        .buffer
        .read_line(&mut line)
        .map_err(|_| ReadTrajError::FrameNotFound)?;

    if line.len() < 44 {
        return Err(ReadTrajError::FrameNotFound);
    }

    let mut position = [0.0f32; 3];
    for (i, item) in position.iter_mut().enumerate() {
        let curr = 20 + i * 8;
        *item = line[curr..curr + 8]
            .trim()
            .parse::<f32>()
            .map_err(|_| ReadTrajError::FrameNotFound)?;
    }

    let velocity = if line.trim_end().len() >= 68 {
        let mut velocity = [0.0f32; 3];

        for (i, item) in velocity.iter_mut().enumerate() {
            let curr = 44 + i * 8;
            *item = line[curr..curr + 8]
                .trim()
                .parse::<f32>()
                .map_err(|_| ReadTrajError::FrameNotFound)?;
        }

        Some(velocity)
    } else {
        None
    };

    Ok((position, velocity))
}

/// Attempt to read the simulation box for the frame.
fn read_box(reader: &mut GroFile) -> Result<SimBox, ReadTrajError> {
    let mut line = String::new();
    reader
        .buffer
        .read_line(&mut line)
        .map_err(|_| ReadTrajError::FrameNotFound)?;

    super::line_as_box(&line).map_err(|_| ReadTrajError::FrameNotFound)
}

impl FrameData for GroFrameData {
    type TrajFile = GroFile;

    fn from_frame(
        traj_file: &mut Self::TrajFile,
        system: &System,
    ) -> Option<Result<Self, crate::errors::ReadTrajError>>
    where
        Self: Sized,
    {
        let (time, step) = match read_time_step(traj_file) {
            Ok(Some((x, y))) => (Some(x), Some(y)),
            Ok(None) => (None, None),
            Err(e) => return Some(Err(e)),
        };

        let n_atoms = match super::get_natoms(&mut traj_file.buffer, traj_file.filename.clone()) {
            Ok(x) => x,
            Err(_) => return Some(Err(ReadTrajError::FrameNotFound)),
        };

        if n_atoms != system.get_n_atoms() {
            return Some(Err(ReadTrajError::AtomsNumberMismatch(
                traj_file.filename.clone(),
            )));
        }

        let mut positions = Vec::with_capacity(n_atoms);
        let mut velocities = Vec::with_capacity(n_atoms);

        for _ in 0..n_atoms {
            match read_position_velocity(traj_file) {
                Ok((pos, vel)) => {
                    positions.push(pos);
                    velocities.push(vel);
                }
                Err(e) => return Some(Err(e)),
            }
        }

        let simbox = match read_box(traj_file) {
            Ok(x) => x,
            Err(e) => return Some(Err(e)),
        };

        Some(Ok(GroFrameData {
            time,
            step,
            simbox,
            positions,
            velocities,
        }))
    }

    fn update_system(self, system: &mut System) {
        for (i, atom) in system.get_atoms_mut().iter_mut().enumerate() {
            let pos = Vector3D::from(unsafe { *self.positions.get_unchecked(i) });
            atom.set_position(pos);

            let vel = unsafe { *self.velocities.get_unchecked(i) };
            match vel {
                Some(x) => atom.set_velocity(Vector3D::from(x)),
                None => atom.reset_velocity(),
            }
        }

        // update the system
        if let Some(x) = self.step {
            system.set_simulation_step(x);
        }

        if let Some(x) = self.time {
            system.set_simulation_time(x);
        }

        system.set_box(self.simbox);
    }
}

impl<'a> TrajRead<'a> for GroReader<'a> {
    type FrameData = GroFrameData;

    fn get_system(&mut self) -> *mut System {
        self.system
    }

    fn get_file_handle(
        &mut self,
    ) -> &mut <<Self as TrajRead<'a>>::FrameData as FrameData>::TrajFile {
        &mut self.gro
    }
}

impl<'a> TrajReadOpen<'a> for GroReader<'a> {
    fn new(system: &'a mut System, filename: impl AsRef<Path>) -> Result<Self, ReadTrajError>
    where
        Self: Sized,
    {
        let file = File::open(&filename)
            .map_err(|_| ReadTrajError::FileNotFound(Box::from(filename.as_ref())))?;

        let buffer = BufReader::new(file);

        Ok(GroReader {
            system: system as *mut System,
            gro: GroFile {
                buffer,
                filename: Box::from(filename.as_ref()),
            },
            phantom: PhantomData,
        })
    }
}

impl<'a> TrajStepRead<'a> for GroReader<'a> {
    fn skip_frame(&mut self) -> Result<bool, ReadTrajError> {
        // read title
        let mut buf = String::new();
        if self
            .gro
            .buffer
            .read_line(&mut buf)
            .map_err(|_| ReadTrajError::SkipFailed)?
            == 0
        {
            return Ok(false);
        }

        let n_atoms = super::get_natoms(&mut self.gro.buffer, self.gro.filename.clone())
            .map_err(|_| ReadTrajError::SkipFailed)?;

        // currently, this works even if the frames we skip over contain inconsistent number of atoms
        // do we want this?

        for _ in 0..(n_atoms + 1) {
            buf.clear();
            if self
                .gro
                .buffer
                .read_line(&mut buf)
                .map_err(|_| ReadTrajError::SkipFailed)?
                == 0
            {
                return Ok(false);
            }
        }

        Ok(true)
    }
}

impl System {
    pub fn gro_iter(
        &mut self,
        filename: impl AsRef<Path>,
    ) -> Result<TrajReader<GroReader>, ReadTrajError> {
        Ok(TrajReader::wrap_traj(GroReader::new(self, filename)?))
    }
}

#[cfg(test)]
mod tests {
    use float_cmp::assert_approx_eq;

    use crate::test_utilities::utilities::{
        compare_atoms, compare_atoms_trr_with_gro, compare_box_low_precision,
    };

    use super::*;

    #[test]
    fn gro_iter() {
        let mut system = System::from_file("test_files/protein_trajectory.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        let expected_times = [
            0.0, 100.0, 200.0, 300.0, 300.0, 500.0, 500.0, 700.0, 800.0, 900.0, 1000.0,
        ];
        let expected_steps = [
            0, 5000, 10000, 15000, 15000, 25000, 25000, 35000, 40000, 45000, 50000,
        ];

        for (i, (frame1, frame2)) in system
            .gro_iter("test_files/protein_trajectory.gro")
            .unwrap()
            .zip(system2.xtc_iter("test_files/short_trajectory.xtc").unwrap())
            .enumerate()
        {
            let frame1 = frame1.unwrap();
            let frame2 = frame2.unwrap();

            assert_approx_eq!(f32, frame1.get_simulation_time(), expected_times[i]);
            assert_eq!(frame1.get_simulation_step(), expected_steps[i]);
            compare_box_low_precision(frame1.get_box().unwrap(), frame2.get_box().unwrap());

            for (a1, a2) in frame1
                .atoms_iter()
                .take(61)
                .zip(frame2.atoms_iter().take(61))
            {
                compare_atoms(a1, a2);
            }
        }
    }

    #[test]
    fn gro_iter_velocities() {
        let mut system = System::from_file("test_files/protein.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        let expected_times = [0.0, 0.0, 480.0];
        let expected_steps = [0, 0, 24000];

        let no_vel_atoms = [vec![30], vec![0, 9], vec![59, 60]];

        for (i, (frame1, frame2)) in system
            .gro_iter("test_files/protein_trajectory_velocities.gro")
            .unwrap()
            .zip(
                system2
                    .trr_iter("test_files/short_trajectory.trr")
                    .unwrap()
                    .with_step(3)
                    .unwrap(),
            )
            .enumerate()
        {
            let frame1 = frame1.unwrap();
            let frame2 = frame2.unwrap();

            assert_approx_eq!(f32, frame1.get_simulation_time(), expected_times[i]);
            assert_eq!(frame1.get_simulation_step(), expected_steps[i]);
            compare_box_low_precision(frame1.get_box().unwrap(), frame2.get_box().unwrap());

            for (a1, a2) in frame1
                .atoms_iter()
                .take(61)
                .zip(frame2.atoms_iter().take(61))
            {
                let mut is_exception = false;
                for exception in &no_vel_atoms[i] {
                    if a1.get_index() == *exception {
                        assert_approx_eq!(
                            f32,
                            a1.get_position().unwrap().x,
                            a2.get_position().unwrap().x,
                            epsilon = 1e-3
                        );

                        assert_approx_eq!(
                            f32,
                            a1.get_position().unwrap().y,
                            a2.get_position().unwrap().y,
                            epsilon = 1e-3
                        );

                        assert_approx_eq!(
                            f32,
                            a1.get_position().unwrap().z,
                            a2.get_position().unwrap().z,
                            epsilon = 1e-3
                        );

                        assert!(!a1.has_velocity());
                        is_exception = true;
                        break;
                    }
                }

                if !is_exception {
                    compare_atoms_trr_with_gro(a1, a2);
                }
            }
        }
    }

    #[test]
    fn gro_iter_with_step() {
        let steps = [2, 3, 4, 5, 7];

        let mut system = System::from_file("test_files/protein_trajectory.gro").unwrap();
        let mut system2 = System::from_file("test_files/example.gro").unwrap();

        for step in steps.into_iter() {
            for (frame1, frame2) in system
                .gro_iter("test_files/protein_trajectory.gro")
                .unwrap()
                .with_step(step)
                .unwrap()
                .zip(
                    system2
                        .xtc_iter("test_files/short_trajectory.xtc")
                        .unwrap()
                        .with_step(step)
                        .unwrap(),
                )
            {
                let frame1 = frame1.unwrap();
                let frame2 = frame2.unwrap();

                compare_box_low_precision(frame1.get_box().unwrap(), frame2.get_box().unwrap());

                for (a1, a2) in frame1
                    .atoms_iter()
                    .take(61)
                    .zip(frame2.atoms_iter().take(61))
                {
                    compare_atoms(a1, a2);
                }
            }
        }
    }

    #[test]
    fn gro_iter_missing_box() {
        let mut system = System::from_file("test_files/protein.gro").unwrap();

        match system
            .gro_iter("test_files/protein_trajectory_missing_box.gro")
            .unwrap()
            .next()
        {
            Some(Ok(_)) => panic!("Function should have failed."),
            Some(Err(ReadTrajError::FrameNotFound)) => (),
            Some(Err(e)) => panic!("Unexpected error type `{}` returned.", e),
            None => panic!("Iterator is empty."),
        }
    }

    #[test]
    fn gro_iter_missing_natoms() {
        let mut system = System::from_file("test_files/protein.gro").unwrap();

        match system
            .gro_iter("test_files/protein_trajectory_missing_natoms.gro")
            .unwrap()
            .skip(1)
            .next()
        {
            Some(Ok(_)) => panic!("Function should have failed."),
            Some(Err(ReadTrajError::FrameNotFound)) => (),
            Some(Err(e)) => panic!("Unexpected error type `{}` returned.", e),
            None => panic!("Iterator is empty."),
        }
    }

    #[test]
    fn gro_iter_invalid_atom_numbers() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        match system
            .gro_iter("test_files/protein_trajectory.gro")
            .unwrap()
            .next()
        {
            Some(Ok(_)) => panic!("Function should have failed."),
            Some(Err(ReadTrajError::AtomsNumberMismatch(x))) => {
                assert_eq!(x.to_str().unwrap(), "test_files/protein_trajectory.gro")
            }
            Some(Err(e)) => panic!("Unexpected error type `{}` returned.", e),
            None => panic!("Iterator is empty."),
        }
    }

    #[test]
    fn gro_iter_missing_title() {
        let mut system = System::from_file("test_files/protein.gro").unwrap();

        match system
            .gro_iter("test_files/protein_trajectory_missing_title.gro")
            .unwrap()
            .skip(1)
            .next()
        {
            Some(Ok(_)) => panic!("Function should have failed."),
            Some(Err(ReadTrajError::FrameNotFound)) => (),
            Some(Err(e)) => panic!("Unexpected error type `{}` returned.", e),
            None => panic!("Iterator is empty."),
        }
    }

    #[test]
    fn gro_iter_incomplete_line() {
        let mut system = System::from_file("test_files/protein.gro").unwrap();

        match system
            .gro_iter("test_files/protein_trajectory_incomplete_line.gro")
            .unwrap()
            .skip(1)
            .next()
        {
            Some(Ok(_)) => panic!("Function should have failed."),
            Some(Err(ReadTrajError::FrameNotFound)) => (),
            Some(Err(e)) => panic!("Unexpected error type `{}` returned.", e),
            None => panic!("Iterator is empty."),
        }
    }
}
