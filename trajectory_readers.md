## Guide: Implementing Trajectory Readers in `groan_rs` (v0.7.0+)

### Overview
This guide shows you how to read custom trajectory formats in `groan_rs`. Check `src/io/traj_io.rs` for method signatures and `src/io/{xtc_io.rs, trr_io.rs}` for examples. A "template" for creating your own custom trajectory reader is also provided further down in this guide.

### Steps

1. **Trajectory File Wrapper**
   - Create a wrapper struct for your trajectory file.
   - Implement the empty `TrajFile` trait:
     ```rust
     impl TrajFile for YourTrajectoryFile {}
     ```

2. **Frame Data Holder**
   - Define a struct to hold frame data.
   - Implement `FrameData` trait, specifying your `TrajFile` type.
   - Include:
     - `from_frame`: Reads and stores data from a frame.
     - `update_system`: Transfers data to the `System` struct.

3. **Trajectory Reader**
   - Make a struct that reads the trajectory file and holds a mutable `System` pointer.
   - Implement `TrajRead` and `TrajReadOpen` traits:
     - `TrajRead` requires:
       - `get_system`: Returns a mutable `System` pointer.
       - `get_file_handle`: Gets a reference to the trajectory file.
     - `TrajReadOpen` requires:
       - `new`: Opens the trajectory file and initializes the handler.

   - Use `System::traj_iter::<YourTrajectoryReader>` to iterate frames. Add a `ProgressPrinter` for progress updates.

### Advanced Iteration

- **Iteration in Range**
  - Extend `FrameData` with `FrameDataTime` for time-based data handling.
  - Implement `TrajRangeRead` in your reader for range-specific iterations.
  - Use `with_range` for time range iterations.

- **Iteration with Step**
  - For stepping through frames, implement `TrajStepRead` in your reader.
  - Include `skip_frame` method for frame skipping.
  - Use `with_step` for stepping through frames.

### Concatenating Trajectories

- Implement `FrameDataTime` to enable trajectory concatenation via `System::traj_cat_iter::<YourTrajectoryReader>`.
- For range-limited iteration, use `with_range`.
- To step through concatenated trajectories, implement `TrajStepTimeRead` with `skip_frame_time`, allowing `with_step` usage.

### Template for a Custom Trajectory Reader
```rust

use std::path::Path;

use crate::{errors::ReadTrajError, system::general::System};

use crate::io::traj_io::*;

/*******************************************/
/*              BASIC ITERATION            */
/*******************************************/

/* Trajectory File Wrapper */

struct TrajFilePlaceholder { /* some fields */ }

impl TrajFile for TrajFilePlaceholder {}

/* Frame Data Holder */

struct FrameDataPlaceholder { /* some fields */ }

impl FrameData for FrameDataPlaceholder {
    type TrajFile = TrajFilePlaceholder;

    fn from_frame(
        traj_file: &mut Self::TrajFile,
        system: &System,
    ) -> Option<Result<Self, ReadTrajError>> {
        todo!();
    }

    fn update_system(self, system: &mut System) {
        todo!();
    }
}

/* Trajectory Reader */

struct TrajectoryReaderPlaceholder { /* some fields */ }

impl<'a> TrajRead<'a> for TrajectoryReaderPlaceholder {
    type FrameData = FrameDataPlaceholder;
    
    fn get_system(&mut self) -> *mut System {
        todo!();
    }
    
    fn get_file_handle(
        &mut self,
    ) -> &mut <<Self as TrajRead<'a>>::FrameData as FrameData>::TrajFile {
        todo!();
    }
}

impl<'a> TrajReadOpen<'a> for TrajectoryReaderPlaceholder {
    fn new(system: &'a mut System, filename: impl AsRef<Path>) -> Result<Self, ReadTrajError> {
        todo!();
    }
}

// If you implement the above methods, you should be able to create a trajectory iterator:
let iterator = system.traj_iter::<TrajectoryReaderPlaceholder>("trajectory_placeholder");

/*******************************************/
/*              ITERATION IN RANGE         */
/*******************************************/

/* Frame Data Extension */

impl FrameDataTime for FrameDataPlaceholder {
    fn get_time(&self) -> f32 {
        todo!();
    }
}

/* Trajectory Reader Extension */

impl<'a> TrajRangeRead<'a> for TrajectoryReaderPlaceholder {
    fn jump_to_start(&mut self, start_time: f32) -> Result<(), ReadTrajError> {
        todo!();
    }
}

// If you implement the above methods, you should be able to apply the `with_range` method to your iterator:
let iterator = system
    .traj_iter::<TrajectoryReaderPlaceholder>("trajectory_placeholder")
    .unwrap()
    .with_range(start_time, end_time)
    .unwrap();


// Furthermore, you can now concatenate your trajectories (and apply `with_range` method to them):
let concatenated_iterator = system
    .traj_cat_iter::<TrajectoryReaderPlaceholder>(&vec!["traj1", "traj2", "traj3"])
    .unwrap()
    .with_range(start_time, end_time)
    .unwrap();


/*******************************************/
/*             ITERATION WITH STEP         */
/*******************************************/

/* Trajectory Reader Extension */

impl<'a> TrajStepRead<'a> for TrajectoryReaderPlaceholder {
    fn skip_frame(&mut self) -> Result<bool, ReadTrajError> {
        todo!();
    }
}

// If you implement the above method, you should be able to apply the `with_step` method to your iterator:
let iterator = system
    .traj_iter::<TrajectoryReaderPlaceholder>("trajectory_placeholder")
    .unwrap()
    .with_step(step)
    .unwrap();

/*******************************************/
/*       CONCATENATING TRAJECTORIES        */
/*******************************************/

// If you want to use the `with_step` method for your concatenated trajectories, you must implement the following method:

impl<'a >TrajStepTimeRead<'a> for TrajectoryReaderPlaceholder {
    fn skip_frame_time(&mut self) -> Result<Option<f32>, ReadTrajError> {
        todo!();
    }
}

// Now you should be able to apply the `with_step` method to your concatenated trajectories:
let concatenated_iterator = system
    .traj_cat_iter::<TrajectoryReaderPlaceholder>(&vec!["traj1", "traj2", "traj3"])
    .unwrap()
    .with_step(step)
    .unwrap();
```