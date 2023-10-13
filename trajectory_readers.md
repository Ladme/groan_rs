## Implementing Trajectory Readers for the `groan_rs` Library

**Note**: This guide is applicable to `groan_rs` version 0.4.0 onwards.

For method signatures, refer to `src/io/traj_io.rs`, and see examples of their implementation in `src/io/xtc_io.rs` and `src/io/trr_io.rs`.

### Step 1: Basic Iteration

To read a custom trajectory format in the `groan_rs` library, you need to define three structures:

1. **Trajectory File Wrapper**:
    - A structure that acts as a wrapper around the trajectory file.
    - It should implement the `TrajFile` trait, which is an empty trait:
    ```rust
    impl TrajFile for YourTrajectoryFile {}
    ```

2. **Frame Data Holder**:
    - A structure that stores data from individual frames of the trajectory.
    - It must implement the `FrameData` trait.
    - You need to provide the type of `TrajFile` that the reader will use.
    - Implement the following methods:
        - `from_frame`: This method defines how data from a frame in the trajectory file should be read and subsequently stored within the `FrameData` structure.
        - `update_system`: This method defines how the data stored within the `FrameData` structure should be transferred to the central `System` structure.

3. **Trajectory File Handler**:
    - A structure that manages the trajectory file and maintains a mutable pointer to the `System` structure.
    - It must implement the `TrajRead` trait.
    - Provide the type of `FrameData` your reader will associate with.
    - Implement the following methods:
        - `new`: This method initializes and opens the trajectory file. It also creates a valid instance of the structure implementing the `TrajRead` trait.
        - `get_system`: Returns a mutable pointer to the `System` structure. Using trajectory iterators inherently involves using `unsafe` operations.
        - `get_file_handle`: Returns a mutable reference to the trajectory file currently being read.

After defining and implementing these traits, you can invoke the `System::traj_iter::<YourTrajectoryReader>` method to iterate over the frames of your trajectory.

### Step 2a: Range-Based Iteration

For iterations based on specific time ranges:

1. **Frame Data Extension**:
    - The `FrameData` structure should also implement the `FrameDataTime` trait.
    - Implement the `get_time` method which retrieves the simulation time from the frame data. Ensure your frames are equipped with simulation time data for this functionality.

2. **Trajectory Reader Extension**:
    - The `TrajRead` structure should also implement the `TrajRangeRead` trait.
    - Implement the `jump_to_start` method. It's responsible for efficiently navigating through the trajectory file to reach the starting time of a specified range.

After adding these functionalities, you can use the `with_range` method on your trajectory reader iterator to iterate over specific time ranges.

### Step 2b: Step-Based Iteration

For iterations that process every Nth frame:

1. **Trajectory Reader Step Extension**:
    - Your `TrajRead` structure should also implement the `TrajStepRead` trait.
    - Implement the `skip_frame` method. It's tasked with efficiently skipping over single frames in the trajectory.

Once you've implemented this trait, you can use the `with_step` method on your trajectory reader iterator. Importantly, there's no requirement to implement the `TrajRangeRead` trait (from Step 2a) to leverage the `with_step` method.
