// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of ProgressPrinter structure for printing the progress of trajectory reading.

use colored::{ColoredString, Colorize};
use std::{io::Write, num::NonZeroUsize};

mod progress_stream;

use progress_stream::ProgressStream;

/// Progress of trajectory reading.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProgressStatus {
    /// Trajectory reading is in progress.
    Running,
    /// Trajectory has been read completely.
    Completed,
    /// Trajectory reading failed.
    Failed,
    /// Trajectory reader is currently jumping towards the start of the iteration.
    Jumping,
}

/// String that can be used inside `ProgressPrinter`.
#[derive(Debug, Clone, PartialEq)]
struct ProgressMessage {
    msg: ColoredString,
}

impl ProgressMessage {
    /// Create new `ProgressMessage`.
    ///
    /// ## Panics
    /// Panics if the string is longer than 9 characters.
    fn new(string: ColoredString) -> Self {
        if string.chars().count() > 9 {
            panic!("FATAL GROAN ERROR | ProgressMessage::new | `ProgressMessage` can not be longer than 9 characters.");
        }

        ProgressMessage { msg: string }
    }

    /// Print formatted `ProgressMessage`.
    fn print(&self, out: &ProgressStream, colored: bool) {
        let message = if colored {
            format!("[{: ^9}]   ", self.msg)
        } else {
            format!("[{: ^9}]   ", self.msg.as_ref() as &str)
        };

        out.write(&message);
    }
}

impl From<&str> for ProgressMessage {
    fn from(string: &str) -> Self {
        ProgressMessage::new(string.into())
    }
}

impl From<String> for ProgressMessage {
    fn from(string: String) -> Self {
        ProgressMessage::new(string.as_str().into())
    }
}

impl From<ColoredString> for ProgressMessage {
    fn from(string: ColoredString) -> Self {
        ProgressMessage::new(string)
    }
}

/// Structure handling printing of progress of reading a trajectory file.
/// Constructed using `ProgressPrinter::new()` and associated with the given
/// trajectory reader using `TrajMasterRead::print_progress()`.
#[derive(Clone)]
pub struct ProgressPrinter {
    /// Stream to write the progress info to.
    output: ProgressStream,
    /// Current status of reading. Default: ProgressStatus::Running.
    status: ProgressStatus,
    /// Frequency of printing. Print every `print_freq`th frame. Default: 100 frames.
    print_freq: NonZeroUsize,
    /// If true, the output will be colored. Default: true.
    colored: bool,
    /// String to be printed with the current simulation step. Default: "Step".cyan().
    step_msg: ColoredString,
    /// String to be printed with the current simulation time. Default: "Time".bright_purple().
    time_msg: ColoredString,
    /// String to be printed when the trajectory reading is in progress. Default: "RUNNING".yellow().
    running_msg: ProgressMessage,
    /// String to be printed when the trajectory reading is completed. Default: "COMPLETED".green().
    completed_msg: ProgressMessage,
    /// String to be printed when the trajectory reading failed. Default: "FAILED!".red().
    failed_msg: ProgressMessage,
    /// String to be printed when the trajectory reader is jumping to the start of iteration. Default: "JUMPING".bright_purple().
    jumping_msg: ProgressMessage,
    /// String terminating the progress message. Default: `\r` (carriage return).
    terminating: String,
    /// Should newline character be printed at the end of the iteration? Default: `true`.
    newline_at_end: bool,
}

#[cfg(feature = "parallel")]
unsafe impl Send for ProgressPrinter {}
#[cfg(feature = "parallel")]
unsafe impl Sync for ProgressPrinter {}

impl ProgressPrinter {
    /// Create an instance of `ProgressPrinter` with default parameters.
    ///
    /// The default values of the `ProgressPrinter` parameters.
    /// - `output`: `std::io::stdout()` (stream to write the progress info to)
    /// - `status`: `ProgressStatus::Running` (current status of trajectory file reading)
    /// - `print_freq`: `100` (progress info will be printed out every 100 trajectory frames read)
    /// - `colored`: `true` (should the output be colored?)
    /// - `step_msg`: `"Step".cyan()` (string associated with information about the simulation step of the frame)
    /// - `time_msg`: `"Time".bright_purple()` (string associated with information about the time of the frame)
    /// - `running_msg`: `"RUNNING".yellow()` (string printed when the trajectory reading is running)
    /// - `completed_msg`: `"COMPLETED".green()` (string printed when the trajectory reading is completed)
    /// - `failed_msg`: `"FAILED!".red()` (string printed when the trajectory reading failed)
    /// - `jumping_msg`: `"JUMPING".bright_purple()` (string printed when the trajectory reader is jumping to the iteration start)
    /// - `terminating`: `\r` (string terminating the progress message; useful to set to `\n` when printing to a file)
    /// - `newline_at_end`: `true` (should newline be printed once iteration is finished or failed?;
    ///    useful to set to `false` when reading multiple trajectories sequentially)
    ///
    /// You can set custom values for any of the parameters by using `with_%PARAMETER()` method
    /// when constructing the `ProgressPrinter`.
    ///
    /// ## Examples
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// # use colored::Colorize;
    /// #
    /// // create progress printer
    /// // print info every 200th trajectory frame read
    /// // yellow string "ANALYZING" will be printed while the trajectory is being read
    /// // blue string "DONE" will be printed once reading is completed
    /// // other parameters of progress printing are default
    /// let printer = ProgressPrinter::new()
    ///     .with_print_freq(200)
    ///     .with_running_msg("ANALYZING".yellow())
    ///     .with_completed_msg("DONE".blue());
    ///
    /// // load gro file
    /// let mut system = System::from_file("system.gro").unwrap();
    ///
    /// // read the trajectory
    /// for frame in system
    ///     .xtc_iter("trajectory.xtc")
    ///     .unwrap()
    ///     // associate the `ProgressPrinter` with the trajectory iterator
    ///     .print_progress(printer)
    /// {
    ///     let frame = frame.unwrap();
    ///     // analyze the frame
    /// }
    /// ```
    ///
    /// By default, `ProgressPrinter` prints to standard output.
    /// However, you can also let it print into a file.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let file = std::fs::File::create("progress.log").unwrap();
    /// let printer = ProgressPrinter::new()
    ///     .with_output(Box::from(file))
    ///     // turning off colored output which does not make sense for a file
    ///     .with_colored(false)
    ///     // setting terminating string to `newline` character
    ///     .with_terminating("\n");
    /// ```
    pub fn new() -> Self {
        ProgressPrinter {
            output: ProgressStream::default(),
            status: ProgressStatus::Running,
            print_freq: NonZeroUsize::new(100).unwrap(),
            colored: true,
            step_msg: "Step".cyan(),
            time_msg: "Time".bright_purple(),
            running_msg: ProgressMessage::new("RUNNING".yellow()),
            completed_msg: ProgressMessage::new("COMPLETED".green()),
            failed_msg: ProgressMessage::new("FAILED!".red()),
            jumping_msg: ProgressMessage::new("JUMPING".bright_purple()),
            terminating: String::from("\r"),
            newline_at_end: true,
        }
    }

    /// Create new `ProgressPrinter` with specific `output` stream.
    #[inline(always)]
    pub fn with_output(mut self, stream: Box<dyn Write>) -> Self {
        self.output = ProgressStream::new(stream);
        self
    }

    /// Create new `ProgressPrinter` with specific value for `status`.
    #[inline(always)]
    pub fn with_status(mut self, status: ProgressStatus) -> Self {
        self.status = status;
        self
    }

    /// Set new status to an already constructed `ProgressPrinter`.
    #[inline(always)]
    pub fn set_status(&mut self, status: ProgressStatus) {
        self.status = status;
    }

    /// Create new `ProgressPrinter` with specific value for `print_freq`.
    ///
    /// ## Panics
    /// Panics if `print_freq` is zero.
    #[inline(always)]
    pub fn with_print_freq(mut self, print_freq: usize) -> Self {
        self.print_freq = NonZeroUsize::new(print_freq).expect(
            "FATAL GROAN ERROR | ProgressPrinter::with_print_freq | `print_freq` must be non-zero.",
        );
        self
    }

    /// Create new `ProgressPrinter` with specific value for `colored`.
    #[inline(always)]
    pub fn with_colored(mut self, colored: bool) -> Self {
        self.colored = colored;
        self
    }

    /// Create new `ProgressPrinter` with specific value for `step_msg`.
    #[inline(always)]
    pub fn with_step_msg(mut self, step_msg: ColoredString) -> Self {
        self.step_msg = step_msg;
        self
    }

    /// Create new `ProgressPrinter` with specific value for `time_msg`.
    #[inline(always)]
    pub fn with_time_msg(mut self, time_msg: ColoredString) -> Self {
        self.time_msg = time_msg;
        self
    }

    /// Create new `ProgressPrinter` with specific value for `running_msg`.
    ///
    /// ## Panics
    /// Panics if the `running_msg` is longer than 9 characters.
    #[inline(always)]
    pub fn with_running_msg(mut self, running_msg: ColoredString) -> Self {
        self.running_msg = ProgressMessage::new(running_msg);
        self
    }

    /// Create new `ProgressPrinter` with specific value for `completed_msg`.
    ///
    /// ## Panics
    /// Panics if the `completed_msg` is longer than 9 characters.
    #[inline(always)]
    pub fn with_completed_msg(mut self, completed_msg: ColoredString) -> Self {
        self.completed_msg = ProgressMessage::new(completed_msg);
        self
    }

    /// Create new `ProgressPrinter` with specific value for `failed_msg`.
    ///
    /// ## Panics
    /// Panics if the `failed_msg` is longer than 9 characters.
    #[inline(always)]
    pub fn with_failed_msg(mut self, failed_msg: ColoredString) -> Self {
        self.failed_msg = ProgressMessage::new(failed_msg);
        self
    }

    /// Create new `ProgressPrinter` with specific value for `jumping_msg`.
    ///
    /// ## Panics
    /// Panics if the `jumping_msg` is longer than 9 characters.
    #[inline(always)]
    pub fn with_jumping_msg(mut self, jumping_msg: ColoredString) -> Self {
        self.jumping_msg = ProgressMessage::new(jumping_msg);
        self
    }

    /// Create new `ProgressPrinter` with specific value for `terminating`.
    #[inline(always)]
    pub fn with_terminating(mut self, string: &str) -> Self {
        self.terminating = string.to_string();
        self
    }

    /// Create new `ProgressPrinter` that will (not) print newline once iteration is finished.
    #[inline(always)]
    pub fn with_newline_at_end(mut self, newline: bool) -> Self {
        self.newline_at_end = newline;
        self
    }

    /// Prints message about jumping to the start of iteration.
    #[inline(always)]
    fn print_jumping_msg(&mut self) {
        self.jumping_msg.print(&self.output, self.colored);
        let message = format!(
            "Jumping to the start of the iteration...{}",
            self.terminating
        );
        self.output.write_flush(&message);
    }

    /// Print formatted information about the current step and simulation time.
    #[inline(always)]
    fn print_steptime(&mut self, sim_step: u64, sim_time: f32) {
        let message = if self.colored {
            format!(
                "{} {:12} | {} {:12} ps{}",
                self.step_msg, sim_step, self.time_msg, sim_time as u64, self.terminating
            )
        } else {
            format!(
                "{} {:12} | {} {:12} ps{}",
                self.step_msg.as_ref() as &str,
                sim_step,
                self.time_msg.as_ref() as &str,
                sim_time as u64,
                self.terminating
            )
        };

        self.output.write_flush(&message);
    }

    /// Print newline at the end of iteration.
    #[inline(always)]
    fn print_newline_at_end(&mut self) {
        self.output.write_flush("\n");
    }

    /// Print progress info about trajectory reading.
    pub fn print(&mut self, frame_number: usize, sim_step: u64, sim_time: f32) {
        if self.status != ProgressStatus::Running || frame_number % self.print_freq == 0 {
            match self.status {
                ProgressStatus::Running => self.running_msg.print(&self.output, self.colored),
                ProgressStatus::Completed => self.completed_msg.print(&self.output, self.colored),
                ProgressStatus::Failed => self.failed_msg.print(&self.output, self.colored),
                ProgressStatus::Jumping => {
                    self.print_jumping_msg();
                    return;
                }
            }

            self.print_steptime(sim_step, sim_time);

            match self.status {
                ProgressStatus::Running | ProgressStatus::Jumping => (),
                ProgressStatus::Completed | ProgressStatus::Failed => {
                    if self.newline_at_end {
                        self.print_newline_at_end();
                    }
                }
            }
        }
    }
}

impl Default for ProgressPrinter {
    fn default() -> Self {
        Self::new()
    }
}

/******************************/
/*         UNIT TESTS         */
/******************************/

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use tempfile::{tempfile, NamedTempFile};

    #[test]
    fn new() {
        let printer = ProgressPrinter::new();

        assert_eq!(printer.status, ProgressStatus::Running);
        assert_eq!(printer.print_freq.get(), 100);
        assert!(printer.colored);
        assert_eq!(printer.step_msg, "Step".cyan());
        assert_eq!(printer.time_msg, "Time".bright_purple());
        assert_eq!(
            printer.running_msg,
            ProgressMessage::new("RUNNING".yellow())
        );
        assert_eq!(
            printer.completed_msg,
            ProgressMessage::new("COMPLETED".green())
        );
        assert_eq!(printer.failed_msg, ProgressMessage::new("FAILED!".red()));
        assert_eq!(
            printer.jumping_msg,
            ProgressMessage::new("JUMPING".bright_purple())
        );
    }

    #[test]
    fn default() {
        let printer = ProgressPrinter::default();

        assert_eq!(printer.status, ProgressStatus::Running);
        assert_eq!(printer.print_freq.get(), 100);
        assert!(printer.colored);
        assert_eq!(printer.step_msg, "Step".cyan());
        assert_eq!(printer.time_msg, "Time".bright_purple());
        assert_eq!(
            printer.running_msg,
            ProgressMessage::new("RUNNING".yellow())
        );
        assert_eq!(
            printer.completed_msg,
            ProgressMessage::new("COMPLETED".green())
        );
        assert_eq!(printer.failed_msg, ProgressMessage::new("FAILED!".red()));
        assert_eq!(
            printer.jumping_msg,
            ProgressMessage::new("JUMPING".bright_purple())
        );
    }

    #[test]
    fn set_status() {
        let mut printer = ProgressPrinter::new();

        printer.set_status(ProgressStatus::Failed);
        assert_eq!(printer.status, ProgressStatus::Failed);

        printer.set_status(ProgressStatus::Completed);
        assert_eq!(printer.status, ProgressStatus::Completed);

        printer.set_status(ProgressStatus::Running);
        assert_eq!(printer.status, ProgressStatus::Running);

        printer.set_status(ProgressStatus::Jumping);
        assert_eq!(printer.status, ProgressStatus::Jumping);
    }

    #[test]
    fn new_complex() {
        let tmp_file = tempfile().unwrap();

        let printer = ProgressPrinter::new()
            .with_output(Box::from(tmp_file))
            .with_status(ProgressStatus::Jumping)
            .with_print_freq(200)
            .with_colored(false)
            .with_step_msg("STEP".into())
            .with_time_msg("time".yellow())
            .with_running_msg("ANALYZING".red())
            .with_completed_msg("DONE".green())
            .with_failed_msg("FAILURE".on_bright_red())
            .with_jumping_msg("JUMP".underline());

        assert_eq!(printer.status, ProgressStatus::Jumping);
        assert_eq!(printer.print_freq.get(), 200);
        assert!(!printer.colored);
        assert_eq!(printer.step_msg, "STEP".into());
        assert_eq!(printer.time_msg, "time".yellow());
        assert_eq!(printer.running_msg, ProgressMessage::new("ANALYZING".red()));
        assert_eq!(printer.completed_msg, ProgressMessage::new("DONE".green()));
        assert_eq!(
            printer.failed_msg,
            ProgressMessage::new("FAILURE".on_bright_red())
        );
        assert_eq!(
            printer.jumping_msg,
            ProgressMessage::new("JUMP".underline())
        );
    }

    #[test]
    #[should_panic(
        expected = "FATAL GROAN ERROR | ProgressMessage::new | `ProgressMessage` can not be longer than 9 characters."
    )]
    fn progress_message_panic() {
        let _msg = ProgressMessage::new("SHOULD_PANIC".red());
    }

    #[test]
    fn print() {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_owned();

        let mut printer = ProgressPrinter::new()
            .with_output(Box::from(output))
            .with_colored(false);

        printer.set_status(ProgressStatus::Jumping);
        printer.print(0, 0, 0.0);
        printer.set_status(ProgressStatus::Running);
        printer.print(0, 0, 0.0);
        printer.print(1, 10, 10.0);
        printer.print(2, 20, 20.0);
        printer.print(5, 50, 50.0);
        printer.print(95, 950, 950.0);
        printer.print(100, 1000, 1000.0);
        printer.print(101, 1010, 1010.0);
        printer.print(200, 2000, 2000.0);
        printer.print(300, 3000, 3000.0);
        printer.set_status(ProgressStatus::Completed);
        printer.print(400, 4000, 4000.0);
        printer.set_status(ProgressStatus::Failed);
        printer.print(500, 5000, 5000.0);

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/progress_expected.txt").unwrap();
        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn print_with_newline() {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_owned();

        let mut printer = ProgressPrinter::new()
            .with_output(Box::from(output))
            .with_colored(false)
            .with_terminating("\n");

        printer.set_status(ProgressStatus::Jumping);
        printer.print(0, 0, 0.0);
        printer.set_status(ProgressStatus::Running);
        printer.print(0, 0, 0.0);
        printer.print(1, 10, 10.0);
        printer.print(2, 20, 20.0);
        printer.print(5, 50, 50.0);
        printer.print(95, 950, 950.0);
        printer.print(100, 1000, 1000.0);
        printer.print(101, 1010, 1010.0);
        printer.print(200, 2000, 2000.0);
        printer.print(300, 3000, 3000.0);
        printer.set_status(ProgressStatus::Completed);
        printer.print(400, 4000, 4000.0);
        printer.set_status(ProgressStatus::Failed);
        printer.print(500, 5000, 5000.0);

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/progress_expected_newline.txt").unwrap();
        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn print_with_terminating() {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_owned();

        let mut printer = ProgressPrinter::new()
            .with_output(Box::from(output))
            .with_colored(false)
            .with_terminating("  |  ");

        printer.set_status(ProgressStatus::Jumping);
        printer.print(0, 0, 0.0);
        printer.set_status(ProgressStatus::Running);
        printer.print(0, 0, 0.0);
        printer.print(1, 10, 10.0);
        printer.print(2, 20, 20.0);
        printer.print(5, 50, 50.0);
        printer.print(95, 950, 950.0);
        printer.print(100, 1000, 1000.0);
        printer.print(101, 1010, 1010.0);
        printer.print(200, 2000, 2000.0);
        printer.print(300, 3000, 3000.0);
        printer.set_status(ProgressStatus::Completed);
        printer.print(400, 4000, 4000.0);
        printer.set_status(ProgressStatus::Failed);
        printer.print(500, 5000, 5000.0);

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/progress_expected_terminating.txt").unwrap();
        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    fn print_with_newline_at_end_false() {
        let output = NamedTempFile::new().unwrap();
        let path_to_output = output.path().to_owned();

        let mut printer = ProgressPrinter::new()
            .with_output(Box::from(output))
            .with_colored(false)
            .with_newline_at_end(false);

        printer.set_status(ProgressStatus::Jumping);
        printer.print(0, 0, 0.0);
        printer.set_status(ProgressStatus::Running);
        printer.print(0, 0, 0.0);
        printer.print(1, 10, 10.0);
        printer.print(2, 20, 20.0);
        printer.print(5, 50, 50.0);
        printer.print(95, 950, 950.0);
        printer.print(100, 1000, 1000.0);
        printer.print(101, 1010, 1010.0);
        printer.print(200, 2000, 2000.0);
        printer.print(300, 3000, 3000.0);
        printer.set_status(ProgressStatus::Completed);
        printer.print(400, 4000, 4000.0);
        printer.set_status(ProgressStatus::Failed);
        printer.print(500, 5000, 5000.0);

        let mut result = File::open(path_to_output).unwrap();
        let mut expected = File::open("test_files/progress_expected_no_newline.txt").unwrap();
        assert!(file_diff::diff_files(&mut result, &mut expected));
    }

    #[test]
    #[should_panic(
        expected = "FATAL GROAN ERROR | ProgressPrinter::with_print_freq | `print_freq` must be non-zero."
    )]
    fn print_freq_zero_panic() {
        let _printer = ProgressPrinter::new().with_print_freq(0);
    }
}
