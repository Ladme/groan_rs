// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of ProgressPrinter structure for printing the progress of trajectory reading.

use colored::{ColoredString, Colorize};
use std::io::{self, Write};

/// Progress of trajectory reading.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ProgressStatus {
    /// Trajectory reading is in progress.
    Running,
    /// Trajectory has been read completely.
    Completed,
    /// Trajectory reading failed.
    Failed,
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
            panic!("Groan error. ProgressMessage can not be longer than 9 characters.");
        }

        ProgressMessage { msg: string }
    }

    /// Print formatted `ProgressMessage`.
    fn print(&self, colored: bool) {
        if colored {
            print!("[{: ^9}]   ", self.msg)
        } else {
            print!("[{: ^9}]   ", self.msg.as_ref() as &str)
        }
    }
}

/// Structure handling printing of progress of reading a trajectory file.
#[derive(Debug, PartialEq)]
pub struct ProgressPrinter {
    /// Current status of reading. Default: ProgressStatus::Running.
    status: ProgressStatus,
    /// Frequency of printing. Print every `print_freq`th step. Default: 500,000 steps.
    print_freq: u64,
    /// If true, the output will be colored. Default: true.
    colored: bool,
    /// String to be printed with the current simulation step. Default: "Step".cyan().
    step_msg: ColoredString,
    /// String to be printed with the current simulation step. Default: "Time".bright_purple().
    time_msg: ColoredString,
    /// String to be printed when the trajectory reading is in progress. Default: "RUNNING".yellow().
    running_msg: ProgressMessage,
    /// String to be printed when the trajectory reading is completed. Default: "COMPLETED".green().
    completed_msg: ProgressMessage,
    /// String to be printed when the trajectory reading failed. Default: "FAILED".red().
    failed_msg: ProgressMessage,
}

impl ProgressPrinter {
    /// Create an instance of `ProgressPrinter` with default values.
    pub fn new() -> Self {
        ProgressPrinter {
            status: ProgressStatus::Running,
            print_freq: 500_000,
            colored: true,
            step_msg: "Step".cyan(),
            time_msg: "Time".bright_purple(),
            running_msg: ProgressMessage::new("RUNNING".yellow()),
            completed_msg: ProgressMessage::new("COMPLETED".green()),
            failed_msg: ProgressMessage::new("FAILED".red()),
        }
    }

    /// Create new `ProgressPrinter` with specific value for `status`.
    pub fn with_status(mut self, status: ProgressStatus) -> Self {
        self.status = status;
        self
    }

    /// Set new status to an already constructed `ProgressPrinter`.
    pub fn set_status(&mut self, status: ProgressStatus) {
        self.status = status;
    }

    /// Create new `ProgressPrinter` with specific value for `print_freq`.
    pub fn with_print_freq(mut self, print_freq: u64) -> Self {
        self.print_freq = print_freq;
        self
    }

    /// Create new `ProgressPrinter` with specific value for `colored`.
    pub fn with_colored(mut self, colored: bool) -> Self {
        self.colored = colored;
        self
    }

    /// Create new `ProgressPrinter` with specific value for `step_msg`.
    pub fn with_step_msg(mut self, step_msg: ColoredString) -> Self {
        self.step_msg = step_msg;
        self
    }

    /// Create new `ProgressPrinter` with specific value for `time_msg`.
    pub fn with_time_msg(mut self, time_msg: ColoredString) -> Self {
        self.time_msg = time_msg;
        self
    }

    /// Create new `ProgressPrinter` with specific value for `running_msg`.
    pub fn with_running_msg(mut self, running_msg: ColoredString) -> Self {
        self.running_msg = ProgressMessage::new(running_msg);
        self
    }

    /// Create new `ProgressPrinter` with specific value for `completed_msg`.
    pub fn with_completed_msg(mut self, completed_msg: ColoredString) -> Self {
        self.completed_msg = ProgressMessage::new(completed_msg);
        self
    }

    /// Create new `ProgressPrinter` with specific value for `failed_msg`.
    pub fn with_failed_msg(mut self, failed_msg: ColoredString) -> Self {
        self.failed_msg = ProgressMessage::new(failed_msg);
        self
    }

    /// Print progress info about trajectory reading.
    pub fn print(&self, sim_step: u64, sim_time: u64) {
        if self.status != ProgressStatus::Running || sim_step % self.print_freq == 0 {
            match self.status {
                ProgressStatus::Running => self.running_msg.print(self.colored),
                ProgressStatus::Completed => self.completed_msg.print(self.colored),
                ProgressStatus::Failed => self.failed_msg.print(self.colored),
            }

            if self.colored {
                print!(
                    "{} {:12} | {} {:12} ps\r",
                    self.step_msg, sim_step, self.time_msg, sim_time
                );
            } else {
                print!(
                    "{} {:12} | {} {:12} ps\r",
                    self.step_msg.as_ref() as &str,
                    sim_step,
                    self.time_msg.as_ref() as &str,
                    sim_time
                );
            }

            match self.status {
                ProgressStatus::Running => (),
                ProgressStatus::Completed | ProgressStatus::Failed => println!(),
            }

            io::stdout().flush().unwrap();
        }
    }
}

impl Default for ProgressPrinter {
    fn default() -> Self {
        Self::new()
    }
}
