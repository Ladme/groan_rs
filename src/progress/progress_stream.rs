// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of ProgressStream handling writing for single-threaded and multi-threaded applications.

#[cfg(feature = "parallel")]
pub(super) use multi_threaded::ProgressStream;
#[cfg(not(feature = "parallel"))]
pub(super) use single_threaded::ProgressStream;

#[cfg(not(feature = "parallel"))]
mod single_threaded {
    use std::{cell::RefCell, io::Write, rc::Rc};

    #[derive(Clone)]
    pub(in crate::progress) struct ProgressStream(Rc<RefCell<dyn Write>>);

    impl ProgressStream {
        #[inline(always)]
        pub(in crate::progress) fn new(writer: Box<dyn Write>) -> Self {
            Self(Rc::from(RefCell::from(writer)))
        }

        /// ## Panics
        /// Panics if the writing fails.
        #[inline(always)]
        pub(in crate::progress) fn write(&self, message: &str) {
            write!(self.0.borrow_mut(), "{}", message).expect("FATAL GROAN ERROR | single_threaded::ProgressStream::write | Could not write to stream.");
        }

        /// ## Panics
        /// Panics if the writing fails.
        #[inline(always)]
        pub(in crate::progress) fn write_flush(&self, message: &str) {
            let mut writer = self.0.borrow_mut();
            write!(writer, "{}", message).expect("FATAL GROAN ERROR | single_threaded::ProgressStream::write_flush | Could not write to stream.");
            writer.flush().expect("FATAL GROAN ERROR | single_threaded::ProgressStream::write_flush | Could not flush stream.");
        }
    }

    impl Default for ProgressStream {
        #[inline(always)]
        fn default() -> Self {
            Self(Rc::from(RefCell::from(std::io::stdout())))
        }
    }
}

#[cfg(feature = "parallel")]
mod multi_threaded {
    use parking_lot::Mutex;
    use std::{io::Write, sync::Arc};

    #[derive(Clone)]
    pub(in crate::progress) struct ProgressStream(Arc<Mutex<dyn Write>>);

    impl ProgressStream {
        #[inline(always)]
        pub(in crate::progress) fn new(writer: Box<dyn Write>) -> Self {
            Self(Arc::from(Mutex::from(writer)))
        }

        /// ## Panics
        /// Panics if the writing fails.
        #[inline(always)]
        pub(in crate::progress) fn write(&self, message: &str) {
            let mut writer = self.0.lock();
            write!(writer, "{}", message).expect("FATAL GROAN ERROR | multi_threaded::ProgressStream::write | Could not write to stream.");
        }

        /// ## Panics
        /// Panics if the writing or flushing fails.
        #[inline(always)]
        pub(in crate::progress) fn write_flush(&self, message: &str) {
            let mut writer = self.0.lock();
            write!(writer, "{}", message).expect("FATAL GROAN ERROR | multi_threaded::ProgressStream::write_flush | Could not write to stream.");
            writer.flush().expect("FATAL GROAN ERROR | multi_threaded::ProgressStream::write_flush | Could not flush stream.");
        }
    }

    impl Default for ProgressStream {
        #[inline(always)]
        fn default() -> Self {
            Self(Arc::from(Mutex::from(std::io::stdout())))
        }
    }
}
