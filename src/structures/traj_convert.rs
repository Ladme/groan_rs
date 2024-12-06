// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of trajectory converters and analyzers.

use std::marker::PhantomData;

use crate::{
    errors::{TrajAnalysisError, TrajConvertAnalysisError, TrajConvertError},
    prelude::TrajMasterRead,
    system::System,
};

/// A generic structure for trajectory converters.
///
/// Trajectory converter can be built on top of any trajectory iterator.
/// Converter provides an interface for iterating through a trajectory, modifying
/// each frame in a defined way, and returning the modified frames.
pub struct TrajConverter<'a, Reader, Converter>
where
    Reader: ConvertableTrajRead<'a>,
    Converter: FrameConvert,
{
    reader: Reader,
    converter: Converter,
    _phantom: PhantomData<&'a Reader>,
}

/// Trait implemented by structures that can be used as trajectory converters.
pub trait FrameConvert {
    type Error: std::error::Error + Send + Sync + Eq + PartialEq + 'static;

    /// Modify the `System` using the converter.
    /// Converter can also be modified by this method.
    fn convert(&mut self, system: &mut System) -> Result<(), Self::Error>;
}

impl<'a, Reader, Converter> Iterator for TrajConverter<'a, Reader, Converter>
where
    Reader: ConvertableTrajRead<'a>,
    Converter: FrameConvert,
{
    type Item = Result<&'a mut System, TrajConvertError<Converter::Error>>;

    /// Read next frame of the trajectory and convert it using the associated converter.
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let frame = match self.reader.next()? {
            Ok(x) => x,
            Err(e) => return Some(Err(TrajConvertError::ReadingError(e))),
        };

        match self.converter.convert(frame) {
            Ok(_) => Some(Ok(frame)),
            Err(e) => Some(Err(TrajConvertError::ConversionError(e))),
        }
    }
}

/// A generic structure for trajectory analyzers.
///
/// Trajectory analyzer can be built on top of any trajectory iterator.
/// Analyzer allows iterating through a trajectory while analyzing each frame
/// and returning both the frame and a computed property associated with it.
pub struct TrajAnalyzer<'a, Reader, Analyzer>
where
    Reader: ConvertableTrajRead<'a>,
    Analyzer: FrameAnalyze,
{
    reader: Reader,
    analyzer: Analyzer,
    _phantom: PhantomData<&'a Reader>,
}

/// Trait implemented by structures that can be used as trajectory analyzers.
pub trait FrameAnalyze {
    type Error: std::error::Error + Send + Sync + Eq + PartialEq + 'static;
    type AnalysisResult;

    /// Analyze the `System` using the analyzer.
    /// Analyzer can be modified by this method.
    fn analyze(&mut self, system: &System) -> Result<Self::AnalysisResult, Self::Error>;
}

impl<'a, Reader, Analyzer> Iterator for TrajAnalyzer<'a, Reader, Analyzer>
where
    Reader: ConvertableTrajRead<'a>,
    Analyzer: FrameAnalyze,
{
    type Item =
        Result<(&'a mut System, Analyzer::AnalysisResult), TrajAnalysisError<Analyzer::Error>>;

    /// Read next frame of the trajectory and analyze it using the associated analyzer.
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let frame = match self.reader.next()? {
            Ok(x) => x,
            Err(e) => return Some(Err(TrajAnalysisError::ReadingError(e))),
        };

        match self.analyzer.analyze(frame) {
            Ok(x) => Some(Ok((frame, x))),
            Err(e) => Some(Err(TrajAnalysisError::AnalysisError(e))),
        }
    }
}

/// A generic structure for trajectory converter-analyzers.
///
/// Trajectory converter-analyzers combine the functionality of both
/// trajectory converters and analyzers. They allow modifying each frame
/// in a trajectory and simultaneously analyzing the frame to return
/// both the modified frame and the associated analysis result.
pub struct TrajConverterAnalyzer<'a, Reader, ConverterAnalyzer>
where
    Reader: ConvertableTrajRead<'a>,
    ConverterAnalyzer: FrameConvertAnalyze,
{
    reader: Reader,
    converter_analyzer: ConverterAnalyzer,
    _phantom: PhantomData<&'a Reader>,
}

/// Trait implemented by structures that can be used as trajectory converter-analyzers.
pub trait FrameConvertAnalyze {
    type Error: std::error::Error + Send + Sync + Eq + PartialEq + 'static;
    type AnalysisResult;

    /// Convert the `System` and analyze it (in any order).
    fn convert_analyze(&mut self, system: &mut System)
        -> Result<Self::AnalysisResult, Self::Error>;
}

impl<'a, Reader, ConverterAnalyzer> Iterator
    for TrajConverterAnalyzer<'a, Reader, ConverterAnalyzer>
where
    Reader: ConvertableTrajRead<'a>,
    ConverterAnalyzer: FrameConvertAnalyze,
{
    type Item = Result<
        (&'a mut System, ConverterAnalyzer::AnalysisResult),
        TrajConvertAnalysisError<ConverterAnalyzer::Error>,
    >;

    /// Read next frame of the trajectory, convert it and analyze it using the associated converter-analyzer.
    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let frame = match self.reader.next()? {
            Ok(x) => x,
            Err(e) => return Some(Err(TrajConvertAnalysisError::ReadingError(e))),
        };

        match self.converter_analyzer.convert_analyze(frame) {
            Ok(x) => Some(Ok((frame, x))),
            Err(e) => Some(Err(TrajConvertAnalysisError::ConversionAnalysisError(e))),
        }
    }
}

/// Implemented by all trajectory readers that can be turned into converters/analyzers.
pub trait ConvertableTrajRead<'a>: TrajMasterRead<'a> + Sized {
    /// Transform the trajectory reader into trajectory converter.
    #[inline(always)]
    fn convert<C>(self, converter: C) -> TrajConverter<'a, Self, C>
    where
        C: FrameConvert,
    {
        TrajConverter {
            reader: self,
            converter,
            _phantom: PhantomData,
        }
    }

    /// Transform the trajectory reader into trajectory analyzer.
    #[inline(always)]
    fn analyze<A>(self, analyzer: A) -> TrajAnalyzer<'a, Self, A>
    where
        A: FrameAnalyze,
    {
        TrajAnalyzer {
            reader: self,
            analyzer,
            _phantom: PhantomData,
        }
    }

    /// Transform the trajectory reader into trajectory converter-analyzer.
    #[inline(always)]
    fn convert_and_analyze<CA>(self, converter_analyzer: CA) -> TrajConverterAnalyzer<'a, Self, CA>
    where
        CA: FrameConvertAnalyze,
    {
        TrajConverterAnalyzer {
            reader: self,
            converter_analyzer,
            _phantom: PhantomData,
        }
    }
}

/// Blanket implementation of `ConvertableTrajRead` trait for all trajectory readers.
impl<'a, T> ConvertableTrajRead<'a> for T where T: TrajMasterRead<'a> {}
