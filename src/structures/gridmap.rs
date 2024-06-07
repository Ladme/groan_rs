// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of a higher-level utility GridMap structure for use in `groan_rs` programs.

use ndarray::Array2;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::num::Wrapping;
use std::path::Path;
use std::{fmt::Display, io::Write};

use crate::errors::GridMapError;
use crate::structures::simbox::SimBox;

/// A map of values for tiles in the xy plane.
/// Useful for analyzing membrane simulations.
///
/// The structure requires two generic parameters:
/// - `RawValue` is the type of the value that will actually be stored inside the grid map.
/// - `VisValue` is the type of the value that will be written out after calling `write_map`.
///
/// The `converter` function handles the conversion from `RawValue` to `VisValue`.
/// The function is automatically called inside `write_map`.
///
/// ## Notes
/// - The actual span of the map is not `span.0` to
///   `span.1`, but `(span.0 - tile_dim / 2.0)` to `(span.1 + tile_dim / 2.0)`.`
///
/// - Consider a span of 2 to 8 nm with a grid size of 2 nm. Such map will actually
///   capture values between 1 to 9 nm. I.e., such grid map will contain
///   four grid tiles: A (1-3 nm), B (3-5 nm), C (5-7 nm), D (7-9 nm).
/// ```text
///  ... ... ... ...
/// : A : B : C : D :
/// | |   |   |   | |
/// 1 2   4   6   8 9
///
/// ```
#[derive(Debug, Clone)]
pub struct GridMap<RawValue: Default + Clone + std::fmt::Debug, VisValue: Display> {
    /// Span of the map in the x-dimension.
    span_x: (f32, f32),
    /// Span of the map in the y-dimension.
    span_y: (f32, f32),
    /// Dimensions of a single grid tile.
    tile_dim: (f32, f32),
    /// All values of the individual grid tiles.
    values: Array2<RawValue>,
    /// Function for converting from `RawValue` to `VisValue`.
    converter: fn(&RawValue) -> VisValue,
}

impl<RawValue: Default + Clone + std::fmt::Debug, VisValue: Display> GridMap<RawValue, VisValue> {
    /// Create an grid map structure of specified span with specified size of the grid tile.
    /// The map will be filled with default values for the `RawValue`.
    ///
    /// The `converter` function handles conversion from internal data inside the map
    /// to data that should be written out in case `write_map` is called.
    ///
    /// ## Example
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// // converter function that takes average of a vector
    /// fn average(vec: &Vec<f32>) -> f32 {
    ///     let sum: f32 = vec.iter().sum();
    ///     sum / vec.len() as f32
    /// }
    ///
    /// // creates a grid map capturing a rectangular area
    /// // between -1.5 and 12.5 nm along the x-axis
    /// // and between 4.0 and 5.0 nm along the y-axis
    /// // each tile of the grid will have a size of 0.1 × 0.02 nm
    /// let map = GridMap::new(
    ///     (-1.5, 12.5), // span along the x-axis
    ///     (4.0, 5.0),   // span along the y-axis
    ///     (0.1, 0.02),  // dimensions of each grid tile (x, y)
    ///     average);     // converter function
    ///
    /// // the map vill be filled with empty vectors
    /// ```
    pub fn new(
        span_x: (f32, f32),
        span_y: (f32, f32),
        tile_dim: (f32, f32),
        converter: fn(&RawValue) -> VisValue,
    ) -> Result<GridMap<RawValue, VisValue>, GridMapError> {
        let values = Array2::<RawValue>::default((
            Self::get_len(span_x, tile_dim.0)?,
            Self::get_len(span_y, tile_dim.1)?,
        ));

        Ok(GridMap {
            span_x,
            span_y,
            tile_dim,
            values,
            converter,
        })
    }

    /// Get the width/height of the map from span and grid tile width/height.
    ///
    /// Performs sanity checks.
    fn get_len(span: (f32, f32), tile: f32) -> Result<usize, GridMapError> {
        let diff = span.1 - span.0;
        if diff < 0.0 {
            return Err(GridMapError::InvalidSpan);
        }

        if tile > diff || tile == 0.0 {
            return Err(GridMapError::InvalidGridTile);
        }

        Ok((diff / tile).round() as usize + 1)
    }

    /// Create a new grid map structure spanning the entire simulation box.
    /// The map will be filled with default values of `RawValue`.
    ///
    /// ## Warning
    /// Only supports orthogonal simulation boxes!
    pub fn from_box(
        simbox: &SimBox,
        tile_dim: (f32, f32),
        converter: fn(&RawValue) -> VisValue,
    ) -> Result<GridMap<RawValue, VisValue>, GridMapError> {
        let span_x = (0.0, simbox.x);
        let span_y = (0.0, simbox.y);

        GridMap::new(span_x, span_y, tile_dim, converter)
    }

    /// Creates a new grid map from an input file.
    ///
    /// ## Expected File Format
    /// Each line of the file should contain three segments:
    /// 1. The x-coordinate value
    /// 2. The y-coordinate value
    /// 3. The serialized `RawValue`
    ///
    /// All possible x-coordinate values must be listed before the y-coordinate changes.
    /// The spacing between coordinates must be consistent. The coordinates should
    /// be ordered and going from the lowest to the highest.
    ///
    /// ### Example:
    /// ```text
    /// 0.0 0.0 1.784562
    /// 1.0 0.0 2.943324
    /// 2.0 0.0 4.213443
    /// 0.0 1.0 0.432224
    /// 1.0 1.0 2.434545
    /// 2.0 1.0 3.947432
    /// 0.0 2.0 4.347842
    /// 1.0 2.0 0.432443
    /// 2.0 2.0 3.024943
    /// ```
    ///
    /// The file may contain comment lines, which must start with a substring specified in the `comments` list.
    /// Empty lines are skipped by default.
    ///
    /// ## Arguments
    /// - `filename` - Path to the file to read the grid map from.
    /// - `converter` - Function to convert from `RawValue` to `VisValue` in the map.
    /// - `split` - Characters used to separate the x, y, and value in each line.
    /// - `parse` - Function specifying how the `RawValue` should be parsed.
    /// - `comments` - List of strings, each of which initiates a comment line that will be skipped.
    pub fn from_file(
        filename: impl AsRef<Path>,
        converter: fn(&RawValue) -> VisValue,
        split: &[char],
        parser: fn(&str) -> Option<RawValue>,
        comments: &[&str],
    ) -> Result<GridMap<RawValue, VisValue>, GridMapError> {
        let file = File::open(&filename)
            .map_err(|_| GridMapError::FileNotFound(Box::from(filename.as_ref())))?;
        let reader = BufReader::new(file);

        let mut x_coords = Vec::new();
        let mut y_coords = Vec::new();
        let mut values = Vec::new();

        for line in reader.lines() {
            let line =
                line.map_err(|_| GridMapError::CouldNotReadLine(Box::from(filename.as_ref())))?;

            if line.trim().is_empty() {
                continue;
            }

            if comments.iter().any(|char| line.starts_with(char)) {
                continue;
            }

            let split: Vec<&str> = line.split(split).collect();
            let x = Self::parse_coordinate(&split, 0, &line, &filename)?;
            let y = Self::parse_coordinate(&split, 1, &line, &filename)?;
            let z = match split.get(2) {
                Some(x) => parser(x).ok_or(GridMapError::CouldNotParseLine(
                    line.to_owned(),
                    Box::from(filename.as_ref()),
                ))?,
                None => {
                    return Err(GridMapError::CouldNotParseLine(
                        line.to_owned(),
                        Box::from(filename.as_ref()),
                    ))
                }
            };

            values.push(z);

            if !x_coords.contains(&x) {
                x_coords.push(x);
            }
            if !y_coords.contains(&y) {
                y_coords.push(y);
            }
        }

        if x_coords.len() < 2 || y_coords.len() < 2 {
            return Err(GridMapError::EmptyGridMap);
        }

        let (span_x, tile_x) = Self::properties_from_coords(&x_coords)?;
        let (span_y, tile_y) = Self::properties_from_coords(&y_coords)?;

        let mut map = GridMap::from_vec(span_y, span_x, (tile_y, tile_x), values, converter)?;
        map.transpose();
        Ok(map)
    }

    /// Get map span and grid tile height/width from coordinates.
    fn properties_from_coords(coords: &[f32]) -> Result<((f32, f32), f32), GridMapError> {
        let min = *coords.first().expect("FATAL GROAN ERROR | GridMap::properties_from_coords | First coordinate value should exist.");
        let second = *coords.get(1).expect("FATAL GROAN ERROR | GridMap::properties_from_coords | Second coordinate value should exist.");
        let max = *coords.last().expect("FATAL GROAN ERROR | GridMap::properties_from_coords | Last coordinate value should exist.");

        let tile = second - min;

        Ok(((min, max), tile))
    }

    /// Creates a new grid map from grid map span, grid tile dimensions, and a vector of values.
    /// The values are mapped in the same way as in `Array2::from_shape_vec`.
    ///
    /// The `converter` function handles conversion from internal data inside the map
    /// to data that should be written out in case `write_map` is called.
    ///
    /// Returns an error in case the provided number of values does not match the expected number of values.
    pub fn from_vec(
        span_x: (f32, f32),
        span_y: (f32, f32),
        tile_dim: (f32, f32),
        values: Vec<RawValue>,
        converter: fn(&RawValue) -> VisValue,
    ) -> Result<GridMap<RawValue, VisValue>, GridMapError> {
        let len_x = Self::get_len(span_x, tile_dim.0)?;
        let len_y = Self::get_len(span_y, tile_dim.1)?;
        let n_values = values.len();
        if len_x * len_y != n_values {
            return Err(GridMapError::InvalidMapDimensions(len_x * len_y, n_values));
        }

        let values = Array2::<RawValue>::from_shape_vec((len_x, len_y), values).expect(
            "FATAL GROAN ERROR | GridMap::from_vec | Could not construct 2D array of values.",
        );

        Ok(GridMap {
            span_x,
            span_y,
            tile_dim,
            values,
            converter,
        })
    }

    /// Get the total number of tiles in the grid map.
    #[inline(always)]
    pub fn n_tiles(&self) -> usize {
        self.values.len()
    }

    /// Get the number of tiles along the x-dimension.
    #[inline(always)]
    pub fn n_tiles_x(&self) -> usize {
        self.values.nrows()
    }

    /// Get the number of tiles along the y-dimension.
    #[inline(always)]
    pub fn n_tiles_y(&self) -> usize {
        self.values.ncols()
    }

    /// Get the unconverted value of the map at target coordinates.
    ///
    /// Returns `None` in case the coordinates are out of the span of the map.
    #[inline(always)]
    pub fn get_at(&self, x: f32, y: f32) -> Option<&RawValue> {
        self.values.get((self.x2index(x), self.y2index(y)))
    }

    /// Get mutable reference to the unconverted value of the map at target coordinates.
    ///
    /// Returns `None` in case the coordinates are out of the span of the map.
    #[inline(always)]
    pub fn get_mut_at(&mut self, x: f32, y: f32) -> Option<&mut RawValue> {
        self.values.get_mut((self.x2index(x), self.y2index(y)))
    }

    /// Get the value of the map at target coordinates
    /// and convert this value to `VisValue`.
    /// The original value inside the map is not modified.
    ///
    /// Returns `None` in case the coordinates are out of the span of the map.
    #[inline(always)]
    pub fn get_at_convert(&self, x: f32, y: f32) -> Option<VisValue> {
        Some((self.converter)(self.get_at(x, y)?))
    }

    /// Write the map with RAW UNCONVERTED values into the provided writer structure.
    pub fn write_map_raw(&self, writer: &mut impl Write) -> Result<(), GridMapError> {
        self.extract_raw().try_for_each(|x| {
            writeln!(writer, "{:10.6} {:10.6} {:?}", x.0, x.1, x.2)
                .map_err(|_| GridMapError::CouldNotWrite)
        })?;

        Ok(())
    }

    /// Write the map into the provided writer structure. Convert the values of the map.
    pub fn write_map(&self, writer: &mut impl Write) -> Result<(), GridMapError> {
        self.extract_convert().try_for_each(|x| {
            writeln!(writer, "{:10.6} {:10.6} {}", x.0, x.1, x.2)
                .map_err(|_| GridMapError::CouldNotWrite)
        })?;

        Ok(())
    }

    /// Extract the map into an iterator over (f32, f32, &RawValue).
    ///
    /// Each element of the returned iterator corresponds to
    /// a grid tile of the map with its associated value.
    /// In other words, each element is a tuple of x-coordinate, y-coordinate and the value.
    ///
    /// ## Notes
    /// - The iterator first loops through all x-coordinates before the y-coordinate is changed.
    pub fn extract_raw(&self) -> impl Iterator<Item = (f32, f32, &RawValue)> + '_ {
        (0..self.n_tiles_y())
            .flat_map(move |y| {
                let y_coord = self.index2y(y);
                (0..self.n_tiles_x()).map(move |x| {
                    let x_coord = self.index2x(x);
                    (x_coord, y_coord, self.get_at(x_coord, y_coord)
                        .expect("FATAL GROAN ERROR | GridMap::extract_raw | Value at target coordinates must exist."))
                })
            })
    }

    /// Extract the map into an iterator over (f32, f32, VisValue).
    ///
    /// Each element of the returned iterator corresponds to a grid tile
    /// of the map with its associated CONVERTED value.
    /// In other words, each element is a tuple of x-coordinate, y-coordinate and the converted value.
    ///
    /// ## Notes
    /// - The iterator first loops through all x-coordinates before the y-coordinate is changed.
    pub fn extract_convert(&self) -> impl Iterator<Item = (f32, f32, VisValue)> + '_ {
        self.extract_raw()
            .map(|x| (x.0, x.1, (self.converter)(x.2)))
    }

    /// Transpose (rotate) the map. Used for reading gridmap files.
    fn transpose(&mut self) {
        (self.span_x, self.span_y) = (self.span_y, self.span_x);
        (self.tile_dim.0, self.tile_dim.1) = (self.tile_dim.1, self.tile_dim.0);
        self.values = self.values.t().to_owned();
    }

    /// Convert an x-coordinate to an index in the map.
    /// Always returns an index even if the coordinate is outside the span of the map.
    #[inline(always)]
    fn x2index(&self, x: f32) -> usize {
        Wrapping(((x - self.span_x.0) / self.tile_dim.0).round() as isize).0 as usize
    }

    /// Convert a y-coordinate to an index in the map.
    /// Always returns an index even if the coordinate is outside the span of the map.
    #[inline(always)]
    fn y2index(&self, y: f32) -> usize {
        Wrapping(((y - self.span_y.0) / self.tile_dim.1).round() as isize).0 as usize
    }

    /// Convert an index of the map to an x coordinate.
    /// Always returns a coordinate even if the index is outside the span of the map.
    #[inline(always)]
    fn index2x(&self, index: usize) -> f32 {
        (index as f32 * self.tile_dim.0) + self.span_x.0
    }

    /// Convert an index of the map to a y coordinate.
    /// Always returns a coordinate even if the index is outside the span of the map.
    #[inline(always)]
    fn index2y(&self, index: usize) -> f32 {
        (index as f32 * self.tile_dim.1) + self.span_y.0
    }

    /// Parse an x- or y-coordinate.
    fn parse_coordinate(
        split: &[&str],
        index: usize,
        line: &str,
        filename: impl AsRef<Path>,
    ) -> Result<f32, GridMapError> {
        split
            .get(index)
            .ok_or_else(|| {
                GridMapError::CouldNotParseLine(line.to_owned(), Box::from(filename.as_ref()))
            })
            .and_then(|s| {
                s.trim().parse::<f32>().map_err(|_| {
                    GridMapError::CouldNotParseLine(line.to_owned(), Box::from(filename.as_ref()))
                })
            })
    }
}

/******************************/
/*        UNIT TESTS          */
/******************************/

#[cfg(test)]
mod tests {
    use std::io::BufWriter;

    use float_cmp::assert_approx_eq;

    use super::*;

    fn sum(vec: &Vec<f32>) -> f32 {
        vec.iter().sum()
    }

    #[test]
    fn new() {
        let gridmap =
            GridMap::<Vec<f32>, f32>::new((-2.0, 7.0), (3.0, 6.0), (0.15, 0.20), sum).unwrap();

        assert_eq!(gridmap.span_x, (-2.0, 7.0));
        assert_eq!(gridmap.span_y, (3.0, 6.0));
        assert_eq!(gridmap.tile_dim, (0.15, 0.20));
        assert_eq!(gridmap.n_tiles(), 61 * 16);
        assert_eq!(gridmap.n_tiles_x(), 61);
        assert_eq!(gridmap.n_tiles_y(), 16);
        assert_eq!(gridmap.values[(5, 8)], vec![]);
    }

    #[test]
    fn new_failures() {
        assert!(matches!(
            GridMap::<Vec<f32>, f32>::new((2.0, -2.0), (3.0, 6.0), (0.15, 0.20), sum),
            Err(GridMapError::InvalidSpan)
        ));
        assert!(matches!(
            GridMap::<Vec<f32>, f32>::new((-2.0, 7.0), (3.0, -6.0), (0.15, 0.20), sum),
            Err(GridMapError::InvalidSpan)
        ));
        assert!(matches!(
            GridMap::<Vec<f32>, f32>::new((-2.0, 7.0), (3.0, 6.0), (0.0, 0.20), sum),
            Err(GridMapError::InvalidGridTile)
        ));
        assert!(matches!(
            GridMap::<Vec<f32>, f32>::new((-2.0, 7.0), (3.0, 6.0), (0.15, 0.0), sum),
            Err(GridMapError::InvalidGridTile)
        ));
        assert!(matches!(
            GridMap::<Vec<f32>, f32>::new((-2.0, 7.0), (3.0, 6.0), (10.5, 0.20), sum),
            Err(GridMapError::InvalidGridTile)
        ));
        assert!(matches!(
            GridMap::<Vec<f32>, f32>::new((-2.0, 7.0), (3.0, 6.0), (0.15, 3.10), sum),
            Err(GridMapError::InvalidGridTile)
        ));
    }

    #[test]
    fn from_box() {
        let simbox = SimBox::from([10.0, 20.0, 14.0]);
        let gridmap = GridMap::<Vec<f32>, f32>::from_box(&simbox, (0.15, 0.20), sum).unwrap();

        assert_eq!(gridmap.span_x, (0.0, 10.0));
        assert_eq!(gridmap.span_y, (0.0, 20.0));
        assert_eq!(gridmap.tile_dim, (0.15, 0.20));
        assert_eq!(gridmap.n_tiles(), 68 * 101);
        assert_eq!(gridmap.n_tiles_x(), 68);
        assert_eq!(gridmap.n_tiles_y(), 101);
        assert_eq!(gridmap.values[(5, 8)], vec![]);
    }

    #[test]
    fn from_vec() {
        let values = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let gridmap = GridMap::<usize, usize>::from_vec(
            (1.0, 2.0),
            (1.0, 2.5),
            (1.0, 0.5),
            values,
            usize::to_owned,
        )
        .unwrap();

        assert_eq!(gridmap.n_tiles_x(), 2);
        assert_eq!(gridmap.n_tiles_y(), 4);
        assert_eq!(gridmap.n_tiles(), 8);

        assert_eq!(*gridmap.get_at(1.0, 1.0).unwrap(), 1);
        assert_eq!(*gridmap.get_at(2.0, 1.0).unwrap(), 5);
        assert_eq!(*gridmap.get_at(1.0, 1.5).unwrap(), 2);
        assert_eq!(*gridmap.get_at(2.0, 1.5).unwrap(), 6);
        assert_eq!(*gridmap.get_at(1.0, 2.0).unwrap(), 3);
        assert_eq!(*gridmap.get_at(2.0, 2.0).unwrap(), 7);
        assert_eq!(*gridmap.get_at(1.0, 2.5).unwrap(), 4);
        assert_eq!(*gridmap.get_at(2.0, 2.5).unwrap(), 8);
    }

    fn sum_usize(vec: &Vec<usize>) -> usize {
        vec.iter().sum()
    }

    fn parse_vec(string: &str) -> Option<Vec<usize>> {
        let split: Vec<&str> = string.split_whitespace().collect();
        Some(
            split
                .into_iter()
                .map(|x| x.parse::<usize>().unwrap())
                .collect::<Vec<usize>>(),
        )
    }

    #[test]
    fn from_file() {
        let gridmap = GridMap::<Vec<usize>, usize>::from_file(
            "test_files/gridmaps/map.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        )
        .unwrap();

        assert_eq!(*gridmap.get_at(0.0, 0.0).unwrap(), vec![10]);
        assert_eq!(*gridmap.get_at(1.0, 0.0).unwrap(), vec![5, 4]);
        assert_eq!(*gridmap.get_at(2.0, 0.0).unwrap(), vec![43, 23, 21]);
        assert_eq!(*gridmap.get_at(0.0, 1.0).unwrap(), vec![4, 8, 12]);
        assert_eq!(*gridmap.get_at(1.0, 1.0).unwrap(), vec![2, 5]);
        assert_eq!(*gridmap.get_at(2.0, 1.0).unwrap(), vec![]);
    }

    #[test]
    fn coord2index() {
        let gridmap =
            GridMap::<Vec<f32>, f32>::new((-2.0, 7.0), (3.0, 6.0), (0.15, 0.20), sum).unwrap();

        assert_eq!(gridmap.x2index(-2.0), 0);
        assert_approx_eq!(f32, gridmap.index2x(0), -2.0);

        assert_eq!(gridmap.x2index(-2.07), 0);
        assert_eq!(gridmap.x2index(-2.08), usize::MAX);
        assert_eq!(gridmap.x2index(-2.32), usize::MAX - 1);

        assert_eq!(gridmap.x2index(7.0), 60);
        assert_approx_eq!(f32, gridmap.index2x(60), 7.0);

        assert_eq!(gridmap.x2index(7.07), 60);
        assert_eq!(gridmap.x2index(7.08), 61);
        assert_eq!(gridmap.x2index(7.32), 62);

        assert_eq!(gridmap.y2index(3.0), 0);
        assert_approx_eq!(f32, gridmap.index2y(0), 3.0);

        assert_eq!(gridmap.y2index(2.91), 0);
        assert_eq!(gridmap.y2index(2.89), usize::MAX);
        assert_eq!(gridmap.y2index(2.58), usize::MAX - 1);

        assert_eq!(gridmap.y2index(6.0), 15);
        assert_approx_eq!(f32, gridmap.index2y(15), 6.0);

        assert_eq!(gridmap.y2index(6.09), 15);
        assert_eq!(gridmap.y2index(6.11), 16);
        assert_eq!(gridmap.y2index(6.42), 17);
    }

    #[test]
    fn get_at() {
        let mut gridmap =
            GridMap::<Vec<f32>, f32>::new((-2.0, 7.0), (3.0, 6.0), (0.15, 0.20), sum).unwrap();

        let val = gridmap.get_mut_at(0.45, 3.57).unwrap();
        val.push(10.4);
        val.push(2.8);

        let val = gridmap.get_at(0.45, 3.57).unwrap();
        assert_eq!(val.len(), 2);
        assert_approx_eq!(f32, val[0], 10.4);
        assert_approx_eq!(f32, val[1], 2.8);

        assert_approx_eq!(f32, gridmap.get_at_convert(0.39, 3.52).unwrap(), 13.2);

        // nearby tiles are empty
        assert_eq!(gridmap.get_at(0.56, 3.57).unwrap().len(), 0);
        assert_eq!(gridmap.get_at(0.45, 3.42).unwrap().len(), 0);
    }

    #[test]
    fn extract_raw() {
        let gridmap = GridMap::<Vec<usize>, usize>::from_file(
            "test_files/gridmaps/map.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        )
        .unwrap();

        let vec = gridmap
            .extract_raw()
            .collect::<Vec<(f32, f32, &Vec<usize>)>>();

        assert_eq!(vec[0], (0.0, 0.0, &vec![10]));
        assert_eq!(vec[1], (1.0, 0.0, &vec![5, 4]));
        assert_eq!(vec[2], (2.0, 0.0, &vec![43, 23, 21]));
        assert_eq!(vec[3], (0.0, 1.0, &vec![4, 8, 12]));
        assert_eq!(vec[4], (1.0, 1.0, &vec![2, 5]));
        assert_eq!(vec[5], (2.0, 1.0, &vec![]));
    }

    #[test]
    fn extract_convert() {
        let gridmap = GridMap::<Vec<usize>, usize>::from_file(
            "test_files/gridmaps/map.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        )
        .unwrap();

        let vec = gridmap
            .extract_convert()
            .collect::<Vec<(f32, f32, usize)>>();

        assert_eq!(vec[0], (0.0, 0.0, 10));
        assert_eq!(vec[1], (1.0, 0.0, 9));
        assert_eq!(vec[2], (2.0, 0.0, 87));
        assert_eq!(vec[3], (0.0, 1.0, 24));
        assert_eq!(vec[4], (1.0, 1.0, 7));
        assert_eq!(vec[5], (2.0, 1.0, 0));
    }

    #[test]
    fn write_raw() {
        let gridmap = GridMap::<Vec<usize>, usize>::from_file(
            "test_files/gridmaps/map.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        )
        .unwrap();

        let mut output = BufWriter::new(Vec::new());
        gridmap.write_map_raw(&mut output).unwrap();
        let output_string = String::from_utf8(output.into_inner().unwrap()).unwrap();

        let expected = "  0.000000   0.000000 [10]\n  1.000000   0.000000 [5, 4]\n  2.000000   0.000000 [43, 23, 21]\n  0.000000   1.000000 [4, 8, 12]\n  1.000000   1.000000 [2, 5]\n  2.000000   1.000000 []\n";

        assert_eq!(output_string, expected);
    }

    #[test]
    fn write() {
        let gridmap = GridMap::<Vec<usize>, usize>::from_file(
            "test_files/gridmaps/map.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        )
        .unwrap();

        let mut output = BufWriter::new(Vec::new());
        gridmap.write_map(&mut output).unwrap();
        let output_string = String::from_utf8(output.into_inner().unwrap()).unwrap();

        let expected = "  0.000000   0.000000 10\n  1.000000   0.000000 9\n  2.000000   0.000000 87\n  0.000000   1.000000 24\n  1.000000   1.000000 7\n  2.000000   1.000000 0\n";

        assert_eq!(output_string, expected);
    }
}
