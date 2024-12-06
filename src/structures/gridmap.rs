// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Implementation of a higher-level utility GridMap structure for use in `groan_rs` programs.

use float_cmp::approx_eq;
use getset::CopyGetters;
use ndarray::{Array2, ShapeBuilder};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::num::Wrapping;
use std::path::Path;
use std::{fmt::Display, io::Write};

use crate::errors::GridMapError;
use crate::structures::simbox::SimBox;

/// Describes whether the data is written or should be written
/// in row-major order or column-major order.
/// Since Rust uses row-major ordering, using `RowMajor` is generally advised.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum DataOrder {
    #[default]
    RowMajor,
    ColumnMajor,
}

/// A map of values for tiles in the xy plane.
/// Useful for analyzing membrane simulations.
///
/// The structure requires three generic parameters:
/// - `RawValue` is the type of the value that will actually be stored inside the grid map.
/// - `VisValue` is the type of the value that will be written out after calling `write_map`.
/// - `Converter` is a function/closure that handles the conversion from `RawValue` to `VisValue`.
///    The function is automatically called inside `write_map`, `get_at_convert`, and `extract_convert`.
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
#[derive(Debug, Clone, CopyGetters)]
pub struct GridMap<
    RawValue: Default + Clone + std::fmt::Debug,
    VisValue: Display,
    Converter: Fn(&RawValue) -> VisValue,
> {
    /// Span of the map in the x-dimension.
    #[getset(get_copy = "pub")]
    span_x: (f32, f32),
    /// Span of the map in the y-dimension.
    #[getset(get_copy = "pub")]
    span_y: (f32, f32),
    /// Dimensions of a single grid tile.
    #[getset(get_copy = "pub")]
    tile_dim: (f32, f32),
    /// All values of the individual grid tiles.
    values: Array2<RawValue>,
    /// Function or closure for converting from `RawValue` to `VisValue`.
    converter: Converter,
}

impl<
        RawValue: Default + Clone + std::fmt::Debug,
        VisValue: Display,
        Converter: Fn(&RawValue) -> VisValue,
    > GridMap<RawValue, VisValue, Converter>
{
    /// Create a grid map structure of specified span with specified size of the grid tile.
    /// The map will be filled with default values for the `RawValue`.
    ///
    /// The `converter` function handles conversion from internal data inside the map
    /// to data that should be written out in case `write_map` is called.
    ///
    /// ## Examples
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
    /// // the map will be filled with empty vectors
    /// ```
    ///
    /// Converter function can also be a closure:
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// let mut value = 10.0;
    /// let sum_plus_value = move |raw: &Vec<f32>| {
    ///     raw.iter().sum::<f32>() + value
    /// };
    ///
    /// let map = GridMap::new(
    ///     (-1.5, 12.5),
    ///     (4.0, 5.0),
    ///     (0.1, 0.02),
    ///     sum_plus_value);
    /// ```
    pub fn new(
        span_x: (f32, f32),
        span_y: (f32, f32),
        tile_dim: (f32, f32),
        converter: Converter,
    ) -> Result<GridMap<RawValue, VisValue, Converter>, GridMapError> {
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

    /// Get the number of tiles in a particular dimension
    /// of the map from span and grid tile width/height.
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
        converter: Converter,
    ) -> Result<GridMap<RawValue, VisValue, Converter>, GridMapError> {
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
    /// Elements can be provided in row-major or in column-major order, the function can recognize the order used.
    /// The spacing between coordinates must be consistent. The coordinates must
    /// be ordered and going from the lowest to the highest.
    ///
    /// ### Examples
    /// Row-major order: (faster)
    /// ```text
    /// 0.0 0.0 1.784562
    /// 0.0 1.0 0.432224
    /// 0.0 2.0 4.347842
    /// 1.0 0.0 2.943324
    /// 1.0 1.0 2.434545
    /// 1.0 2.0 0.432443
    /// 2.0 0.0 4.213443
    /// 2.0 1.0 3.947432
    /// 2.0 2.0 3.024943
    /// ```
    ///
    /// Column-major order:
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
        converter: Converter,
        split: &[char],
        parser: impl Fn(&str) -> Option<RawValue>,
        comments: &[&str],
    ) -> Result<GridMap<RawValue, VisValue, Converter>, GridMapError> {
        let file = File::open(&filename)
            .map_err(|_| GridMapError::FileNotFound(Box::from(filename.as_ref())))?;
        let reader = BufReader::new(file);

        let mut x_coords = Vec::new();
        let mut y_coords = Vec::new();
        let mut values = Vec::new();

        let mut data_order = DataOrder::RowMajor;
        for line in reader.lines() {
            let line =
                line.map_err(|_| GridMapError::CouldNotReadLine(Box::from(filename.as_ref())))?;

            if line.trim().is_empty() {
                continue;
            }

            if comments.iter().any(|char| line.starts_with(char)) {
                continue;
            }

            // split the line and remove all empty substrings
            let split: Vec<&str> = line
                .split(split)
                .filter(|substring| !substring.trim().is_empty())
                .collect();

            let x = Self::parse_coordinate(&split, 0, &line, &filename)?;
            let y = Self::parse_coordinate(&split, 1, &line, &filename)?;
            let z = match split.get(2) {
                Some(x) => parser(x).ok_or(GridMapError::CouldNotParseLine(
                    line.to_owned(),
                    Box::from(filename.as_ref()),
                ))?,
                None => parser(" ").ok_or(GridMapError::CouldNotParseLine(
                    line.to_owned(),
                    Box::from(filename.as_ref()),
                ))?,
            };

            values.push(z);

            x_coords.push(x);
            y_coords.push(y);

            // check for column-major order
            if y_coords.len() == 2 && y_coords[0] == y_coords[1] {
                data_order = DataOrder::ColumnMajor;
            }
        }

        if x_coords.len() < 2 || y_coords.len() < 2 {
            return Err(GridMapError::EmptyGridMap);
        }

        let (span_x, span_y, (tile_x, tile_y)) =
            Self::validate_and_get_span(&x_coords, &y_coords, data_order)?;

        GridMap::from_vec(
            span_x,
            span_y,
            (tile_x, tile_y),
            values,
            data_order,
            converter,
        )
    }

    /// Check that the x and y coordinates used for the construction of the GridMap are valid.
    /// Also get the span of the grid map.
    #[allow(clippy::type_complexity)]
    fn validate_and_get_span(
        x_coords: &[f32],
        y_coords: &[f32],
        data_order: DataOrder,
    ) -> Result<((f32, f32), (f32, f32), (f32, f32)), GridMapError> {
        // choose which coordinates are fast or slow changing based on data order
        let (fast_changing_coords, slow_changing_coords) = match data_order {
            DataOrder::ColumnMajor => (x_coords, y_coords),
            DataOrder::RowMajor => (y_coords, x_coords),
        };

        // validate the slow-changing coordinates
        let (min_s, second_s, max_s, block_size) =
            Self::validate_slow_changing_coords(slow_changing_coords)?;

        // validate the fast-changing coordinates
        let (min_f, second_f, max_f) =
            Self::validate_fast_changing_coords(fast_changing_coords, block_size)?;

        // check that the order of fast-changing coordinates is increasing
        Self::check_order_increasing(
            fast_changing_coords,
            slow_changing_coords,
            data_order,
            block_size,
        )?;

        // calculate the tile sizes
        let tile_f = second_f - min_f;
        let tile_s = second_s - min_s;

        // return the appropriate spans based on the data order
        match data_order {
            DataOrder::ColumnMajor => Ok(((min_f, max_f), (min_s, max_s), (tile_f, tile_s))),
            DataOrder::RowMajor => Ok(((min_s, max_s), (min_f, max_f), (tile_s, tile_f))),
        }
    }

    /// Validate slow-changing coordinates, calculate their span, and determine block size.
    fn validate_slow_changing_coords(
        slow_changing_coords: &[f32],
    ) -> Result<(f32, f32, f32, usize), GridMapError> {
        let mut previous_slow = *slow_changing_coords.first().expect(
        "FATAL GROAN ERROR | GridMap::validate_slow_changing_coords | First coordinate value should exist.",
    );
        let mut counter = 0;
        let mut block_size = None;
        let (mut min_s, mut second_s, mut max_s) = (f32::MAX, f32::MAX, f32::MIN);

        Self::update_min_max(previous_slow, &mut min_s, &mut second_s, &mut max_s);

        for &s in slow_changing_coords.iter().skip(1) {
            counter += 1;

            if !approx_eq!(f32, s, previous_slow) {
                if previous_slow > s {
                    return Err(GridMapError::NotIncreasing(
                        s.to_string(),
                        previous_slow.to_string(),
                    ));
                }

                match block_size {
                    Some(block) => match counter.cmp(&block) {
                        std::cmp::Ordering::Less => {
                            return Err(GridMapError::InvalidCoordinates(s.to_string()));
                        }
                        std::cmp::Ordering::Greater => {
                            panic!("FATAL GROAN ERROR | GridMap::validate_slow_changing_coords | Invalid slow-changing coordinate detected but this should have been handled elsewhere.");
                        }
                        std::cmp::Ordering::Equal => {
                            previous_slow = s;
                            Self::update_min_max(
                                previous_slow,
                                &mut min_s,
                                &mut second_s,
                                &mut max_s,
                            );
                        }
                    },
                    None => {
                        block_size = Some(counter);
                        previous_slow = s;
                        Self::update_min_max(previous_slow, &mut min_s, &mut second_s, &mut max_s);
                    }
                }

                counter = 0;
                continue;
            }

            if block_size.is_some() && counter >= block_size.unwrap() {
                return Err(GridMapError::InvalidCoordinates(s.to_string()));
            }
        }

        let unwrapped_block_size = block_size
            .ok_or_else(|| GridMapError::InvalidCoordinates(previous_slow.to_string()))?;

        Ok((min_s, second_s, max_s, unwrapped_block_size))
    }

    /// Validate fast-changing coordinates and calculate their span.
    fn validate_fast_changing_coords(
        fast_changing_coords: &[f32],
        block_size: usize,
    ) -> Result<(f32, f32, f32), GridMapError> {
        let (mut min_f, mut second_f, mut max_f) = (f32::MAX, f32::MAX, f32::MIN);

        for i in 0..block_size {
            let coord_i = *fast_changing_coords.get(i).expect(
            "FATAL GROAN ERROR | GridMap::validate_fast_changing_coords | Access to fast-changing coordinates is out of range. (1)",
        );

            Self::update_min_max(coord_i, &mut min_f, &mut second_f, &mut max_f);

            let mut j = i + block_size;
            while j < fast_changing_coords.len() {
                let coord_j = fast_changing_coords.get(j).expect(
                "FATAL GROAN ERROR | GridMap::validate_fast_changing_coords | Access to fast-changing coordinates is out of range. (2)",
            );

                if !approx_eq!(f32, coord_i, *coord_j) {
                    return Err(GridMapError::InvalidCoordinates(coord_j.to_string()));
                }

                j += block_size;
            }
        }

        Ok((min_f, second_f, max_f))
    }

    /// Check if fast-changing coordinates are in increasing order.
    fn check_order_increasing(
        fast_changing_coords: &[f32],
        slow_changing_coords: &[f32],
        data_order: DataOrder,
        block_size: usize,
    ) -> Result<(), GridMapError> {
        for i in 0..(block_size - 1) {
            let coord1 = fast_changing_coords.get(i).unwrap();
            let coord2 = fast_changing_coords.get(i + 1).unwrap();

            if coord1 > coord2 {
                return Err(GridMapError::NotIncreasing(
                    coord2.to_string(),
                    coord1.to_string(),
                ));
            }

            if approx_eq!(f32, *coord1, *coord2) {
                let slow_changing = slow_changing_coords.get(i).expect(
                "FATAL GROAN ERROR | GridMap::check_order_increasing | Could not find slow-changing coordinate which should exist.",
            );
                let (x, y) = match data_order {
                    DataOrder::ColumnMajor => (coord1, slow_changing),
                    DataOrder::RowMajor => (slow_changing, coord1),
                };

                return Err(GridMapError::PointDefinedMultipleTimes(
                    x.to_string(),
                    y.to_string(),
                ));
            }
        }

        Ok(())
    }

    /// Update the min, second_min, and max values.
    fn update_min_max(current: f32, min: &mut f32, second_min: &mut f32, max: &mut f32) {
        if current < *min {
            *second_min = *min;
            *min = current;
        } else if current < *second_min {
            *second_min = current;
        }

        if current > *max {
            *max = current;
        }
    }

    /// Creates a new grid map from grid map span, grid tile dimensions, and a vector of values.
    /// `data_order` specifies whether the elements provided in the `values` are row-major or column-major ordered.
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
        data_order: DataOrder,
        converter: Converter,
    ) -> Result<GridMap<RawValue, VisValue, Converter>, GridMapError> {
        let len_x = Self::get_len(span_x, tile_dim.0)?;
        let len_y = Self::get_len(span_y, tile_dim.1)?;
        let n_values = values.len();
        if len_x * len_y != n_values {
            return Err(GridMapError::InvalidMapDimensions(len_x * len_y, n_values));
        }

        let values = match data_order {
            DataOrder::RowMajor => Array2::<RawValue>::from_shape_vec((len_x, len_y), values),
            DataOrder::ColumnMajor => {
                Array2::<RawValue>::from_shape_vec((len_x, len_y).f(), values)
            }
        }
        .expect("FATAL GROAN ERROR | GridMap::from_vec | Could not construct 2D array of values.");

        Ok(GridMap {
            span_x,
            span_y,
            tile_dim,
            values,
            converter,
        })
    }

    /// Set all raw values in the gridmap to the default values.
    pub fn clear(&mut self) {
        self.values
            .iter_mut()
            .for_each(|x| *x = RawValue::default())
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

    /// Get the coordinates of tile target point is located in.
    ///
    /// Returns `None` if the point is outside of the map.
    #[inline(always)]
    pub fn get_tile(&self, x: f32, y: f32) -> Option<(f32, f32)> {
        if !self.is_inside(x, y) {
            return None;
        }

        Some((self.index2x(self.x2index(x)), self.index2y(self.y2index(y))))
    }

    /// Check whether the point is inside the grid map.
    ///
    /// - Note that the actual span of the map is not `span.0` to
    ///   `span.1`, but `(span.0 - tile_dim / 2.0)` to `(span.1 + tile_dim / 2.0)`.`
    #[inline(always)]
    pub fn is_inside(&self, x: f32, y: f32) -> bool {
        self.get_at(x, y).is_some()
    }

    /// Get the unconverted value of the map at target coordinates.
    ///
    /// Returns `None` in case the coordinates are outside of the span of the map.
    #[inline(always)]
    pub fn get_at(&self, x: f32, y: f32) -> Option<&RawValue> {
        self.values.get((self.x2index(x), self.y2index(y)))
    }

    /// Get mutable reference to the unconverted value of the map at target coordinates.
    ///
    /// Returns `None` in case the coordinates are outside of the span of the map.
    #[inline(always)]
    pub fn get_mut_at(&mut self, x: f32, y: f32) -> Option<&mut RawValue> {
        self.values.get_mut((self.x2index(x), self.y2index(y)))
    }

    /// Get the value of the map at target coordinates
    /// and convert this value to `VisValue`.
    /// The original value inside the map is not modified.
    ///
    /// Returns `None` in case the coordinates are outside of the span of the map.
    #[inline(always)]
    pub fn get_at_convert(&self, x: f32, y: f32) -> Option<VisValue> {
        Some((self.converter)(self.get_at(x, y)?))
    }

    /// Write the map with RAW UNCONVERTED values into the provided writer structure.
    ///
    /// The map is written in row-major order. For column-major order, see [`write_map_raw_column_major`](`GridMap::write_map_raw_column_major`).
    pub fn write_map_raw(&self, writer: &mut impl Write) -> Result<(), GridMapError> {
        self.extract_raw().try_for_each(|x| {
            writeln!(writer, "{:10.6} {:10.6} {:?}", x.0, x.1, x.2)
                .map_err(|_| GridMapError::CouldNotWrite)
        })?;

        Ok(())
    }

    /// Write the map with RAW UNCONVERTED values into the provided writer structure.
    ///
    /// The map is written in column-major order. For raw-major order, see [`write_map_raw`](`GridMap::write_map_raw`).
    pub fn write_map_raw_column_major(&self, writer: &mut impl Write) -> Result<(), GridMapError> {
        self.extract_raw_column_major().try_for_each(|x| {
            writeln!(writer, "{:10.6} {:10.6} {:?}", x.0, x.1, x.2)
                .map_err(|_| GridMapError::CouldNotWrite)
        })?;

        Ok(())
    }

    /// Write the map into the provided writer structure. Convert the values of the map.
    ///
    /// The map is written in row-major order. For column-major order, see [`write_map_column_major`](`GridMap::write_map_column_major`).
    pub fn write_map(&self, writer: &mut impl Write) -> Result<(), GridMapError> {
        self.extract_convert().try_for_each(|x| {
            writeln!(writer, "{:10.6} {:10.6} {}", x.0, x.1, x.2)
                .map_err(|_| GridMapError::CouldNotWrite)
        })?;

        Ok(())
    }

    /// Write the map into the provided writer structure. Convert the values of the map.
    ///
    /// The map is written in column-major order. For raw-major order, see [`write_map`](`GridMap::write_map`).
    pub fn write_map_column_major(&self, writer: &mut impl Write) -> Result<(), GridMapError> {
        self.extract_convert_column_major().try_for_each(|x| {
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
    /// The data are exported in row-major order. For column-major order, see [`extract_raw_column_major`](`GridMap::extract_raw_column_major`).
    pub fn extract_raw(&self) -> impl Iterator<Item = (f32, f32, &RawValue)> + '_ {
        (0..self.n_tiles_x())
            .flat_map(move |x| {
                let x_coord = self.index2x(x);
                (0..self.n_tiles_y()).map(move |y| {
                    let y_coord = self.index2y(y);
                    (x_coord, y_coord, self.get_at(x_coord, y_coord)
                        .expect("FATAL GROAN ERROR | GridMap::extract_raw | Value at target coordinates must exist."))
                })
            }
        )
    }

    /// Extract the map into an iterator over (f32, f32, &RawValue).
    ///
    /// Each element of the returned iterator corresponds to
    /// a grid tile of the map with its associated value.
    /// In other words, each element is a tuple of x-coordinate, y-coordinate and the value.
    ///
    /// The data are exported in column-major order. For row-major order, see [`extract_raw`](`GridMap::extract_raw`).
    pub fn extract_raw_column_major(&self) -> impl Iterator<Item = (f32, f32, &RawValue)> + '_ {
        (0..self.n_tiles_y())
            .flat_map(move |y| {
                let y_coord = self.index2y(y);
                (0..self.n_tiles_x()).map(move |x| {
                    let x_coord = self.index2y(x);
                    (x_coord, y_coord, self.get_at(x_coord, y_coord)
                        .expect("FATAL GROAN ERROR | GridMap::extract_raw_column_major | Value at target coordinates must exist."))
                })
            }
        )
    }

    /// Extract the map into an iterator over (f32, f32, VisValue).
    ///
    /// Each element of the returned iterator corresponds to a grid tile
    /// of the map with its associated CONVERTED value.
    /// In other words, each element is a tuple of x-coordinate, y-coordinate and the converted value.
    ///
    /// The data are exported in row-major order.
    /// For column-major order, see [`extract_convert_column_major`](`GridMap::extract_convert_column_major`).
    pub fn extract_convert(&self) -> impl Iterator<Item = (f32, f32, VisValue)> + '_ {
        self.extract_raw()
            .map(|x| (x.0, x.1, (self.converter)(x.2)))
    }

    /// Extract the map into an iterator over (f32, f32, VisValue).
    ///
    /// Each element of the returned iterator corresponds to a grid tile
    /// of the map with its associated CONVERTED value.
    /// In other words, each element is a tuple of x-coordinate, y-coordinate and the converted value.
    ///
    /// The data are exported in column-major order.
    /// For row-major order, see [`extract_convert`](`GridMap::extract_convert`).
    pub fn extract_convert_column_major(&self) -> impl Iterator<Item = (f32, f32, VisValue)> + '_ {
        self.extract_raw_column_major()
            .map(|x| (x.0, x.1, (self.converter)(x.2)))
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

    #[allow(clippy::ptr_arg)]
    fn sum(vec: &Vec<f32>) -> f32 {
        vec.iter().sum()
    }

    #[test]
    fn new() {
        let gridmap = GridMap::new((-2.0, 7.0), (3.0, 6.0), (0.15, 0.20), sum).unwrap();

        assert_eq!(gridmap.span_x, (-2.0, 7.0));
        assert_eq!(gridmap.span_y, (3.0, 6.0));
        assert_eq!(gridmap.tile_dim, (0.15, 0.20));
        assert_eq!(gridmap.n_tiles(), 61 * 16);
        assert_eq!(gridmap.n_tiles_x(), 61);
        assert_eq!(gridmap.n_tiles_y(), 16);
        assert_eq!(gridmap.values[(5, 8)], Vec::<f32>::new());
    }

    #[test]
    fn new_failures() {
        assert!(matches!(
            GridMap::<Vec<f32>, f32, fn(&Vec<f32>) -> f32>::new(
                (2.0, -2.0),
                (3.0, 6.0),
                (0.15, 0.20),
                sum
            ),
            Err(GridMapError::InvalidSpan)
        ));
        assert!(matches!(
            GridMap::<Vec<f32>, f32, fn(&Vec<f32>) -> f32>::new(
                (-2.0, 7.0),
                (3.0, -6.0),
                (0.15, 0.20),
                sum
            ),
            Err(GridMapError::InvalidSpan)
        ));
        assert!(matches!(
            GridMap::<Vec<f32>, f32, fn(&Vec<f32>) -> f32>::new(
                (-2.0, 7.0),
                (3.0, 6.0),
                (0.0, 0.20),
                sum
            ),
            Err(GridMapError::InvalidGridTile)
        ));
        assert!(matches!(
            GridMap::<Vec<f32>, f32, fn(&Vec<f32>) -> f32>::new(
                (-2.0, 7.0),
                (3.0, 6.0),
                (0.15, 0.0),
                sum
            ),
            Err(GridMapError::InvalidGridTile)
        ));
        assert!(matches!(
            GridMap::<Vec<f32>, f32, fn(&Vec<f32>) -> f32>::new(
                (-2.0, 7.0),
                (3.0, 6.0),
                (10.5, 0.20),
                sum
            ),
            Err(GridMapError::InvalidGridTile)
        ));
        assert!(matches!(
            GridMap::<Vec<f32>, f32, fn(&Vec<f32>) -> f32>::new(
                (-2.0, 7.0),
                (3.0, 6.0),
                (0.15, 3.10),
                sum
            ),
            Err(GridMapError::InvalidGridTile)
        ));
    }

    #[test]
    fn from_box() {
        let simbox = SimBox::from([10.0, 20.0, 14.0]);
        let gridmap = GridMap::from_box(&simbox, (0.15, 0.20), sum).unwrap();

        assert_eq!(gridmap.span_x, (0.0, 10.0));
        assert_eq!(gridmap.span_y, (0.0, 20.0));
        assert_eq!(gridmap.tile_dim, (0.15, 0.20));
        assert_eq!(gridmap.n_tiles(), 68 * 101);
        assert_eq!(gridmap.n_tiles_x(), 68);
        assert_eq!(gridmap.n_tiles_y(), 101);
        assert_eq!(gridmap.values[(5, 8)], Vec::<f32>::new());
    }

    #[test]
    fn from_vec_column_major() {
        let values = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let gridmap = GridMap::from_vec(
            (1.0, 2.0),
            (1.0, 2.5),
            (1.0, 0.5),
            values,
            DataOrder::ColumnMajor,
            usize::to_owned,
        )
        .unwrap();

        assert_eq!(gridmap.n_tiles_x(), 2);
        assert_eq!(gridmap.n_tiles_y(), 4);
        assert_eq!(gridmap.n_tiles(), 8);

        assert_eq!(*gridmap.get_at(1.0, 1.0).unwrap(), 1);
        assert_eq!(*gridmap.get_at(2.0, 1.0).unwrap(), 2);
        assert_eq!(*gridmap.get_at(1.0, 1.5).unwrap(), 3);
        assert_eq!(*gridmap.get_at(2.0, 1.5).unwrap(), 4);
        assert_eq!(*gridmap.get_at(1.0, 2.0).unwrap(), 5);
        assert_eq!(*gridmap.get_at(2.0, 2.0).unwrap(), 6);
        assert_eq!(*gridmap.get_at(1.0, 2.5).unwrap(), 7);
        assert_eq!(*gridmap.get_at(2.0, 2.5).unwrap(), 8);
    }

    #[test]
    fn from_vec_row_major() {
        let values = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let gridmap = GridMap::from_vec(
            (1.0, 2.0),
            (1.0, 2.5),
            (1.0, 0.5),
            values,
            DataOrder::RowMajor,
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

    #[allow(clippy::ptr_arg)]
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
    fn from_file_column_major() {
        let gridmap = GridMap::from_file(
            "test_files/gridmaps/map_column_major.dat",
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
        assert_eq!(*gridmap.get_at(2.0, 1.0).unwrap(), Vec::<usize>::new());
    }

    #[test]
    fn from_file_column_major_not_increasing() {
        match GridMap::from_file(
            "test_files/gridmaps/map_column_major_decreasing.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        ) {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(GridMapError::NotIncreasing(x, y)) => {
                assert_eq!(x, String::from("1"));
                assert_eq!(y, String::from("2"));
            }
            Err(e) => panic!("Invalid error type returned: '{}'", e),
        }
    }

    #[test]
    fn from_file_column_major_not_fully_increasing() {
        match GridMap::from_file(
            "test_files/gridmaps/map_column_major_not_fully_increasing.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        ) {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(GridMapError::NotIncreasing(x, y)) => {
                assert_eq!(x, String::from("0"));
                assert_eq!(y, String::from("1"));
            }
            Err(e) => panic!("Invalid error type returned: '{}'", e),
        }
    }

    #[test]
    fn from_file_column_major_x_inconsistency() {
        match GridMap::from_file(
            "test_files/gridmaps/map_column_major_x_inconsistency.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        ) {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(GridMapError::InvalidCoordinates(x)) => {
                assert_eq!(x, String::from("2"));
            }
            Err(e) => panic!("Invalid error type returned: '{}'", e),
        }
    }

    #[test]
    fn from_file_column_major_y_inconsistency() {
        match GridMap::from_file(
            "test_files/gridmaps/map_column_major_y_inconsistency.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        ) {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(GridMapError::InvalidCoordinates(x)) => {
                assert_eq!(x, String::from("3.5"));
            }
            Err(e) => panic!("Invalid error type returned: '{}'", e),
        }
    }

    #[test]
    fn from_file_column_major_redefinition() {
        match GridMap::from_file(
            "test_files/gridmaps/map_column_major_redefinition.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        ) {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(GridMapError::PointDefinedMultipleTimes(x, y)) => {
                assert_eq!(x, String::from("0"));
                assert_eq!(y, String::from("0"));
            }
            Err(e) => panic!("Invalid error type returned: '{}'", e),
        }
    }

    #[test]
    fn from_file_row_major() {
        let gridmap = GridMap::from_file(
            "test_files/gridmaps/map_row_major.dat",
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
        assert_eq!(*gridmap.get_at(2.0, 1.0).unwrap(), Vec::<usize>::new());
    }

    #[test]
    fn from_file_row_major_not_increasing() {
        match GridMap::from_file(
            "test_files/gridmaps/map_row_major_decreasing.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        ) {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(GridMapError::NotIncreasing(x, y)) => {
                assert_eq!(x, String::from("1"));
                assert_eq!(y, String::from("2"));
            }
            Err(e) => panic!("Invalid error type returned: '{}'", e),
        }
    }

    #[test]
    fn from_file_row_major_x_inconsistency() {
        match GridMap::from_file(
            "test_files/gridmaps/map_row_major_x_inconsistency.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        ) {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(GridMapError::InvalidCoordinates(x)) => {
                assert_eq!(x, String::from("1"));
            }
            Err(e) => panic!("Invalid error type returned: '{}'", e),
        }
    }

    #[test]
    fn from_file_row_major_y_inconsistency() {
        match GridMap::from_file(
            "test_files/gridmaps/map_row_major_y_inconsistency.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        ) {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(GridMapError::InvalidCoordinates(x)) => {
                assert_eq!(x, String::from("0"));
            }
            Err(e) => panic!("Invalid error type returned: '{}'", e),
        }
    }

    #[test]
    fn from_file_row_major_redefinition() {
        match GridMap::from_file(
            "test_files/gridmaps/map_row_major_redefinition.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        ) {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(GridMapError::PointDefinedMultipleTimes(x, y)) => {
                assert_eq!(x, String::from("0"));
                assert_eq!(y, String::from("1"));
            }
            Err(e) => panic!("Invalid error type returned: '{}'", e),
        }
    }

    #[test]
    fn all_coordinates_same() {
        match GridMap::from_file(
            "test_files/gridmaps/all_coordinates_same.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        ) {
            Ok(_) => panic!("Should have failed, but succeeded."),
            Err(GridMapError::InvalidCoordinates(x)) => {
                assert_eq!(x, String::from("0"));
            }
            Err(e) => panic!("Invalid error type returned: '{}'", e),
        }
    }

    #[test]
    fn clear() {
        let values = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let mut gridmap = GridMap::from_vec(
            (1.0, 2.0),
            (1.0, 2.5),
            (1.0, 0.5),
            values,
            DataOrder::ColumnMajor,
            usize::to_owned,
        )
        .unwrap();

        gridmap.clear();
        assert!(gridmap.values.iter().all(|&x| x == 0));

        let mut gridmap = GridMap::from_file(
            "test_files/gridmaps/map_column_major.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        )
        .unwrap();

        gridmap.clear();
        assert!(gridmap.values.iter().all(|x| x.is_empty()));
    }

    #[test]
    fn coord2index() {
        let gridmap = GridMap::new((-2.0, 7.0), (3.0, 6.0), (0.15, 0.20), sum).unwrap();

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
        let mut gridmap = GridMap::new((-2.0, 7.0), (3.0, 6.0), (0.15, 0.20), sum).unwrap();

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
        let gridmap = GridMap::from_file(
            "test_files/gridmaps/map_row_major.dat",
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
        assert_eq!(vec[1], (0.0, 1.0, &vec![4, 8, 12]));
        assert_eq!(vec[2], (1.0, 0.0, &vec![5, 4]));
        assert_eq!(vec[3], (1.0, 1.0, &vec![2, 5]));
        assert_eq!(vec[4], (2.0, 0.0, &vec![43, 23, 21]));
        assert_eq!(vec[5], (2.0, 1.0, &vec![]));
    }

    #[test]
    fn extract_raw_column_major() {
        let gridmap = GridMap::from_file(
            "test_files/gridmaps/map_column_major.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        )
        .unwrap();

        let vec = gridmap
            .extract_raw_column_major()
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
        let gridmap = GridMap::from_file(
            "test_files/gridmaps/map_column_major.dat",
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
        assert_eq!(vec[1], (0.0, 1.0, 24));
        assert_eq!(vec[2], (1.0, 0.0, 9));
        assert_eq!(vec[3], (1.0, 1.0, 7));
        assert_eq!(vec[4], (2.0, 0.0, 87));
        assert_eq!(vec[5], (2.0, 1.0, 0));
    }

    #[test]
    fn extract_convert_column_major() {
        let gridmap = GridMap::from_file(
            "test_files/gridmaps/map_row_major.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        )
        .unwrap();

        let vec = gridmap
            .extract_convert_column_major()
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
        let gridmap = GridMap::from_file(
            "test_files/gridmaps/map_row_major.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        )
        .unwrap();

        let mut output = BufWriter::new(Vec::new());
        gridmap.write_map_raw(&mut output).unwrap();
        let output_string = String::from_utf8(output.into_inner().unwrap()).unwrap();

        let expected = "  0.000000   0.000000 [10]\n  0.000000   1.000000 [4, 8, 12]\n  1.000000   0.000000 [5, 4]\n  1.000000   1.000000 [2, 5]\n  2.000000   0.000000 [43, 23, 21]\n  2.000000   1.000000 []\n";

        assert_eq!(output_string, expected);
    }

    #[test]
    fn write_raw_column_major() {
        let gridmap = GridMap::from_file(
            "test_files/gridmaps/map_column_major.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        )
        .unwrap();

        let mut output = BufWriter::new(Vec::new());
        gridmap.write_map_raw_column_major(&mut output).unwrap();
        let output_string = String::from_utf8(output.into_inner().unwrap()).unwrap();

        let expected = "  0.000000   0.000000 [10]\n  1.000000   0.000000 [5, 4]\n  2.000000   0.000000 [43, 23, 21]\n  0.000000   1.000000 [4, 8, 12]\n  1.000000   1.000000 [2, 5]\n  2.000000   1.000000 []\n";

        assert_eq!(output_string, expected);
    }

    #[test]
    fn write() {
        let gridmap = GridMap::from_file(
            "test_files/gridmaps/map_column_major.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        )
        .unwrap();

        let mut output = BufWriter::new(Vec::new());
        gridmap.write_map(&mut output).unwrap();
        let output_string = String::from_utf8(output.into_inner().unwrap()).unwrap();

        let expected = "  0.000000   0.000000 10\n  0.000000   1.000000 24\n  1.000000   0.000000 9\n  1.000000   1.000000 7\n  2.000000   0.000000 87\n  2.000000   1.000000 0\n";

        assert_eq!(output_string, expected);
    }

    #[test]
    fn write_column_major() {
        let gridmap = GridMap::from_file(
            "test_files/gridmaps/map_row_major.dat",
            sum_usize,
            &['|'],
            parse_vec,
            &["@", "#"],
        )
        .unwrap();

        let mut output = BufWriter::new(Vec::new());
        gridmap.write_map_column_major(&mut output).unwrap();
        let output_string = String::from_utf8(output.into_inner().unwrap()).unwrap();

        let expected = "  0.000000   0.000000 10\n  1.000000   0.000000 9\n  2.000000   0.000000 87\n  0.000000   1.000000 24\n  1.000000   1.000000 7\n  2.000000   1.000000 0\n";

        assert_eq!(output_string, expected);
    }

    #[test]
    fn is_inside() {
        let gridmap = GridMap::new((-2.0, 7.0), (3.0, 6.0), (0.15, 0.20), sum).unwrap();

        assert!(gridmap.is_inside(-1.8, 4.5));
        assert!(gridmap.is_inside(-2.05, 3.09));
        assert!(gridmap.is_inside(7.05, 6.09));
        assert!(!gridmap.is_inside(7.11, 4.5));
        assert!(!gridmap.is_inside(0.11, 2.8));
    }

    #[test]
    fn get_tile() {
        let gridmap = GridMap::new((-2.0, 7.0), (3.0, 6.0), (0.15, 0.20), sum).unwrap();

        assert_eq!(gridmap.get_tile(-1.8, 4.5).unwrap(), (-1.85, 4.6));
        assert_eq!(gridmap.get_tile(5.8, 4.2).unwrap(), (5.8, 4.2));
        assert_eq!(gridmap.get_tile(5.7843, 4.12374).unwrap(), (5.8, 4.2));
        assert!(gridmap.get_tile(7.11, 4.5).is_none());
        assert!(gridmap.get_tile(0.11, 2.8).is_none());
    }

    #[test]
    fn using_closure_as_converter() {
        let number = 17;
        let sum_and_add_number = move |raw: &Vec<usize>| raw.iter().sum::<usize>() + number;

        let gridmap = GridMap::from_file(
            "test_files/gridmaps/map_row_major.dat",
            sum_and_add_number,
            &['|'],
            parse_vec,
            &["@", "#"],
        )
        .unwrap();

        let vec = gridmap
            .extract_convert()
            .collect::<Vec<(f32, f32, usize)>>();
        assert_eq!(vec[0], (0.0, 0.0, 27));
        assert_eq!(vec[1], (0.0, 1.0, 41));
        assert_eq!(vec[2], (1.0, 0.0, 26));
        assert_eq!(vec[3], (1.0, 1.0, 24));
        assert_eq!(vec[4], (2.0, 0.0, 104));
        assert_eq!(vec[5], (2.0, 1.0, 17));
    }
}
