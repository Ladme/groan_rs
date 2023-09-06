// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of some higher level functions for `groan_rs` programs.

use crate::dimension::Dimension;
use crate::errors::{GroupError, ParseNdxError};
use crate::system::System;
use crate::vector3d::Vector3D;

impl System {
    /// Read ndx file using the provided default if index is None.
    ///
    /// ## Returns
    /// `Ok` if parsing was successful or if `index` is `None` and
    /// a) the ndx file does not exist, or
    /// b) the ndx file does not match the gro file.
    /// Also returns `Ok` if a warning has been detected, but prints it to `stderr`.
    /// Else returns `ParseNdxError`.
    ///
    /// ## Notes
    /// - Prints any warning it detects, but still returns `Ok`.
    ///
    pub fn read_ndx_with_default(
        &mut self,
        index: &Option<String>,
        default: &str,
    ) -> Result<(), ParseNdxError> {
        let mut def = true;
        let filename = match index {
            Some(x) => {
                def = false;
                x.to_string()
            }
            None => default.to_string(),
        };

        match self.read_ndx(filename) {
            Ok(_) => Ok(()),
            Err(ParseNdxError::DuplicateGroupsWarning(e)) => {
                eprintln!("{}", ParseNdxError::DuplicateGroupsWarning(e));
                Ok(())
            }

            Err(ParseNdxError::InvalidNamesWarning(e)) => {
                eprintln!("{}", ParseNdxError::InvalidNamesWarning(e));
                Ok(())
            }

            // this error is ignored if using default
            Err(ParseNdxError::InvalidAtomIndex(e)) => {
                if !def {
                    Err(ParseNdxError::InvalidAtomIndex(e))
                } else {
                    Ok(())
                }
            }

            // this error is ignored if using default
            Err(ParseNdxError::FileNotFound(e)) => {
                if !def {
                    Err(ParseNdxError::FileNotFound(e))
                } else {
                    Ok(())
                }
            }
            Err(e) => Err(e),
        }
    }

    /// Place the center of geometry of the reference group to the center of a simulation box.
    /// This shifts all other atoms accordingly and wraps them into the simulation box.
    /// In other words, this function performs an operation commonly called 'centering'.
    /// As this function uses Bai-Breen method for calculating the center of geometry in periodic systems,
    /// the reference group can be of any size and atom distribution.
    ///
    /// ## Arguments
    /// - `reference` - Name of the group to center.
    /// - `dimension` - Dimension enum for specifying centering dimensions.
    ///
    /// ## Returns
    /// `Ok` if centering was successful. `GroupError` if the group does not exist.
    ///
    /// ## Warning
    /// - This function currently only supports orthogonal simulation boxes.
    ///
    /// ## Example
    /// Center the group `Protein` in the xy-plane.
    /// The positions of the atoms along the z-dimension will be unchanged.
    /// ```no_run
    /// use groan_rs::prelude::*;
    ///
    /// // load system from a gro file
    /// let mut system = System::from_file("system.gro").unwrap();
    ///  
    /// // create a group Protein by autodetecting protein atoms
    /// system.group_create("Protein", "@protein").unwrap();
    ///
    /// // center the group 'Protein'
    /// system.atoms_center("Protein", Dimension::XY).unwrap();
    /// ```
    pub fn atoms_center(
        &mut self,
        reference: &str,
        dimension: Dimension,
    ) -> Result<(), GroupError> {
        let reference_center = self.group_get_center(reference)?;

        let box_center = self.get_box_center();
        let mut shift: Vector3D = [
            box_center.x - reference_center.x,
            box_center.y - reference_center.y,
            box_center.z - reference_center.z,
        ]
        .into();

        shift.filter(dimension);
        self.atoms_translate(&shift);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::assert_approx_eq;
    use std::path::Path;

    #[test]
    fn provided() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(
            system.read_ndx_with_default(&Some("test_files/index.ndx".to_string()), "not_used.ndx"),
            Ok(())
        );

        assert!(system.group_exists("System"));
    }

    #[test]
    fn provided_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(
            system.read_ndx_with_default(&Some("nonexistent.ndx".to_string()), "not_used.ndx"),
            Err(ParseNdxError::FileNotFound(Box::from(Path::new(
                "nonexistent.ndx"
            ))))
        );

        assert!(!system.group_exists("System"));
    }

    #[test]
    fn provided_parsing_fails() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(
            system.read_ndx_with_default(
                &Some("test_files/index_invalid_line.ndx".to_string()),
                "not_used.ndx"
            ),
            Err(ParseNdxError::ParseLineErr(
                "  16   17   18   19   20   21   -22   23   24   25   26   27   28   29   30"
                    .to_string()
            ))
        );

        assert!(!system.group_exists("System"));
    }

    #[test]
    fn provided_parsing_warning_duplicate() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(
            system.read_ndx_with_default(
                &Some("test_files/index_duplicate_groups.ndx".to_string()),
                "not_used.ndx"
            ),
            Ok(()),
        );

        assert!(system.group_exists("System"));
    }

    #[test]
    fn provided_parsing_warning_invalid() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(
            system.read_ndx_with_default(
                &Some("test_files/index_invalid_names.ndx".to_string()),
                "not_used.ndx"
            ),
            Ok(())
        );

        assert!(system.group_exists("System"));
    }

    #[test]
    fn default() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(
            system.read_ndx_with_default(&None, "test_files/index.ndx"),
            Ok(())
        );

        assert!(system.group_exists("System"));
    }

    #[test]
    fn default_nonexistent() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(
            system.read_ndx_with_default(&None, "nonexistent.ndx"),
            Ok(())
        );

        assert!(!system.group_exists("System"));
    }

    #[test]
    fn default_parsing_fails() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(
            system.read_ndx_with_default(&None, "test_files/index_invalid_line.ndx"),
            Err(ParseNdxError::ParseLineErr(
                "  16   17   18   19   20   21   -22   23   24   25   26   27   28   29   30"
                    .to_string()
            ))
        );

        assert!(!system.group_exists("System"));
    }

    #[test]
    fn default_parsing_warning_duplicate() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(
            system.read_ndx_with_default(&None, "test_files/index_duplicate_groups.ndx"),
            Ok(())
        );

        assert!(system.group_exists("System"));
    }

    #[test]
    fn default_parsing_warning_invalid() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(
            system.read_ndx_with_default(&None, "test_files/index_invalid_names.ndx"),
            Ok(())
        );

        assert!(system.group_exists("System"));
    }

    #[test]
    fn default_parsing_invalid_atom_index() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(
            system.read_ndx_with_default(&None, "test_files/index_invalid_index1.ndx"),
            Ok(())
        );

        assert!(!system.group_exists("System"));
    }

    #[test]
    fn atoms_center_none() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::None).unwrap();

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.x, 9.497);
        assert_approx_eq!(f32, atom1.y, 1.989);
        assert_approx_eq!(f32, atom1.z, 7.498);

        assert_approx_eq!(f32, atom2.x, 8.829);
        assert_approx_eq!(f32, atom2.y, 11.186);
        assert_approx_eq!(f32, atom2.z, 2.075);
    }

    #[test]
    fn atoms_center_x() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::X).unwrap();

        assert_approx_eq!(
            f32,
            system.group_get_center("Protein").unwrap().x,
            system.get_box_center().x
        );

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.x, 6.1465545);
        assert_approx_eq!(f32, atom1.y, 1.989);
        assert_approx_eq!(f32, atom1.z, 7.498);

        assert_approx_eq!(f32, atom2.x, 5.478555);
        assert_approx_eq!(f32, atom2.y, 11.186);
        assert_approx_eq!(f32, atom2.z, 2.075);
    }

    #[test]
    fn atoms_center_y() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::Y).unwrap();

        assert_approx_eq!(
            f32,
            system.group_get_center("Protein").unwrap().y,
            system.get_box_center().y
        );

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.x, 9.497);
        assert_approx_eq!(f32, atom1.y, 6.033055);
        assert_approx_eq!(f32, atom1.z, 7.498);

        assert_approx_eq!(f32, atom2.x, 8.829);
        assert_approx_eq!(f32, atom2.y, 2.2167444);
        assert_approx_eq!(f32, atom2.z, 2.075);
    }

    #[test]
    fn atoms_center_z() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::Z).unwrap();

        assert_approx_eq!(
            f32,
            system.group_get_center("Protein").unwrap().z,
            system.get_box_center().z
        );

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.x, 9.497);
        assert_approx_eq!(f32, atom1.y, 1.989);
        assert_approx_eq!(f32, atom1.z, 7.6634398);

        assert_approx_eq!(f32, atom2.x, 8.829);
        assert_approx_eq!(f32, atom2.y, 11.186);
        assert_approx_eq!(f32, atom2.z, 2.2404397);
    }

    #[test]
    fn atoms_center_xy() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::XY).unwrap();

        let group_center = system.group_get_center("Protein").unwrap();
        assert_approx_eq!(f32, group_center.x, system.get_box_center().x);
        assert_approx_eq!(f32, group_center.y, system.get_box_center().y);

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.x, 6.1465545);
        assert_approx_eq!(f32, atom1.y, 6.033055);
        assert_approx_eq!(f32, atom1.z, 7.498);

        assert_approx_eq!(f32, atom2.x, 5.478555);
        assert_approx_eq!(f32, atom2.y, 2.2167444);
        assert_approx_eq!(f32, atom2.z, 2.075);
    }

    #[test]
    fn atoms_center_xz() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::XZ).unwrap();

        let group_center = system.group_get_center("Protein").unwrap();
        assert_approx_eq!(f32, group_center.x, system.get_box_center().x);
        assert_approx_eq!(f32, group_center.z, system.get_box_center().z);

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.x, 6.1465545);
        assert_approx_eq!(f32, atom1.y, 1.989);
        assert_approx_eq!(f32, atom1.z, 7.6634398);

        assert_approx_eq!(f32, atom2.x, 5.478555);
        assert_approx_eq!(f32, atom2.y, 11.186);
        assert_approx_eq!(f32, atom2.z, 2.2404397);
    }

    #[test]
    fn atoms_center_yz() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::YZ).unwrap();

        let group_center = system.group_get_center("Protein").unwrap();
        assert_approx_eq!(f32, group_center.y, system.get_box_center().y);
        assert_approx_eq!(f32, group_center.z, system.get_box_center().z);

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.x, 9.497);
        assert_approx_eq!(f32, atom1.y, 6.033055);
        assert_approx_eq!(f32, atom1.z, 7.6634398);

        assert_approx_eq!(f32, atom2.x, 8.829);
        assert_approx_eq!(f32, atom2.y, 2.2167444);
        assert_approx_eq!(f32, atom2.z, 2.2404397);
    }

    #[test]
    fn atoms_center_xyz() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::XYZ).unwrap();

        let group_center = system.group_get_center("Protein").unwrap();
        assert_approx_eq!(f32, group_center.x, system.get_box_center().x);
        assert_approx_eq!(f32, group_center.y, system.get_box_center().y);
        assert_approx_eq!(f32, group_center.z, system.get_box_center().z);

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.x, 6.1465545);
        assert_approx_eq!(f32, atom1.y, 6.033055);
        assert_approx_eq!(f32, atom1.z, 7.6634398);

        assert_approx_eq!(f32, atom2.x, 5.478555);
        assert_approx_eq!(f32, atom2.y, 2.2167444);
        assert_approx_eq!(f32, atom2.z, 2.2404397);
    }
}
