// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Implementation of some higher level functions for `groan_rs` programs.

use std::ops::Deref;

use crate::errors::{AtomError, GroupError, ParseNdxError};
use crate::structures::{dimension::Dimension, vector3d::Vector3D};
use crate::system::System;

/// ## Methods providing higher level utilities to the `System` structure.
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
    /// - `Ok` if centering was successful.
    /// - `GroupError::NotFound` if the group does not exist.
    /// - `GroupError::EmptyGroup` if the group contains no atoms.
    /// - `GroupError::InvalidSimBox` if the system has no simulation box or the simulation box is not orthogonal.
    /// - `GroupError::InvalidPosition` if any of the atoms in the system has an undefined position.
    ///
    /// ## Example
    /// Center the group `Protein` in the xy-plane.
    /// The positions of the atoms along the z-dimension will be unchanged.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// // load system from a gro file
    /// let mut system = System::from_file("system.gro").unwrap();
    ///  
    /// // create a group Protein by autodetecting protein atoms
    /// system.group_create("Protein", "@protein").unwrap();
    ///
    /// // center the group 'Protein'
    /// system.atoms_center("Protein", Dimension::XY).unwrap();
    /// ```
    ///
    /// ## Notes
    /// - The Bai-Breen method may introduce some (typically minor) error into the calculated center of geometry.
    ///   This method should be used for visual centering, not if you require the group to be _precisely_ in the center of the box!
    pub fn atoms_center(
        &mut self,
        reference: &str,
        dimension: Dimension,
    ) -> Result<(), GroupError> {
        let reference_center = self.group_estimate_center(reference)?;

        let box_center = self.get_box_center().map_err(GroupError::InvalidSimBox)?;
        let mut shift = Vector3D(box_center.deref() - reference_center.deref());

        shift.filter(dimension);
        match self.atoms_translate(&shift) {
            Ok(_) => Ok(()),
            Err(AtomError::InvalidSimBox(e)) => Err(GroupError::InvalidSimBox(e)),
            Err(AtomError::InvalidPosition(e)) => Err(GroupError::InvalidPosition(e)),
            _ => panic!("FATAL GROAN ERROR | System::atoms_center | Invalid error type returned from `System::atoms_translate.`"),
        }
    }

    /// Place the center of mass of the reference group to the center of a simulation box.
    /// This shifts all other atoms accordingly and wraps them into the simulation box.
    /// In other words, this function performs an operation commonly called 'centering'.
    /// As this function uses Bai-Breen method for calculating the center of geometry in periodic systems,
    /// the reference group can be of any size and atom distribution.
    ///
    /// Unlike `System::atoms_center`, this function uses center of *mass* of the reference group,
    /// instead of center of geometry.
    ///
    /// ## Arguments
    /// - `reference` - Name of the group to center.
    /// - `dimension` - Dimension enum for specifying centering dimensions.
    ///
    /// ## Returns
    /// - `Ok` if centering was successful.
    /// - `GroupError::NotFound` if the group does not exist.
    /// - `GroupError::EmptyGroup` if the group contains no atoms.
    /// - `GroupError::InvalidSimBox` if the system has no simulation box or the simulation box is not orthogonal.
    /// - `GroupError::InvalidPosition` if any of the atoms in the system has an undefined position.
    /// - `GroupError::InvalidMass` if any of the atoms in the reference group has no mass.
    ///
    /// ## Example
    /// Center the group `Protein` in the xy-plane. Use center of mass as the reference point.
    /// The positions of the atoms along the z-dimension will be unchanged.
    /// ```no_run
    /// # use groan_rs::prelude::*;
    /// #
    /// // load system from a gro file
    /// let mut system = System::from_file("system.gro").unwrap();
    ///  
    /// // create a group Protein by autodetecting protein atoms
    /// system.group_create("Protein", "@protein").unwrap();
    ///
    /// // center the group 'Protein'
    /// system.atoms_center_mass("Protein", Dimension::XY).unwrap();
    /// ```
    ///
    /// ## Notes
    /// - The Bai-Breen method may introduce some (typically minor) error into the calculated center of geometry.
    ///   This method should be used for visual centering, not if you require the group to be _precisely_ in the center of the box!
    pub fn atoms_center_mass(
        &mut self,
        reference: &str,
        dimension: Dimension,
    ) -> Result<(), GroupError> {
        let reference_com = self.group_get_com(reference)?;

        let box_center = self.get_box_center().map_err(GroupError::InvalidSimBox)?;
        let mut shift = Vector3D(box_center.deref() - reference_com.deref());

        shift.filter(dimension);
        match self.atoms_translate(&shift) {
            Ok(_) => Ok(()),
            Err(AtomError::InvalidSimBox(e)) => Err(GroupError::InvalidSimBox(e)),
            Err(AtomError::InvalidPosition(e)) => Err(GroupError::InvalidPosition(e)),
            _ => panic!("FATAL GROAN ERROR | System::atoms_center_mass | Invalid error type returned from `System::atoms_translate.`"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        errors::{MassError, PositionError, SimBoxError},
        structures::element::Elements,
    };
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

        assert_approx_eq!(f32, atom1.unwrap().x, 9.497);
        assert_approx_eq!(f32, atom1.unwrap().y, 1.989);
        assert_approx_eq!(f32, atom1.unwrap().z, 7.498);

        assert_approx_eq!(f32, atom2.unwrap().x, 8.829);
        assert_approx_eq!(f32, atom2.unwrap().y, 11.186);
        assert_approx_eq!(f32, atom2.unwrap().z, 2.075);
    }

    #[test]
    fn atoms_center_x() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::X).unwrap();

        assert_approx_eq!(
            f32,
            system.group_get_center("Protein").unwrap().x,
            system.get_box_center().unwrap().x
        );

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.unwrap().x, 6.1465545);
        assert_approx_eq!(f32, atom1.unwrap().y, 1.989);
        assert_approx_eq!(f32, atom1.unwrap().z, 7.498);

        assert_approx_eq!(f32, atom2.unwrap().x, 5.478555);
        assert_approx_eq!(f32, atom2.unwrap().y, 11.186);
        assert_approx_eq!(f32, atom2.unwrap().z, 2.075);
    }

    #[test]
    fn atoms_center_y() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::Y).unwrap();

        assert_approx_eq!(
            f32,
            system.group_get_center("Protein").unwrap().y,
            system.get_box_center().unwrap().y
        );

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.unwrap().x, 9.497);
        assert_approx_eq!(f32, atom1.unwrap().y, 6.033055);
        assert_approx_eq!(f32, atom1.unwrap().z, 7.498);

        assert_approx_eq!(f32, atom2.unwrap().x, 8.829);
        assert_approx_eq!(f32, atom2.unwrap().y, 2.2167444);
        assert_approx_eq!(f32, atom2.unwrap().z, 2.075);
    }

    #[test]
    fn atoms_center_z() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::Z).unwrap();

        assert_approx_eq!(
            f32,
            system.group_get_center("Protein").unwrap().z,
            system.get_box_center().unwrap().z
        );

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.unwrap().x, 9.497);
        assert_approx_eq!(f32, atom1.unwrap().y, 1.989);
        assert_approx_eq!(f32, atom1.unwrap().z, 7.6634398);

        assert_approx_eq!(f32, atom2.unwrap().x, 8.829);
        assert_approx_eq!(f32, atom2.unwrap().y, 11.186);
        assert_approx_eq!(f32, atom2.unwrap().z, 2.2404397);
    }

    #[test]
    fn atoms_center_xy() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::XY).unwrap();

        let group_center = system.group_get_center("Protein").unwrap();
        assert_approx_eq!(f32, group_center.x, system.get_box_center().unwrap().x);
        assert_approx_eq!(f32, group_center.y, system.get_box_center().unwrap().y);

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.unwrap().x, 6.1465545);
        assert_approx_eq!(f32, atom1.unwrap().y, 6.033055);
        assert_approx_eq!(f32, atom1.unwrap().z, 7.498);

        assert_approx_eq!(f32, atom2.unwrap().x, 5.478555);
        assert_approx_eq!(f32, atom2.unwrap().y, 2.2167444);
        assert_approx_eq!(f32, atom2.unwrap().z, 2.075);
    }

    #[test]
    fn atoms_center_xz() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::XZ).unwrap();

        let group_center = system.group_get_center("Protein").unwrap();
        assert_approx_eq!(f32, group_center.x, system.get_box_center().unwrap().x);
        assert_approx_eq!(f32, group_center.z, system.get_box_center().unwrap().z);

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.unwrap().x, 6.1465545);
        assert_approx_eq!(f32, atom1.unwrap().y, 1.989);
        assert_approx_eq!(f32, atom1.unwrap().z, 7.6634398);

        assert_approx_eq!(f32, atom2.unwrap().x, 5.478555);
        assert_approx_eq!(f32, atom2.unwrap().y, 11.186);
        assert_approx_eq!(f32, atom2.unwrap().z, 2.2404397);
    }

    #[test]
    fn atoms_center_yz() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::YZ).unwrap();

        let group_center = system.group_get_center("Protein").unwrap();
        assert_approx_eq!(f32, group_center.y, system.get_box_center().unwrap().y);
        assert_approx_eq!(f32, group_center.z, system.get_box_center().unwrap().z);

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.unwrap().x, 9.497);
        assert_approx_eq!(f32, atom1.unwrap().y, 6.033055);
        assert_approx_eq!(f32, atom1.unwrap().z, 7.6634398);

        assert_approx_eq!(f32, atom2.unwrap().x, 8.829);
        assert_approx_eq!(f32, atom2.unwrap().y, 2.2167444);
        assert_approx_eq!(f32, atom2.unwrap().z, 2.2404397);
    }

    #[test]
    fn atoms_center_xyz() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        system.atoms_center("Protein", Dimension::XYZ).unwrap();

        let group_center = system.group_get_center("Protein").unwrap();
        assert_approx_eq!(f32, group_center.x, system.get_box_center().unwrap().x);
        assert_approx_eq!(f32, group_center.y, system.get_box_center().unwrap().y);
        assert_approx_eq!(f32, group_center.z, system.get_box_center().unwrap().z);

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.unwrap().x, 6.1465545);
        assert_approx_eq!(f32, atom1.unwrap().y, 6.033055);
        assert_approx_eq!(f32, atom1.unwrap().z, 7.6634398);

        assert_approx_eq!(f32, atom2.unwrap().x, 5.478555);
        assert_approx_eq!(f32, atom2.unwrap().y, 2.2167444);
        assert_approx_eq!(f32, atom2.unwrap().z, 2.2404397);
    }

    #[test]
    fn atoms_center_fail_group() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        match system.atoms_center("Nonexistent", Dimension::XYZ) {
            Ok(_) => panic!("Function should have failed."),
            Err(GroupError::NotFound(x)) => assert_eq!(x, "Nonexistent"),
            Err(e) => panic!(
                "Function failed successfully, but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn atoms_center_fail_empty() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        match system.atoms_center("Backbone", Dimension::XYZ) {
            Ok(_) => panic!("Function should have failed."),
            Err(GroupError::EmptyGroup(x)) => assert_eq!(x, "Backbone"),
            Err(e) => panic!(
                "Function failed successfully, but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn atoms_center_fail_simbox() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();
        system.reset_box();

        match system.atoms_center("Protein", Dimension::XYZ) {
            Ok(_) => panic!("Function should have failed."),
            Err(GroupError::InvalidSimBox(SimBoxError::DoesNotExist)) => (),
            Err(e) => panic!(
                "Function failed successfully, but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn atoms_center_fail_position() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();
        system.get_atom_mut(15).unwrap().reset_position();

        match system.atoms_center("Protein", Dimension::XYZ) {
            Ok(_) => panic!("Function should have failed."),
            Err(GroupError::InvalidPosition(PositionError::NoPosition(x))) => assert_eq!(x, 15),
            Err(e) => panic!(
                "Function failed successfully, but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn atoms_center_mass_x() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.guess_elements(Elements::default()).unwrap();
        system.group_create("Protein", "@protein").unwrap();

        system.atoms_center_mass("Protein", Dimension::X).unwrap();

        assert_approx_eq!(
            f32,
            system.group_get_com("Protein").unwrap().x,
            system.get_box_center().unwrap().x
        );

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.unwrap().x, 3.456437);
        assert_approx_eq!(f32, atom1.unwrap().y, 3.899);
        assert_approx_eq!(f32, atom1.unwrap().z, 4.993);

        assert_approx_eq!(f32, atom2.unwrap().x, 2.0444372);
        assert_approx_eq!(f32, atom2.unwrap().y, 3.823);
        assert_approx_eq!(f32, atom2.unwrap().z, 0.378);
    }

    #[test]
    fn atoms_center_mass_y() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.guess_elements(Elements::default()).unwrap();
        system.group_create("Protein", "@protein").unwrap();

        system.atoms_center_mass("Protein", Dimension::Y).unwrap();

        assert_approx_eq!(
            f32,
            system.group_get_com("Protein").unwrap().y,
            system.get_box_center().unwrap().y
        );

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.unwrap().x, 4.322);
        assert_approx_eq!(f32, atom1.unwrap().y, 3.475028);
        assert_approx_eq!(f32, atom1.unwrap().z, 4.993);

        assert_approx_eq!(f32, atom2.unwrap().x, 2.910);
        assert_approx_eq!(f32, atom2.unwrap().y, 3.399028);
        assert_approx_eq!(f32, atom2.unwrap().z, 0.378);
    }

    #[test]
    fn atoms_center_mass_z() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.guess_elements(Elements::default()).unwrap();
        system.group_create("Protein", "@protein").unwrap();

        system.atoms_center_mass("Protein", Dimension::Z).unwrap();

        assert_approx_eq!(
            f32,
            system.group_get_com("Protein").unwrap().z,
            system.get_box_center().unwrap().z
        );

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.unwrap().x, 4.322);
        assert_approx_eq!(f32, atom1.unwrap().y, 3.899);
        assert_approx_eq!(f32, atom1.unwrap().z, 5.4376106);

        assert_approx_eq!(f32, atom2.unwrap().x, 2.910);
        assert_approx_eq!(f32, atom2.unwrap().y, 3.823);
        assert_approx_eq!(f32, atom2.unwrap().z, 0.82261086);
    }

    #[test]
    fn atoms_center_mass_xy() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.guess_elements(Elements::default()).unwrap();
        system.group_create("Protein", "@protein").unwrap();

        system.atoms_center_mass("Protein", Dimension::XY).unwrap();

        assert_approx_eq!(
            f32,
            system.group_get_com("Protein").unwrap().x,
            system.get_box_center().unwrap().x
        );

        assert_approx_eq!(
            f32,
            system.group_get_com("Protein").unwrap().y,
            system.get_box_center().unwrap().y
        );

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.unwrap().x, 3.456437);
        assert_approx_eq!(f32, atom1.unwrap().y, 3.475028);
        assert_approx_eq!(f32, atom1.unwrap().z, 4.993);

        assert_approx_eq!(f32, atom2.unwrap().x, 2.0444372);
        assert_approx_eq!(f32, atom2.unwrap().y, 3.399028);
        assert_approx_eq!(f32, atom2.unwrap().z, 0.378);
    }

    #[test]
    fn atoms_center_mass_xz() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.guess_elements(Elements::default()).unwrap();
        system.group_create("Protein", "@protein").unwrap();

        system.atoms_center_mass("Protein", Dimension::XZ).unwrap();

        assert_approx_eq!(
            f32,
            system.group_get_com("Protein").unwrap().x,
            system.get_box_center().unwrap().x
        );

        assert_approx_eq!(
            f32,
            system.group_get_com("Protein").unwrap().z,
            system.get_box_center().unwrap().z
        );

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.unwrap().x, 3.456437);
        assert_approx_eq!(f32, atom1.unwrap().y, 3.899);
        assert_approx_eq!(f32, atom1.unwrap().z, 5.4376106);

        assert_approx_eq!(f32, atom2.unwrap().x, 2.0444372);
        assert_approx_eq!(f32, atom2.unwrap().y, 3.823);
        assert_approx_eq!(f32, atom2.unwrap().z, 0.82261086);
    }

    #[test]
    fn atoms_center_mass_yz() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.guess_elements(Elements::default()).unwrap();
        system.group_create("Protein", "@protein").unwrap();

        system.atoms_center_mass("Protein", Dimension::YZ).unwrap();

        assert_approx_eq!(
            f32,
            system.group_get_com("Protein").unwrap().y,
            system.get_box_center().unwrap().y
        );

        assert_approx_eq!(
            f32,
            system.group_get_com("Protein").unwrap().z,
            system.get_box_center().unwrap().z
        );

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.unwrap().x, 4.322);
        assert_approx_eq!(f32, atom1.unwrap().y, 3.475028);
        assert_approx_eq!(f32, atom1.unwrap().z, 5.4376106);

        assert_approx_eq!(f32, atom2.unwrap().x, 2.910);
        assert_approx_eq!(f32, atom2.unwrap().y, 3.399028);
        assert_approx_eq!(f32, atom2.unwrap().z, 0.82261086);
    }

    #[test]
    fn atoms_center_mass_xyz() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.guess_elements(Elements::default()).unwrap();
        system.group_create("Protein", "@protein").unwrap();

        system.atoms_center_mass("Protein", Dimension::XYZ).unwrap();

        assert_approx_eq!(
            f32,
            system.group_get_com("Protein").unwrap().x,
            system.get_box_center().unwrap().x
        );

        assert_approx_eq!(
            f32,
            system.group_get_com("Protein").unwrap().y,
            system.get_box_center().unwrap().y
        );

        assert_approx_eq!(
            f32,
            system.group_get_com("Protein").unwrap().z,
            system.get_box_center().unwrap().z
        );

        let atom1 = system.atoms_iter().next().unwrap().get_position();
        let atom2 = system.atoms_iter().last().unwrap().get_position();

        assert_approx_eq!(f32, atom1.unwrap().x, 3.456437);
        assert_approx_eq!(f32, atom1.unwrap().y, 3.475028);
        assert_approx_eq!(f32, atom1.unwrap().z, 5.4376106);

        assert_approx_eq!(f32, atom2.unwrap().x, 2.0444372);
        assert_approx_eq!(f32, atom2.unwrap().y, 3.399028);
        assert_approx_eq!(f32, atom2.unwrap().z, 0.82261086);
    }

    #[test]
    fn atoms_center_mass_fail_group() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.guess_elements(Elements::default()).unwrap();

        match system.atoms_center_mass("Protein", Dimension::XYZ) {
            Ok(_) => panic!("Function should have failed."),
            Err(GroupError::NotFound(x)) => assert_eq!(x, "Protein"),
            Err(e) => panic!(
                "Function failed successfully, but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn atoms_center_mass_fail_empty() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.guess_elements(Elements::default()).unwrap();
        system.group_create("Empty", "resname X").unwrap();

        match system.atoms_center_mass("Empty", Dimension::XYZ) {
            Ok(_) => panic!("Function should have failed."),
            Err(GroupError::EmptyGroup(x)) => assert_eq!(x, "Empty"),
            Err(e) => panic!(
                "Function failed successfully, but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn atoms_center_mass_fail_simbox() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.guess_elements(Elements::default()).unwrap();
        system.group_create("Protein", "@protein").unwrap();
        system.reset_box();

        match system.atoms_center_mass("Protein", Dimension::XYZ) {
            Ok(_) => panic!("Function should have failed."),
            Err(GroupError::InvalidSimBox(SimBoxError::DoesNotExist)) => (),
            Err(e) => panic!(
                "Function failed successfully, but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn atoms_center_mass_fail_position() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.guess_elements(Elements::default()).unwrap();
        system.group_create("Protein", "@protein").unwrap();
        system.get_atom_mut(15).unwrap().reset_position();

        match system.atoms_center_mass("Protein", Dimension::XYZ) {
            Ok(_) => panic!("Function should have failed."),
            Err(GroupError::InvalidPosition(PositionError::NoPosition(x))) => assert_eq!(x, 15),
            Err(e) => panic!(
                "Function failed successfully, but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }

    #[test]
    fn atoms_center_mass_fail_mass() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.gro").unwrap();
        system.group_create("Protein", "@protein").unwrap();

        match system.atoms_center_mass("Protein", Dimension::XYZ) {
            Ok(_) => panic!("Function should have failed."),
            Err(GroupError::InvalidMass(MassError::NoMass(x))) => assert_eq!(x, 0),
            Err(e) => panic!(
                "Function failed successfully, but incorrect error type `{:?}` was returned.",
                e
            ),
        }
    }
}
