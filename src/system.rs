// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

//! Implementation of the System structure and associated methods.

use std::collections::HashMap;

use crate::atom::Atom;
use crate::errors::{GroupError};
use crate::group::Group;
use crate::simbox::SimBox;
use crate::vector3d::Vector3D;


#[derive(Debug)]
pub struct System {
    name: String,
    atoms: Vec<Atom>,
    simulation_box: SimBox,
    groups: HashMap<String, Group>,
}

impl System {

    /// Create new empty System structure.
    /// ## Notes
    /// - The returned System structure will contain two default groups "all" and "All"
    /// consisting of all the atoms in the system.
    pub fn new(name: &str, atoms: Vec<Atom>, simulation_box: SimBox) -> Self {
        let mut system = System { name: name.to_string(),
                 atoms,
                 simulation_box,
                 groups: HashMap::new()};

        match system.group_create_all() {
            Err(_) => panic!("Internal error. Group `all` or `All` already exists as System is created."),
            Ok(_) => return system,
        }
    }


    /// Get the name of the molecular system.
    pub fn get_name(&self) -> &str { &self.name }

    /// Get immutable reference to the atoms in the system.
    pub fn get_atoms_as_ref(&self) -> &Vec<Atom> { &self.atoms }

    /// Get mutable reference to the atoms in the system.
    pub unsafe fn get_atoms_as_ref_mut(&mut self) -> &mut Vec<Atom> { &mut self.atoms }

    /// Get copy of the atoms in the system.
    pub fn get_atoms_copy(&self) -> Vec<Atom> { self.atoms.clone() }

    /// Get immutable refernece to the groups in the system.
    pub fn get_groups_as_ref(&self) -> &HashMap<String, Group> { &self.groups }

    /// Get mutable reference to the groups in the system.
    pub unsafe fn get_groups_as_ref_mut(&mut self) -> &mut HashMap<String, Group> { &mut self.groups }

    /// Get copy of the groups in the system.
    pub fn get_groups_copy(&self) -> HashMap<String, Group> { self.groups.clone() }

    /// Get immutable reference to the simulation box.
    pub fn get_box_as_ref(&self) -> &SimBox { &self.simulation_box }

    /// Get mutable reference to the simulation box.
    pub unsafe fn get_box_as_ref_mut(&mut self) -> &mut SimBox { &mut self.simulation_box }

    /// Get copy of the simulation box.
    pub fn get_box_copy(&self) -> SimBox { self.simulation_box.clone() }

    /// Get the number of atoms in the system.
    pub fn get_n_atoms(&self) -> u64 { self.atoms.len() as u64 }

    /// Create two groups each containing all atoms in the system: "all" and "All".
    /// ## Returns
    /// Ok if both groups were created or GroupError in case any group with the same name already exists.
    fn group_create_all(&mut self) -> Result<(), GroupError> {
        self.group_create_from_ranges("all", vec![(0, self.get_n_atoms())])?;
        self.group_create_from_ranges("All", vec![(0, self.get_n_atoms())])?;
        Ok(())
    }

    /// PLACEHOLDER. Make a group with a given name from the given query. PLACEHOLDER.
    pub fn group_create(&mut self, _name: &str, _query: &str) -> Result<(), GroupError> {
        panic!("Not yet implemented.")
    }


    /// Create a group with a given name from the provided atom indices.
    /// ## Returns
    /// Ok if the group was created or GroupError in case a group with the same name already exists.
    /// ## Example
    /// Creating a group "My Group" containing the atoms indexed as 0, 1, 2, 3, 10, 11, and 12.
    /// ```no_run
    /// use groan_rs::system::System;
    /// 
    /// let mut system = System::from_file("system.gro").unwrap();
    /// 
    /// match system.group_create_from_indices("My Group", vec![0, 1, 2, 3, 12, 10, 11]) {
    ///     Err(e) => eprintln!("{}", e),
    ///     Ok(_) => (),
    /// }
    /// ```
    pub fn group_create_from_indices(&mut self, name: &str, atom_indices: Vec<u64>) -> Result<(), GroupError> {
        if self.groups.contains_key(name) {
            return Err(GroupError::AlreadyExists(name.to_string()));
        }

        let group = Group::from_indices(atom_indices, self.get_n_atoms());
        self.groups.insert(name.to_string(), group);
        Ok(())
    }
    

    /// Create a group with a given name from the provided atom ranges.
    /// ## Returns
    /// Ok if the group was created or GroupError in case a group with the same name already exists.
    /// ## Example
    /// Creating a group "My Group" containing the atoms of the system with indices 0-25 and 50-75 (including).
    /// ```no_run
    /// use groan_rs::system::System;
    /// 
    /// let mut system = System::from_file("system.gro").unwrap();
    /// 
    /// match system.group_create_from_ranges("My Group", vec![(0, 25), (50, 75)]) {
    ///     Err(e) => eprintln!("{}", e),
    ///     Ok(_) => (),
    /// }
    /// ```
    pub fn group_create_from_ranges(&mut self, name: &str, atom_ranges: Vec<(u64, u64)>) -> Result<(), GroupError> {
        if self.groups.contains_key(name) {
            return Err(GroupError::AlreadyExists(name.to_string()));
        }
        
        let group = Group::from_ranges(atom_ranges, self.get_n_atoms());
        self.groups.insert(name.to_string(), group);
        Ok(())
    }


    /// Check whether a group with a given name exists in the system.
    pub fn group_exists(&self, name: &str) -> bool {
        self.groups.contains_key(name)
    }


    /// Create an iterator over a group of atoms. The atoms are immutable.
    /// ## Returns
    /// An immutable iterator over the atoms in the group. GroupError in case the group does not exist.
    /// 
    /// ## Example
    /// Printing the atoms of group "Protein".
    /// ```no_run
    /// use groan_rs::system::System;
    /// 
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("Protein", "resid 1 to 29");
    /// 
    /// match system.group_iter("Protein") {
    ///     Ok(iterator) => {
    ///         for atom in iterator {
    ///             println!("{:?}", atom);
    ///         }
    ///     },
    ///     Err(e) => eprintln!("{}", e),
    /// };
    /// ```
    pub fn group_iter(&self, name: &str) -> Result<impl Iterator<Item = &Atom>, GroupError> {
        let group = self.groups.get(name)
        .ok_or(GroupError::NotFound(name.to_string()))?;

        Ok(group.atom_ranges
            .iter()
            .flat_map(move |(start, end)| self.atoms.get(*start as usize ..=*end as usize))
            .flatten()
        )
    }


    /// Create an iterator over a group of atoms. The atoms are mutable.
    /// ## Returns
    /// A mutable iterator over atoms in the group. GroupError in case the group does not exist.
    /// 
    /// ## Example
    /// Translating the atoms of the group "Protein". 
    /// Note that using `group_translate` may be faster.
    /// ```no_run
    /// use groan_rs::system::System;
    /// 
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("Protein", "resid 1 to 29");
    /// let simulation_box = system.get_box_copy();
    /// 
    /// match system.group_iter_mut("Protein") {
    ///     Ok(iterator) => {
    ///         for atom in iterator {
    ///             atom.translate(&[1.0, 2.0, -1.0].into(), &simulation_box);
    ///         }
    ///     },
    ///     Err(e) => eprintln!("{}", e),
    /// };
    /// ```
    pub fn group_iter_mut<'a>(&'a mut self, name: &str) -> Result<impl Iterator<Item = &'a mut Atom>, GroupError> {
        let group = self.groups.get_mut(name)
        .ok_or(GroupError::NotFound(name.to_string()))?;

        unsafe {
            let atoms_ptr = self.atoms.as_mut_ptr();

            Ok(group.atom_ranges
                .iter()
                .flat_map(move |(start, end)| {
                    let start_ptr = atoms_ptr.add(*start as usize);
                    let end_ptr = atoms_ptr.add(*end as usize + 1); // Add 1 to include the end index
                    std::slice::from_raw_parts_mut(start_ptr, end_ptr.offset_from(start_ptr) as usize)
                }))
        }
    }


    /// Create an iterator over all atoms in the system. The atoms are immutable.
    /// ## Example
    /// ```no_run
    /// use groan_rs::system::System;
    /// 
    /// let mut system = System::from_file("system.gro").unwrap();
    /// 
    /// for atom in system.atoms_iter() {
    ///     println!("{:?}", atom);
    /// }
    /// ```
    pub fn atoms_iter(&self) -> impl Iterator<Item = &Atom> {
        self.atoms.iter()
    }


    /// Create an iterator over all atoms in the system. The atoms are mutable.
    /// ## Example
    /// Translating all the atoms in the system by a specified vector. 
    /// Note that using `atoms_translate` may be faster.
    /// ```no_run
    /// use groan_rs::system::System;
    /// 
    /// let mut system = System::from_file("system.gro").unwrap();
    /// let simulation_box = system.get_box_copy();
    /// 
    /// for atom in system.atoms_iter_mut() {
    ///     atom.translate(&[1.0, -1.0, 2.5].into(), &simulation_box);
    /// }
    /// ```
    pub fn atoms_iter_mut(&mut self) -> impl Iterator<Item = &mut Atom> {
        self.atoms.iter_mut()
    }


    /// Translate all atoms of a group by target vector.
    /// ## Returns
    /// Ok or GroupError in case the group does not exist.
    /// ## Example
    /// Translating the atoms of the group "Protein".
    /// ```no_run
    /// use groan_rs::system::System;
    /// 
    /// let mut system = System::from_file("system.gro").unwrap();
    /// system.group_create("Protein", "resid 1 to 29");
    /// 
    /// match system.group_translate("Protein", &[1.0, 2.0, -1.0].into()) {
    ///     Err(e) => eprintln!("{}", e),
    ///     Ok(_) => (),   
    /// }
    /// ```
    pub fn group_translate(&mut self, name: &str, vector: &Vector3D) -> Result<(), GroupError> {
        let group = self.groups.get_mut(name)
        .ok_or(GroupError::NotFound(name.to_string()))?;

        for (start, end) in group.atom_ranges.iter_mut() {
            for i in *start..=*end {
                unsafe { self.atoms.get_unchecked_mut(i as usize).translate(vector, &self.simulation_box) }
            }
        }

        Ok(())
    }

    /// Translate all atoms in the system by target vector.
    /// ## Example
    /// ```no_run
    /// use groan_rs::system::System;
    /// 
    /// let mut system = System::from_file("system.gro").unwrap();
    /// 
    /// system.atoms_translate(&[1.0, 2.0, -1.0].into());
    /// ```
    pub fn atoms_translate(&mut self, vector: &Vector3D) {
        for atom in self.atoms.iter_mut() {
            atom.translate(vector, &self.simulation_box);
        }
    }

    /// Get the number of atoms in target group. 
    /// ## Returns
    /// The number of atoms or GroupError if the group does not exist.
    /// ## Example:
    /// Get the number of atoms in the group "Protein".
    /// ```no_run
    /// use groan_rs::system::System;
    /// 
    /// let mut system = System::from_file("system.gro").unwrap();
    /// 
    /// let n_atoms = system.group_n_atoms("Protein");
    /// ```
    pub fn group_n_atoms(&self, name: &str) -> Result<u64, GroupError> {
        let group = self.groups.get(name)
        .ok_or(GroupError::NotFound(name.to_string()))?;

        Ok(group.get_n_atoms())        
    }


}




#[cfg(test)]
mod tests {
    use super::*;
    use float_cmp::assert_approx_eq;

    #[test]
    fn test_new_system() {

        let system = System::new("System generated using the `groan_rs` library.", Vec::new(), [1.5, 3.3, 0.8].into());

        assert_eq!(system.get_name(), "System generated using the `groan_rs` library.");
        assert_eq!(system.get_atoms_as_ref().len(), 0);

        assert_approx_eq!(f32, system.get_box_as_ref().v1x, 1.5f32);
        assert_approx_eq!(f32, system.get_box_as_ref().v2y, 3.3f32);
        assert_approx_eq!(f32, system.get_box_as_ref().v3z, 0.8f32);
        assert_eq!(system.get_box_as_ref().v1y, 0.0f32);
        assert_eq!(system.get_box_as_ref().v1z, 0.0f32);
        assert_eq!(system.get_box_as_ref().v2x, 0.0f32);
        assert_eq!(system.get_box_as_ref().v2z, 0.0f32);
        assert_eq!(system.get_box_as_ref().v3x, 0.0f32);
        assert_eq!(system.get_box_as_ref().v3y, 0.0f32);

        assert!(system.group_exists("all"));
    }

    #[test]
    fn test_get_n_atoms() {
        let system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(system.get_n_atoms(), 16844);
    }

    #[test]
    fn test_get_name() {
        let system = System::from_file("test_files/example.gro").unwrap();
        assert_eq!(system.get_name(), "INSANE! Membrane UpperLeaflet>POPC=1 LowerLeaflet>POPC=1");
    }

    #[test]
    fn test_get_box_copy() {
        let system = System::from_file("test_files/example_box9.gro").unwrap();

        let mut simbox = system.get_box_copy();

        assert_approx_eq!(f32, simbox.x, system.get_box_as_ref().x);
        assert_approx_eq!(f32, simbox.y, system.get_box_as_ref().y);
        assert_approx_eq!(f32, simbox.z, system.get_box_as_ref().z);

        assert_approx_eq!(f32, simbox.v1y, system.get_box_as_ref().v1y);
        assert_approx_eq!(f32, simbox.v1z, system.get_box_as_ref().v1z);
        assert_approx_eq!(f32, simbox.v2x, system.get_box_as_ref().v2x);

        assert_approx_eq!(f32, simbox.v2z, system.get_box_as_ref().v2z);
        assert_approx_eq!(f32, simbox.v3x, system.get_box_as_ref().v3x);
        assert_approx_eq!(f32, simbox.v3y, system.get_box_as_ref().v3y);

        simbox.v1x = 10.3;
        simbox.v2x = 8.4;
    }

    #[test]
    fn test_get_atoms_copy() {
        let system = System::from_file("test_files/example.gro").unwrap();

        let mut atoms = system.get_atoms_copy();

        for (extracted_atom, system_atom) in atoms.iter().zip(system.get_atoms_as_ref().iter()) {
            assert_eq!(system_atom.get_gmx_atom_number(), extracted_atom.get_gmx_atom_number());
        }

        let _ = atoms.pop();
        assert_eq!(atoms.len(), 16843);
        assert_eq!(system.get_atoms_as_ref().len(), 16844);
    }


    #[test]
    fn test_get_groups_copy() {
        let system = System::from_file("test_files/example_box9.gro").unwrap();

        let mut groups = system.get_groups_copy();

        assert!(groups.contains_key("all"));

        let new_group = Group::from_indices(vec![1, 3, 6, 8], 1000);
        groups.insert("Test".to_string(), new_group);

        assert!(groups.contains_key("Test"));
        assert!(!system.get_groups_as_ref().contains_key("Test"));
    }

    #[test]
    fn test_group_create_from_indices() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.group_create_from_indices("Test", vec![1, 3, 6, 8]).unwrap();

        assert!(system.group_exists("Test"));
        assert_eq!(system.group_n_atoms("Test").unwrap(), 4);

        match system.group_create_from_indices("Test", vec![1, 10, 11, 12, 13]) {
            Err(GroupError::AlreadyExists(e)) => assert_eq!(e, "Test"),
            Ok(_) => panic!("Creating group should have failed, but it was successful."),
            Err(e) => panic!("Failed successfully but incorrect error type `{:?}` was returned.", e),
        }

        assert!(system.group_exists("Test"));
        assert_eq!(system.group_n_atoms("Test").unwrap(), 4);   
    }


    #[test]
    fn test_group_create_from_ranges() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.group_create_from_ranges("Test", vec![(0, 15), (10, 20)]).unwrap();

        assert!(system.group_exists("Test"));
        assert_eq!(system.group_n_atoms("Test").unwrap(), 21);

        match system.group_create_from_ranges("Test", vec![(10, 15)]) {
            Err(GroupError::AlreadyExists(e)) => assert_eq!(e, "Test"),
            Ok(_) => panic!("Creating group should have failed, but it was successful."),
            Err(e) => panic!("Failed successfully but incorrect error type `{:?}` was returned.", e),
        }

        assert!(system.group_exists("Test"));
        assert_eq!(system.group_n_atoms("Test").unwrap(), 21);   
    }


    #[test]
    fn test_group_iter() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        for (group_atom, system_atom) in system.group_iter("Protein").unwrap().zip(system.get_atoms_as_ref().iter().take(61)) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
        }   
    }

    #[test]
    fn test_group_iter_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();
        system.read_ndx("test_files/index.ndx").unwrap();

        let extracted = system.get_atoms_copy();

        for (group_atom, system_atom) in system.group_iter_mut("Protein").unwrap().zip(extracted.iter().take(61)) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
            
            group_atom.translate_nopbc(&Vector3D::from([0.5, -1.1, 2.4]));
            assert_approx_eq!(f32, group_atom.get_position().x, system_atom.get_position().x + 0.5);
            assert_approx_eq!(f32, group_atom.get_position().y, system_atom.get_position().y - 1.1);
            assert_approx_eq!(f32, group_atom.get_position().z, system_atom.get_position().z + 2.4);
        }   
    }

    #[test]
    fn test_atoms_iter() {
        let system = System::from_file("test_files/example.gro").unwrap();

        let atoms = system.get_atoms_as_ref();

        for (extracted_atom, system_atom) in atoms.iter().zip(system.atoms_iter()) {
            assert_eq!(system_atom.get_gmx_atom_number(), extracted_atom.get_gmx_atom_number());
        }
    }

    #[test]
    fn test_atoms_iter_mut() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        let extracted = system.get_atoms_copy();

        for (group_atom, system_atom) in system.atoms_iter_mut().zip(extracted.iter()) {
            assert_eq!(system_atom.get_gmx_atom_number(), group_atom.get_gmx_atom_number());
            
            group_atom.translate_nopbc(&Vector3D::from([0.5, -1.1, 2.4]));
            assert_approx_eq!(f32, group_atom.get_position().x, system_atom.get_position().x + 0.5);
            assert_approx_eq!(f32, group_atom.get_position().y, system_atom.get_position().y - 1.1);
            assert_approx_eq!(f32, group_atom.get_position().z, system_atom.get_position().z + 2.4);
        }   
    }

    #[test]
    fn test_atoms_translate() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.atoms_translate(&Vector3D::from([3.5, -1.1, 5.4]));

        let first = system.atoms_iter().next().unwrap();
        let last = system.atoms_iter().last().unwrap();

        let first_pos = first.get_position();
        let last_pos = last.get_position();

        assert_approx_eq!(f32, first_pos.x, 12.997);
        assert_approx_eq!(f32, first_pos.y, 0.889);
        assert_approx_eq!(f32, first_pos.z, 1.64453);

        assert_approx_eq!(f32, last_pos.x, 12.329);
        assert_approx_eq!(f32, last_pos.y, 10.086);
        assert_approx_eq!(f32, last_pos.z, 7.475);
    }

    #[test]
    fn test_group_translate() {
        let mut system = System::from_file("test_files/example.gro").unwrap();

        system.group_translate("all", &Vector3D::from([3.5, -1.1, 5.4])).unwrap();

        let first = system.atoms_iter().next().unwrap();
        let last = system.atoms_iter().last().unwrap();

        let first_pos = first.get_position();
        let last_pos = last.get_position();

        assert_approx_eq!(f32, first_pos.x, 12.997);
        assert_approx_eq!(f32, first_pos.y, 0.889);
        assert_approx_eq!(f32, first_pos.z, 1.64453);

        assert_approx_eq!(f32, last_pos.x, 12.329);
        assert_approx_eq!(f32, last_pos.y, 10.086);
        assert_approx_eq!(f32, last_pos.z, 7.475);
    }
}