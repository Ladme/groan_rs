// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Implementation of methods for identification of hydrogen bonds.

use getset::CopyGetters;
use hashbrown::HashSet;
use indexmap::IndexMap;

use crate::{
    errors::{AtomError, CellGridError, HBondError, PositionError, SimBoxError},
    prelude::{CellGrid, CellNeighbors, Dimension, SimBox, Vector3D},
    structures::{
        group::Group,
        traj_convert::{ConvertableTrajRead, FrameAnalyze, TrajAnalyzer},
    },
};

use super::System;

/// Structure defining a chain for the hydrogen bond analysis.
/// Construct using [`HBondChain::new`].
#[derive(Debug, Clone)]
pub struct HBondChain {
    acceptors: String,
    donors: String,
    hydrogens: String,
}

impl HBondChain {
    /// Define a chain for the hydrogen bond analysis.
    ///
    /// ## Parameters
    /// - `acceptors`: groan selection language (GSL) query specifying atoms that can be acceptors of hydrogen bonds
    /// - `donors`: GSL query specifying atoms that can be donors of hydrogen bonds
    /// - `hydrogens`: GSL query specifying hydrogen atoms of the chain
    ///
    /// During the analysis, only hydrogens bound to donor atoms
    /// and only donor atoms with at least one hydrogen will be considered.
    /// You don't have to select them manually.
    pub fn new(acceptors: &str, donors: &str, hydrogens: &str) -> Self {
        Self {
            acceptors: acceptors.to_owned(),
            donors: donors.to_owned(),
            hydrogens: hydrogens.to_owned(),
        }
    }
}

/// Structure describing the identified hydrogen bond.
#[derive(Debug, Clone, CopyGetters)]
pub struct HBond {
    /// Index of the donor atom. Atoms are indexed starting from 0.
    #[getset(get_copy = "pub")]
    donor: usize,
    /// Index of the hydrogen.
    #[getset(get_copy = "pub")]
    hydrogen: usize,
    /// Index of the acceptor atom.
    #[getset(get_copy = "pub")]
    acceptor: usize,
    /// Distance between donor and acceptor (in nm).
    #[getset(get_copy = "pub")]
    distance: f32,
    /// Angle between the donor, hydrogen, and acceptor (in degrees).
    #[getset(get_copy = "pub")]
    angle: f32,
}

impl HBond {
    /// Creates a new `HBond` instance.
    #[allow(dead_code)]
    fn new(donor: usize, hydrogen: usize, acceptor: usize, distance: f32, angle: f32) -> Self {
        HBond {
            donor,
            hydrogen,
            acceptor,
            distance,
            angle,
        }
    }
}

/// Structure storing information needed for the hydrogen bond analysis.
#[derive(Debug, Clone)]
pub struct HBondAnalysis {
    chains: Vec<HBondChainGroups>,
    pairs: Vec<(usize, usize)>,
    max_distance: f32,
    min_angle: f32,
}

/// Structure storing acceptor and donor atoms for a single chain.
#[derive(Debug, Clone)]
struct HBondChainGroups {
    acceptors: Group,
    donors: Vec<Donor>,
}

/// Indices of a single donor atom and its bonded hydrogens.
#[derive(Debug, Clone)]
struct Donor {
    donor: usize,
    hydrogens: Vec<usize>,
}

impl HBondChainGroups {
    /// Create new chains for the hydrogen bond analysis.
    fn new(
        system: &System,
        acceptors: &str,
        donors: &str,
        hydrogens: &str,
    ) -> Result<Self, HBondError> {
        let acc_group = Group::from_query(acceptors, system).map_err(HBondError::SelectError)?;
        let don_group = Group::from_query(donors, system).map_err(HBondError::SelectError)?;
        let hyd_group = Group::from_query(hydrogens, system).map_err(HBondError::SelectError)?;

        let mut donor_structs = Vec::new();
        for donor in don_group.get_atoms().iter() {
            // get hydrogens bonded to the donor atom
            let mut donor_hydrogens = Vec::new();
            for bonded in system.bonded_atoms_iter(donor).unwrap() {
                if hyd_group.get_atoms().isin(bonded.get_index()) {
                    donor_hydrogens.push(bonded.get_index());
                }
            }

            // ignore donor atoms with no bonded hydrogens
            if donor_hydrogens.is_empty() {
                continue;
            }

            donor_structs.push(Donor {
                donor,
                hydrogens: donor_hydrogens,
            })
        }

        if acc_group.get_n_atoms() == 0 && donor_structs.is_empty() {
            return Err(HBondError::EmptyChain);
        }

        Ok(HBondChainGroups {
            acceptors: acc_group,
            donors: donor_structs,
        })
    }
}

/// Identified hydrogen bonds for each pair of chains.
pub type HBondMap = IndexMap<(usize, usize), Vec<HBond>>;

impl FrameAnalyze for HBondAnalysis {
    /// Hydrogen bonds identified for each analyzed pair of chains.
    type AnalysisResult = HBondMap;
    type Error = HBondError;

    /// Find hydrogen bonds in a single frame.
    fn analyze(&mut self, system: &System) -> Result<HBondMap, HBondError> {
        // assign acceptors to cell grids to speed up donor-acceptor searches
        // safety: groups are constructed for the given system
        let acc_cell_grids = self
            .chains
            .iter()
            .map(|chain| unsafe {
                CellGrid::new_from_group(system, &chain.acceptors, self.max_distance)
            })
            .collect::<Result<Vec<_>, CellGridError>>()
            .map_err(HBondError::CellGridError)?;

        // find hydrogen bonds between the specified pairs of chains
        let mut all_hbonds = IndexMap::new();
        for (chain1, chain2) in self.pairs.iter() {
            let hbonds = if chain1 == chain2 {
                HBondAnalysis::analyze_single(
                    system,
                    &acc_cell_grids[*chain1],
                    &self.chains[*chain1].donors,
                    self.max_distance,
                    self.min_angle,
                )?
            } else {
                HBondAnalysis::analyze_pair(
                    system,
                    &acc_cell_grids[*chain1],
                    &acc_cell_grids[*chain2],
                    &self.chains[*chain1].donors,
                    &self.chains[*chain2].donors,
                    self.max_distance,
                    self.min_angle,
                )?
            };

            if all_hbonds.insert((*chain1, *chain2), hbonds).is_some() {
                panic!("FATAL GROAN ERROR | HBondAnalysis::analyze | Hydrogen bonds for pair {}x{} were calculated multiple times.", chain1, chain2);
            }
        }

        Ok(all_hbonds)
    }
}

impl HBondAnalysis {
    /// Find hydrogen bonds between two chains.
    #[inline]
    fn analyze_pair<'a>(
        system: &'a System,
        acc_grid1: &CellGrid<'a>,
        acc_grid2: &CellGrid<'a>,
        donors1: &[Donor],
        donors2: &[Donor],
        max_distance: f32,
        min_angle: f32,
    ) -> Result<Vec<HBond>, HBondError> {
        let mut bonds = Self::analyze_single(system, acc_grid1, donors2, max_distance, min_angle)?;
        bonds.extend(Self::analyze_single(
            system,
            acc_grid2,
            donors1,
            max_distance,
            min_angle,
        )?);

        Ok(bonds)
    }

    /// Find hydrogen bonds inside a single chain or between acceptors of chain A and donors of chain B.
    fn analyze_single<'a>(
        system: &'a System,
        acc_grid: &CellGrid<'a>,
        donors: &[Donor],
        max_distance: f32,
        min_angle: f32,
    ) -> Result<Vec<HBond>, HBondError> {
        let simbox = system
            .get_box()
            .ok_or(HBondError::InvalidSimBox(SimBoxError::DoesNotExist))?;

        let mut hbonds = Vec::new();
        for donor in donors.iter() {
            // safety: donor index must be valid
            let donor_atom = unsafe { system.get_atom_unchecked(donor.donor) };
            let donor_position = donor_atom.get_position().ok_or_else(|| {
                HBondError::AtomError(AtomError::InvalidPosition(PositionError::NoPosition(
                    donor.donor,
                )))
            })?;

            for acceptor in acc_grid
                .neighbors_iter(donor_position.clone(), CellNeighbors::default())
                .filter(|a| a.get_index() != donor.donor)
            {
                let acceptor_position = acceptor.get_position().ok_or_else(|| {
                    HBondError::AtomError(AtomError::InvalidPosition(PositionError::NoPosition(
                        acceptor.get_index(),
                    )))
                })?;

                // donor-acceptor distance must be lower than `max_distance`
                let distance = acceptor_position.distance(donor_position, Dimension::XYZ, simbox);
                if distance > max_distance {
                    continue;
                }

                for hydrogen in donor.hydrogens.iter() {
                    // safety: hydrogen index must be valid
                    let hydrogen_atom = unsafe { system.get_atom_unchecked(*hydrogen) };
                    let hydrogen_position = hydrogen_atom.get_position().ok_or_else(|| {
                        HBondError::AtomError(AtomError::InvalidPosition(
                            PositionError::NoPosition(*hydrogen),
                        ))
                    })?;

                    // donor-hydrogen-acceptor angle must be higher than `min_angle`
                    let angle = Self::calc_angle(
                        donor_position,
                        hydrogen_position,
                        acceptor_position,
                        simbox,
                    );
                    if angle < min_angle {
                        continue;
                    }

                    hbonds.push(HBond {
                        donor: donor.donor,
                        hydrogen: *hydrogen,
                        acceptor: acceptor.get_index(),
                        distance,
                        angle,
                    });
                }
            }
        }

        Ok(hbonds)
    }

    /// Calculate the angle between donor, hydrogen, and acceptor.
    #[inline]
    fn calc_angle(
        donor: &Vector3D,
        hydrogen: &Vector3D,
        acceptor: &Vector3D,
        simbox: &SimBox,
    ) -> f32 {
        let hd_vec = hydrogen.vector_to(donor, simbox);
        let ha_vec = hydrogen.vector_to(acceptor, simbox);
        let angle = hd_vec.angle(&ha_vec).to_degrees();
        if angle.is_nan() {
            // if angle is `nan`, the angle can be either 0° or 180°; we need to check the position of the hydrogen
            Self::handle_nan(donor, hydrogen, acceptor, simbox)
        } else {
            angle
        }
    }

    /// Handle the situation where the angle is `nan`.
    #[cold]
    fn handle_nan(
        donor: &Vector3D,
        hydrogen: &Vector3D,
        acceptor: &Vector3D,
        simbox: &SimBox,
    ) -> f32 {
        // if the hydrogen is closer to the acceptor than the donor is to the acceptor, then the angle is 180°
        // and we have a hydrogen bond
        if hydrogen.distance(acceptor, Dimension::XYZ, simbox)
            < donor.distance(acceptor, Dimension::XYZ, simbox)
        {
            180.0
        // otherwise the hydrogen is actually oriented away from the acceptor, and we have no hydrogen bond
        } else {
            0.0
        }
    }

    /// Check that all chains are part of at least on pair and that all chains specified in pairs exist.
    fn sanity_check_pairs(pairs: &[(usize, usize)], n_chains: usize) -> Result<(), HBondError> {
        let mut pairs_set = HashSet::new();
        let mut used_chains = HashSet::new();

        for (chain1, chain2) in pairs.iter() {
            for ch in [chain1, chain2] {
                if *ch >= n_chains {
                    return Err(HBondError::NonexistentChain(*ch));
                }
            }

            if chain1 != chain2 {
                if !pairs_set.insert((*chain1, *chain2)) || !pairs_set.insert((*chain2, *chain1)) {
                    return Err(HBondError::PairSpecifiedMultipleTimes(*chain1, *chain2));
                }
            } else if !pairs_set.insert((*chain1, *chain2)) {
                return Err(HBondError::PairSpecifiedMultipleTimes(*chain1, *chain2));
            }

            used_chains.insert(*chain1);
            used_chains.insert(*chain2);
        }

        if used_chains.len() != n_chains {
            return Err(HBondError::UnusedChain);
        }

        Ok(())
    }
}

/// Trajectory analyzer-iterator. Reads a trajectory file returning the current frame and the hydrogen bonds identified in this frame.
pub type HBondIterator<'a, Reader> = TrajAnalyzer<'a, Reader, HBondAnalysis>;

pub trait HBondTrajRead<'a>: ConvertableTrajRead<'a> {
    /// Find hydrogen bonds between specified chains for each frame of the trajectory.
    ///
    /// ## Parameters
    /// - `chains`: Chains to work with; see [`HBondChain::new`] for more information.
    /// - `pairs`: Pairs of chains between which hydrogen bonds should be searched.
    /// - `max_distance`: Maximum distance (in nm) between a donor and an acceptor atom for the hydrogen bond
    ///   to be considered (recommended value is 0.3 nm).
    /// - `min_angle`: Minimum angle (in degrees) between donor-hydrogen-acceptor for the hydrogen bond
    ///   to be considered (recommended value is 150°).
    ///
    /// ## Returns
    /// - Returns an `HBondIterator`, which is an iterable trajectory analyzer.
    ///   Each iteration yields the current unmodified frame of the trajectory along with the
    ///   identified hydrogen bonds for the individual pairs of chains.
    /// - In case of an error, returns an `HBondError`.
    ///
    /// ## Example
    /// We have a simulation of four protein chains that form a tetramer. We want to find
    /// hydrogen bonds between chains A and B, B and C, C and D, and D and A.
    /// We are NOT interested in hydrogen bonds between chains A and C or B and D.
    /// However, we are also interested in hydrogen bonds within chains A and C.
    ///
    /// ```no_run
    /// # #[cfg(any(feature = "molly", not(feature = "no-xdrfile")))] {
    /// # use groan_rs::prelude::*;
    /// #
    /// // Returns hydrogen bonds for each frame of the analyzed trajectory.
    /// fn find_hbonds() -> Result<Vec<HBondMap>, Box<dyn std::error::Error + Send + Sync>> {
    ///     // the PDB file must contain connectivity information; you can also use a TPR file
    ///     let mut system = System::from_file("system.pdb")?;
    ///     system.guess_elements(Elements::default())?;
    ///
    ///     let mut hbonds = Vec::new();
    ///
    ///     for result in system.xtc_iter("md.xtc")?.hbonds_analyze(
    ///         // acceptors, donors, and hydrogens of the individual chains
    ///         // only hydrogens bonded to donors and donors with at least one hydrogen will be considered
    ///         vec![
    ///             HBondChain::new(
    ///                 "chain A and element symbol O N",
    ///                 "chain A and element symbol O N",
    ///                 "chain A and element symbol H",
    ///             ),
    ///             HBondChain::new(
    ///                 "chain B and element symbol O N",
    ///                 "chain B and element symbol O N",
    ///                 "chain B and element symbol H",
    ///             ),
    ///             HBondChain::new(
    ///                 "chain C and element symbol O N",
    ///                 "chain C and element symbol O N",
    ///                 "chain C and element symbol H",
    ///             ),
    ///             HBondChain::new(
    ///                 "chain D and element symbol O N",
    ///                 "chain D and element symbol O N",
    ///                 "chain D and element symbol H",
    ///             ),
    ///         ],
    ///         // search for hydrogen bonds between the following pairs of chains: A-B, B-C, C-D, D-A, A-A, C-C
    ///         vec![(0, 1), (1, 2), (2, 3), (3, 0), (0, 0), (2, 2)],
    ///         // maximal donor-acceptor distance (in nm)
    ///         0.3,
    ///         // minimal donor-hydrogen-acceptor angle (in degrees)
    ///         150.0,
    ///     )? {
    ///         let (_, hbonds_for_frame) = result?;
    ///         hbonds.push(hbonds_for_frame);
    ///     }
    ///
    ///     Ok(hbonds)
    /// }
    /// # }
    /// ```
    ///
    /// ## Notes
    /// - Only orthogonal simulation boxes are supported!
    /// - Selections of acceptor/donors should not overlap with the selection for hydrogens!
    ///   Acceptors can overlap with donors.
    /// - The system must contain information about connectivity between atoms (i.e., bonds).
    /// - All specified chains must be used; otherwise, an error is returned.
    /// - Calculation for any pair of chains cannot be requested multiple times.
    fn hbonds_analyze(
        mut self,
        chains: Vec<HBondChain>,
        pairs: Vec<(usize, usize)>,
        max_distance: f32,
        min_angle: f32,
    ) -> Result<HBondIterator<'a, Self>, HBondError> {
        let system = unsafe { &*self.get_system() };

        let chains_groups = chains
            .into_iter()
            .map(|chain| {
                HBondChainGroups::new(system, &chain.acceptors, &chain.donors, &chain.hydrogens)
            })
            .collect::<Result<Vec<_>, HBondError>>()?;

        HBondAnalysis::sanity_check_pairs(&pairs, chains_groups.len())?;

        let analysis = HBondAnalysis {
            chains: chains_groups,
            pairs,
            max_distance,
            min_angle,
        };

        Ok(self.analyze(analysis))
    }
}

/// Blanket implementation of `HBondTrajRead` trait for all convertable trajectory readers.
impl<'a, T> HBondTrajRead<'a> for T where T: ConvertableTrajRead<'a> {}

#[cfg(test)]
mod tests {
    use float_cmp::assert_approx_eq;

    use super::*;

    fn compare_hbonds(bond1: &HBond, bond2: &HBond) {
        assert_eq!(bond1.donor, bond2.donor);
        assert_eq!(bond1.hydrogen, bond2.hydrogen);
        assert_eq!(bond1.acceptor, bond2.acceptor);
        assert_approx_eq!(f32, bond1.distance, bond2.distance, epsilon = 1e-3);
        assert_approx_eq!(f32, bond1.angle, bond2.angle, epsilon = 1e-3);
    }

    #[test]
    #[cfg(any(feature = "molly", not(feature = "no-xdrfile")))]
    fn test_hbonds_analyze_simple_water() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.tpr").unwrap();

        let expected_n_hbonds = [
            4675, 4644, 4629, 4617, 4651, 4532, 4649, 4611, 4621, 4701, 4694, 4650, 4565, 4681,
            4699, 4711, 4652, 4649, 4697, 4652, 4644,
        ];
        let expected_hbonds = [
            HBond::new(17527, 17528, 21100, 0.262, 157.241),
            HBond::new(32782, 32784, 22366, 0.287, 164.069),
            HBond::new(17515, 17516, 19948, 0.283, 163.444),
            HBond::new(32785, 32786, 26968, 0.288, 154.899),
            HBond::new(17521, 17522, 22327, 0.290, 157.527),
            HBond::new(32785, 32786, 24001, 0.273, 157.663),
            HBond::new(17515, 17516, 18742, 0.271, 168.506),
            HBond::new(32785, 32787, 20182, 0.296, 164.179),
            HBond::new(17515, 17516, 30169, 0.278, 175.880),
            HBond::new(32785, 32787, 20692, 0.299, 171.338),
            HBond::new(17521, 17523, 27085, 0.256, 167.116),
            HBond::new(32785, 32787, 26698, 0.272, 165.705),
            HBond::new(17515, 17516, 19909, 0.282, 164.541),
            HBond::new(32785, 32787, 23671, 0.257, 177.534),
            HBond::new(17518, 17519, 27838, 0.276, 173.045),
            HBond::new(32779, 32780, 18964, 0.290, 158.859),
            HBond::new(17515, 17517, 19708, 0.271, 158.164),
            HBond::new(32785, 32787, 23878, 0.262, 174.790),
            HBond::new(17515, 17516, 32596, 0.275, 159.038),
            HBond::new(32782, 32783, 30319, 0.291, 167.926),
            HBond::new(17521, 17523, 21157, 0.277, 174.121),
            HBond::new(32782, 32784, 19798, 0.268, 151.393),
            HBond::new(17515, 17516, 30709, 0.281, 162.790),
            HBond::new(32785, 32787, 29284, 0.266, 156.265),
            HBond::new(17518, 17520, 28738, 0.295, 153.376),
            HBond::new(32785, 32786, 24304, 0.277, 154.977),
            HBond::new(17515, 17517, 19168, 0.266, 155.122),
            HBond::new(32785, 32786, 19582, 0.262, 169.901),
            HBond::new(17515, 17516, 22936, 0.279, 157.140),
            HBond::new(32785, 32787, 19150, 0.267, 176.785),
            HBond::new(17515, 17516, 28546, 0.269, 157.875),
            HBond::new(32785, 32786, 25645, 0.279, 162.166),
            HBond::new(17515, 17516, 30748, 0.271, 170.012),
            HBond::new(32785, 32786, 19261, 0.267, 155.783),
            HBond::new(17515, 17517, 32479, 0.296, 171.556),
            HBond::new(32785, 32787, 23401, 0.292, 165.894),
            HBond::new(17515, 17517, 21205, 0.249, 174.105),
            HBond::new(32779, 32780, 25045, 0.281, 153.504),
            HBond::new(17518, 17519, 26014, 0.296, 176.495),
            HBond::new(32785, 32787, 22585, 0.266, 156.385),
            HBond::new(17530, 17532, 32503, 0.269, 171.010),
            HBond::new(32785, 32786, 31561, 0.269, 168.088),
        ];
        for (frame, result) in system
            .xtc_iter("test_files/aa_membrane_peptide.xtc")
            .unwrap()
            .hbonds_analyze(
                vec![HBondChain::new(
                    "resname SOL and name OW",
                    "resname SOL and name OW",
                    "resname SOL and name HW1 HW2",
                )],
                vec![(0, 0)],
                0.3,
                150.0,
            )
            .unwrap()
            .enumerate()
        {
            let (_, hbonds_frame) = result.unwrap();
            assert_eq!(hbonds_frame.len(), 1);
            assert_eq!(
                hbonds_frame.get(&(0, 0)).unwrap().len(),
                expected_n_hbonds[frame]
            );

            let first_hbond = hbonds_frame.get(&(0, 0)).unwrap().first().unwrap();
            let last_hbond = hbonds_frame.get(&(0, 0)).unwrap().last().unwrap();

            compare_hbonds(first_hbond, &expected_hbonds[frame * 2]);
            compare_hbonds(last_hbond, &expected_hbonds[frame * 2 + 1]);
        }
    }

    #[test]
    #[cfg(any(feature = "molly", not(feature = "no-xdrfile")))]
    fn test_hbonds_analyze_water_various_parameters() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.tpr").unwrap();

        for (distance, angle) in [0.1, 0.2, 0.3, 0.5]
            .into_iter()
            .zip([150.0, 130.0, 160.0, 160.0])
        {
            for result in system
                .xtc_iter("test_files/aa_membrane_peptide.xtc")
                .unwrap()
                .hbonds_analyze(
                    vec![HBondChain::new(
                        "resname SOL and name OW",
                        "resname SOL and name OW",
                        "resname SOL and name HW1 HW2",
                    )],
                    vec![(0, 0)],
                    distance,
                    angle,
                )
                .unwrap()
            {
                let (frame, hbonds_frame) = result.unwrap();
                assert_eq!(hbonds_frame.len(), 1);

                for bond in hbonds_frame.get(&(0, 0)).unwrap().iter() {
                    assert!(bond.distance <= distance);
                    assert_approx_eq!(
                        f32,
                        bond.distance,
                        frame
                            .atoms_distance(bond.acceptor, bond.donor, Dimension::XYZ)
                            .unwrap()
                    );

                    let hydrogen = frame.get_atom(bond.hydrogen).unwrap();
                    let donor = frame.get_atom(bond.donor).unwrap();
                    let acceptor = frame.get_atom(bond.acceptor).unwrap();

                    let hd_vec = hydrogen
                        .get_position()
                        .unwrap()
                        .vector_to(donor.get_position().unwrap(), frame.get_box().unwrap());
                    let ha_vec = hydrogen
                        .get_position()
                        .unwrap()
                        .vector_to(acceptor.get_position().unwrap(), frame.get_box().unwrap());

                    let expected_angle = hd_vec.angle(&ha_vec).to_degrees();
                    if expected_angle.is_nan() {
                        continue;
                    }

                    assert!(bond.angle >= angle);
                    assert_approx_eq!(f32, bond.angle, expected_angle);
                }
            }
        }
    }

    #[test]
    #[cfg(any(feature = "molly", not(feature = "no-xdrfile")))]
    fn test_hbonds_analyze_simple_protein() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.tpr").unwrap();

        let mut expected_hbonds = [
            HBond::new(118, 119, 61, 0.277, 158.384),
            HBond::new(129, 130, 72, 0.299, 155.371),
            HBond::new(193, 194, 132, 0.286, 164.124),
            HBond::new(212, 213, 151, 0.287, 168.456),
            HBond::new(238, 239, 170, 0.282, 158.546),
            HBond::new(291, 292, 241, 0.297, 162.274),
            HBond::new(309, 310, 252, 0.279, 170.511),
            HBond::new(358, 359, 301, 0.263, 172.024),
            HBond::new(92, 93, 50, 0.290, 153.003),
            HBond::new(118, 119, 61, 0.265, 174.305),
            HBond::new(133, 134, 91, 0.280, 157.538),
            HBond::new(171, 172, 121, 0.296, 150.854),
            HBond::new(231, 232, 170, 0.284, 153.273),
            HBond::new(238, 239, 170, 0.265, 172.227),
            HBond::new(272, 273, 230, 0.293, 163.410),
            HBond::new(291, 292, 241, 0.299, 165.239),
            HBond::new(309, 310, 252, 0.274, 162.477),
            HBond::new(358, 359, 301, 0.277, 154.927),
            HBond::new(92, 93, 50, 0.285, 163.483),
            HBond::new(118, 119, 61, 0.262, 161.571),
            HBond::new(212, 213, 151, 0.282, 157.311),
            HBond::new(231, 232, 170, 0.293, 164.043),
            HBond::new(242, 243, 192, 0.299, 168.391),
            HBond::new(272, 273, 230, 0.281, 154.472),
            HBond::new(291, 292, 241, 0.279, 174.785),
            HBond::new(313, 314, 271, 0.296, 170.105),
            HBond::new(332, 333, 290, 0.282, 168.727),
            HBond::new(358, 359, 301, 0.264, 157.403),
            HBond::new(62, 63, 12, 0.293, 156.482),
            HBond::new(92, 93, 50, 0.284, 158.151),
            HBond::new(129, 130, 72, 0.250, 165.311),
            HBond::new(152, 153, 110, 0.278, 152.059),
            HBond::new(212, 213, 151, 0.282, 177.401),
            HBond::new(231, 232, 170, 0.295, 164.278),
            HBond::new(238, 239, 170, 0.279, 161.186),
            HBond::new(272, 273, 230, 0.294, 162.935),
            HBond::new(291, 292, 241, 0.291, 153.813),
            HBond::new(309, 310, 252, 0.293, 157.863),
            HBond::new(313, 314, 271, 0.279, 160.410),
            HBond::new(358, 359, 301, 0.262, 178.920),
            HBond::new(62, 63, 12, 0.281, 156.160),
            HBond::new(92, 93, 50, 0.286, 154.766),
            HBond::new(118, 119, 61, 0.277, 162.976),
            HBond::new(212, 213, 151, 0.281, 154.649),
            HBond::new(238, 239, 170, 0.264, 177.380),
            HBond::new(242, 243, 192, 0.284, 150.486),
            HBond::new(272, 273, 230, 0.293, 173.451),
            HBond::new(309, 310, 252, 0.268, 164.417),
            HBond::new(118, 119, 61, 0.249, 169.915),
            HBond::new(133, 134, 91, 0.292, 151.596),
            HBond::new(152, 153, 110, 0.286, 152.724),
            HBond::new(212, 213, 151, 0.289, 167.275),
            HBond::new(238, 239, 170, 0.280, 156.094),
            HBond::new(291, 292, 241, 0.297, 173.762),
            HBond::new(309, 310, 252, 0.300, 160.401),
            HBond::new(358, 359, 301, 0.266, 152.983),
            HBond::new(69, 70, 12, 0.300, 165.175),
            HBond::new(92, 93, 50, 0.280, 173.466),
            HBond::new(118, 119, 61, 0.286, 159.107),
            HBond::new(212, 213, 151, 0.296, 165.331),
            HBond::new(238, 239, 170, 0.266, 178.577),
            HBond::new(242, 243, 192, 0.288, 154.414),
            HBond::new(272, 273, 230, 0.268, 165.085),
            HBond::new(291, 292, 241, 0.271, 153.922),
            HBond::new(309, 310, 252, 0.269, 175.224),
            HBond::new(358, 359, 301, 0.279, 162.198),
            HBond::new(62, 63, 12, 0.287, 164.914),
            HBond::new(129, 130, 72, 0.278, 157.190),
            HBond::new(212, 213, 151, 0.293, 164.561),
            HBond::new(238, 239, 170, 0.253, 168.866),
            HBond::new(298, 299, 241, 0.267, 160.681),
            HBond::new(309, 310, 252, 0.289, 173.380),
            HBond::new(332, 333, 290, 0.286, 178.049),
            HBond::new(351, 352, 331, 0.291, 154.856),
            HBond::new(92, 93, 50, 0.283, 161.852),
            HBond::new(118, 119, 61, 0.270, 175.245),
            HBond::new(129, 130, 72, 0.299, 171.389),
            HBond::new(171, 172, 121, 0.283, 162.073),
            HBond::new(212, 213, 151, 0.289, 166.929),
            HBond::new(238, 239, 170, 0.275, 161.445),
            HBond::new(242, 243, 192, 0.272, 151.658),
            HBond::new(129, 130, 72, 0.276, 171.190),
            HBond::new(238, 239, 170, 0.282, 154.474),
            HBond::new(242, 243, 192, 0.295, 169.167),
            HBond::new(291, 292, 241, 0.275, 162.751),
            HBond::new(118, 119, 61, 0.275, 156.837),
            HBond::new(122, 123, 72, 0.291, 158.820),
            HBond::new(129, 130, 72, 0.269, 165.649),
            HBond::new(152, 153, 110, 0.295, 159.428),
            HBond::new(171, 172, 121, 0.289, 166.658),
            HBond::new(231, 232, 170, 0.296, 151.064),
            HBond::new(238, 239, 170, 0.274, 168.835),
            HBond::new(272, 273, 230, 0.296, 158.132),
            HBond::new(291, 292, 241, 0.284, 166.551),
            HBond::new(309, 310, 252, 0.286, 164.124),
            HBond::new(118, 119, 61, 0.276, 154.467),
            HBond::new(122, 123, 72, 0.287, 153.002),
            HBond::new(129, 130, 72, 0.274, 172.460),
            HBond::new(152, 153, 110, 0.282, 162.243),
            HBond::new(212, 213, 151, 0.271, 154.546),
            HBond::new(238, 239, 170, 0.275, 169.254),
            HBond::new(69, 70, 12, 0.278, 152.624),
            HBond::new(118, 119, 61, 0.274, 155.652),
            HBond::new(152, 153, 110, 0.291, 169.426),
            HBond::new(212, 213, 151, 0.287, 176.696),
            HBond::new(242, 243, 192, 0.282, 163.440),
            HBond::new(309, 310, 252, 0.275, 164.452),
            HBond::new(62, 63, 12, 0.281, 158.566),
            HBond::new(118, 119, 61, 0.274, 168.244),
            HBond::new(129, 130, 72, 0.284, 156.707),
            HBond::new(152, 153, 110, 0.293, 154.077),
            HBond::new(212, 213, 151, 0.292, 160.418),
            HBond::new(238, 239, 170, 0.272, 155.343),
            HBond::new(309, 310, 252, 0.282, 159.194),
            HBond::new(358, 359, 301, 0.267, 177.397),
            HBond::new(118, 119, 61, 0.276, 154.921),
            HBond::new(129, 130, 72, 0.300, 165.732),
            HBond::new(193, 194, 132, 0.300, 151.595),
            HBond::new(238, 239, 170, 0.285, 173.748),
            HBond::new(272, 273, 230, 0.294, 157.629),
            HBond::new(291, 292, 241, 0.295, 163.543),
            HBond::new(309, 310, 252, 0.251, 169.897),
            HBond::new(69, 70, 12, 0.264, 169.224),
            HBond::new(111, 112, 61, 0.291, 150.213),
            HBond::new(118, 119, 61, 0.262, 158.168),
            HBond::new(129, 130, 72, 0.287, 160.073),
            HBond::new(152, 153, 110, 0.284, 171.070),
            HBond::new(212, 213, 151, 0.285, 170.801),
            HBond::new(238, 239, 170, 0.282, 171.452),
            HBond::new(242, 243, 192, 0.297, 158.289),
            HBond::new(272, 273, 230, 0.285, 156.213),
            HBond::new(291, 292, 241, 0.291, 165.486),
            HBond::new(309, 310, 252, 0.269, 167.058),
            HBond::new(358, 359, 301, 0.272, 152.840),
            HBond::new(92, 93, 50, 0.282, 168.961),
            HBond::new(118, 119, 61, 0.242, 152.040),
            HBond::new(129, 130, 72, 0.269, 154.648),
            HBond::new(238, 239, 170, 0.286, 170.543),
            HBond::new(272, 273, 230, 0.300, 152.106),
            HBond::new(309, 310, 252, 0.294, 170.559),
            HBond::new(358, 359, 301, 0.279, 151.287),
            HBond::new(62, 63, 12, 0.290, 151.161),
            HBond::new(92, 93, 50, 0.283, 167.106),
            HBond::new(111, 112, 61, 0.293, 154.804),
            HBond::new(129, 130, 72, 0.290, 170.316),
            HBond::new(133, 134, 91, 0.278, 165.672),
            HBond::new(171, 172, 121, 0.289, 150.103),
            HBond::new(212, 213, 151, 0.293, 176.738),
            HBond::new(238, 239, 170, 0.257, 159.220),
            HBond::new(242, 243, 192, 0.294, 152.489),
            HBond::new(253, 254, 211, 0.295, 155.806),
            HBond::new(291, 292, 241, 0.297, 160.017),
            HBond::new(309, 310, 252, 0.249, 169.495),
            HBond::new(152, 153, 110, 0.292, 152.904),
            HBond::new(171, 172, 121, 0.296, 151.920),
            HBond::new(193, 194, 132, 0.285, 152.620),
            HBond::new(238, 239, 170, 0.284, 167.530),
            HBond::new(309, 310, 252, 0.265, 170.858),
            HBond::new(358, 359, 301, 0.271, 170.366),
            HBond::new(62, 63, 12, 0.294, 170.336),
            HBond::new(111, 112, 61, 0.299, 164.961),
            HBond::new(118, 119, 61, 0.294, 159.656),
            HBond::new(129, 130, 72, 0.289, 163.824),
            HBond::new(212, 213, 151, 0.287, 164.235),
            HBond::new(231, 232, 170, 0.291, 156.460),
            HBond::new(238, 239, 170, 0.241, 154.949),
            HBond::new(242, 243, 192, 0.283, 150.066),
            HBond::new(291, 292, 241, 0.281, 174.879),
            HBond::new(309, 310, 252, 0.291, 162.452),
            HBond::new(358, 359, 301, 0.260, 153.707),
            HBond::new(62, 63, 12, 0.299, 173.253),
            HBond::new(92, 93, 50, 0.291, 161.575),
            HBond::new(118, 119, 61, 0.298, 159.814),
            HBond::new(129, 130, 72, 0.273, 173.691),
            HBond::new(193, 194, 132, 0.283, 164.727),
            HBond::new(212, 213, 151, 0.283, 154.504),
            HBond::new(238, 239, 170, 0.283, 167.153),
            HBond::new(253, 254, 230, 0.295, 155.103),
            HBond::new(272, 273, 230, 0.283, 174.786),
            HBond::new(309, 310, 252, 0.293, 166.829),
            HBond::new(358, 359, 301, 0.279, 159.971),
        ]
        .into_iter();

        for result in system
            .xtc_iter("test_files/aa_membrane_peptide.xtc")
            .unwrap()
            .hbonds_analyze(
                vec![HBondChain::new(
                    "@protein and elsymbol N O",
                    "@protein and elsymbol N O",
                    "element name hydrogen",
                )],
                vec![(0, 0)],
                0.3,
                150.0,
            )
            .unwrap()
        {
            let (_, hbonds_frame) = result.unwrap();
            assert_eq!(hbonds_frame.len(), 1);

            for bond in hbonds_frame.get(&(0, 0)).unwrap().iter() {
                compare_hbonds(bond, &expected_hbonds.next().unwrap());
            }
        }
    }

    #[test]
    fn test_hbonds_analyze_protein_water() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.tpr").unwrap();

        let hbonds_expected_pp = [
            HBond::new(69, 70, 12, 0.289, 151.553),
            HBond::new(118, 119, 61, 0.278, 162.882),
            HBond::new(129, 130, 72, 0.286, 164.932),
            HBond::new(133, 134, 91, 0.293, 161.708),
            HBond::new(152, 153, 110, 0.280, 173.645),
            HBond::new(193, 194, 132, 0.296, 156.956),
            HBond::new(212, 213, 151, 0.281, 161.991),
            HBond::new(231, 232, 170, 0.288, 150.634),
            HBond::new(238, 239, 170, 0.267, 178.283),
            HBond::new(253, 254, 211, 0.295, 162.803),
            HBond::new(309, 310, 252, 0.287, 174.566),
            HBond::new(332, 333, 290, 0.300, 155.969),
            HBond::new(358, 359, 301, 0.255, 150.333),
        ];

        let hbonds_expected_pw = [
            HBond::new(17725, 17727, 358, 0.287, 161.447),
            HBond::new(19834, 19835, 350, 0.275, 161.685),
            HBond::new(21883, 21885, 361, 0.279, 158.936),
            HBond::new(24292, 24293, 362, 0.274, 165.203),
            HBond::new(26173, 26174, 309, 0.267, 169.959),
            HBond::new(29929, 29930, 331, 0.266, 171.977),
            HBond::new(29929, 29931, 361, 0.267, 152.451),
            HBond::new(30919, 30921, 361, 0.271, 167.135),
            HBond::new(32, 33, 24247, 0.297, 166.586),
            HBond::new(187, 189, 30775, 0.285, 168.139),
        ];

        for result in system
            .gro_iter("test_files/aa_membrane_peptide.gro")
            .unwrap()
            .hbonds_analyze(
                vec![
                    HBondChain::new(
                        "@protein and elsymbol N O",
                        "@protein and elsymbol N O",
                        "element name hydrogen",
                    ),
                    HBondChain::new(
                        "resname SOL and name OW",
                        "resname SOL and name OW",
                        "element name hydrogen",
                    ),
                ],
                vec![(0, 0), (0, 1)],
                0.3,
                150.0,
            )
            .unwrap()
        {
            let (_, hbonds_frame) = result.unwrap();
            assert_eq!(hbonds_frame.len(), 2);

            for (b, bond) in hbonds_frame.get(&(0, 0)).unwrap().iter().enumerate() {
                compare_hbonds(bond, &hbonds_expected_pp[b]);
            }

            for (b, bond) in hbonds_frame.get(&(0, 1)).unwrap().iter().enumerate() {
                compare_hbonds(bond, &hbonds_expected_pw[b]);
            }
        }
    }

    #[test]
    fn test_hbonds_analyze_fail_invalid_pair() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.tpr").unwrap();

        match system
            .gro_iter("test_files/aa_membrane_peptide.gro")
            .unwrap()
            .hbonds_analyze(
                vec![
                    HBondChain::new(
                        "resid 1 to 10 and elsymbol N O",
                        "resid 1 to 10 and elsymbol N O",
                        "resid 1 to 10 and elsymbol H",
                    ),
                    HBondChain::new(
                        "resid 11 to 23 and elsymbol N O",
                        "resid 11 to 23 and elsymbol N O",
                        "resid 11 to 23 and elsymbol H",
                    ),
                ],
                vec![(0, 1), (0, 2)],
                3.0,
                150.0,
            ) {
            Ok(_) => panic!("Function should have failed."),
            Err(HBondError::NonexistentChain(x)) => assert_eq!(x, 2),
            Err(e) => panic!("Unexpected error type `{:?}` returned.", e),
        }
    }

    #[test]
    fn test_hbonds_analyze_fail_unused_chain() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.tpr").unwrap();

        match system
            .gro_iter("test_files/aa_membrane_peptide.gro")
            .unwrap()
            .hbonds_analyze(
                vec![
                    HBondChain::new(
                        "resid 1 to 10 and elsymbol N O",
                        "resid 1 to 10 and elsymbol N O",
                        "resid 1 to 10 and elsymbol H",
                    ),
                    HBondChain::new(
                        "resid 11 to 23 and elsymbol N O",
                        "resid 11 to 23 and elsymbol N O",
                        "resid 11 to 23 and elsymbol H",
                    ),
                ],
                vec![(0, 0)],
                3.0,
                150.0,
            ) {
            Ok(_) => panic!("Function should have failed."),
            Err(HBondError::UnusedChain) => (),
            Err(e) => panic!("Unexpected error type `{:?}` returned.", e),
        }
    }

    #[test]
    fn test_hbonds_analyze_fail_duplicate_pairs_1() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.tpr").unwrap();

        match system
            .gro_iter("test_files/aa_membrane_peptide.gro")
            .unwrap()
            .hbonds_analyze(
                vec![
                    HBondChain::new(
                        "resid 1 to 10 and elsymbol N O",
                        "resid 1 to 10 and elsymbol N O",
                        "resid 1 to 10 and elsymbol H",
                    ),
                    HBondChain::new(
                        "resid 11 to 23 and elsymbol N O",
                        "resid 11 to 23 and elsymbol N O",
                        "resid 11 to 23 and elsymbol H",
                    ),
                ],
                vec![(0, 1), (0, 0), (0, 1)],
                3.0,
                150.0,
            ) {
            Ok(_) => panic!("Function should have failed."),
            Err(HBondError::PairSpecifiedMultipleTimes(x, y)) => {
                assert_eq!(x, 0);
                assert_eq!(y, 1);
            }
            Err(e) => panic!("Unexpected error type `{:?}` returned.", e),
        }
    }

    #[test]
    fn test_hbonds_analyze_fail_duplicate_pairs_2() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.tpr").unwrap();

        match system
            .gro_iter("test_files/aa_membrane_peptide.gro")
            .unwrap()
            .hbonds_analyze(
                vec![
                    HBondChain::new(
                        "resid 1 to 10 and elsymbol N O",
                        "resid 1 to 10 and elsymbol N O",
                        "resid 1 to 10 and elsymbol H",
                    ),
                    HBondChain::new(
                        "resid 11 to 23 and elsymbol N O",
                        "resid 11 to 23 and elsymbol N O",
                        "resid 11 to 23 and elsymbol H",
                    ),
                ],
                vec![(1, 0), (0, 0), (0, 1)],
                3.0,
                150.0,
            ) {
            Ok(_) => panic!("Function should have failed."),
            Err(HBondError::PairSpecifiedMultipleTimes(x, y)) => {
                assert_eq!(x, 0);
                assert_eq!(y, 1);
            }
            Err(e) => panic!("Unexpected error type `{:?}` returned.", e),
        }
    }

    #[test]
    fn test_hbonds_analyze_fail_duplicate_pairs_3() {
        let mut system = System::from_file("test_files/aa_membrane_peptide.tpr").unwrap();

        match system
            .gro_iter("test_files/aa_membrane_peptide.gro")
            .unwrap()
            .hbonds_analyze(
                vec![
                    HBondChain::new(
                        "resid 1 to 10 and elsymbol N O",
                        "resid 1 to 10 and elsymbol N O",
                        "resid 1 to 10 and elsymbol H",
                    ),
                    HBondChain::new(
                        "resid 11 to 23 and elsymbol N O",
                        "resid 11 to 23 and elsymbol N O",
                        "resid 11 to 23 and elsymbol H",
                    ),
                ],
                vec![(0, 0), (1, 0), (0, 0)],
                3.0,
                150.0,
            ) {
            Ok(_) => panic!("Function should have failed."),
            Err(HBondError::PairSpecifiedMultipleTimes(x, y)) => {
                assert_eq!(x, 0);
                assert_eq!(y, 0);
            }
            Err(e) => panic!("Unexpected error type `{:?}` returned.", e),
        }
    }
}
