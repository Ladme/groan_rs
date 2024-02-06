// Released under MIT License.
// Copyright (c) 2023-2024 Ladislav Bartos

//! Enum capturing file types supported by `groan_rs`.

use std::path::Path;

/// Types of files supported by `groan_rs`.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum FileType {
    Unknown,
    GRO,
    PDB,
    XTC,
    NDX,
    TRR,
    PQR,
}

impl FileType {
    /// Identify file type from the name of the file (based on file extension).
    pub fn from_name(filename: impl AsRef<Path>) -> FileType {
        let extension = match filename.as_ref().extension() {
            Some(x) => x,
            None => return FileType::Unknown,
        };

        match extension.to_str() {
            Some("gro") => FileType::GRO,
            Some("pdb") => FileType::PDB,
            Some("xtc") => FileType::XTC,
            Some("ndx") => FileType::NDX,
            Some("trr") => FileType::TRR,
            Some("pqr") => FileType::PQR,
            Some(_) | None => FileType::Unknown,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identify_gro() {
        assert_eq!(FileType::from_name("file.gro"), FileType::GRO);
    }

    #[test]
    fn identify_pdb() {
        assert_eq!(FileType::from_name("file.pdb"), FileType::PDB);
    }

    #[test]
    fn identify_xtc() {
        assert_eq!(FileType::from_name("file.xtc"), FileType::XTC);
    }

    #[test]
    fn identify_ndx() {
        assert_eq!(FileType::from_name("file.ndx"), FileType::NDX);
    }

    #[test]
    fn identify_trr() {
        assert_eq!(FileType::from_name("file.trr"), FileType::TRR);
    }

    #[test]
    fn identity_pqr() {
        assert_eq!(FileType::from_name("file.pqr"), FileType::PQR);
    }

    #[test]
    fn identify_unknown() {
        assert_eq!(FileType::from_name("file.txt"), FileType::Unknown);
    }

    #[test]
    fn identify_noextension() {
        assert_eq!(FileType::from_name("file"), FileType::Unknown);
    }
}
