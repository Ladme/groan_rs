// Released under MIT License.
// Copyright (c) 2023-2025 Ladislav Bartos

//! Enum capturing file types supported by `groan_rs`.

use std::fmt::Display;
use std::path::Path;

/// Types of files supported by `groan_rs`.
#[derive(Debug, PartialEq, Clone, Copy, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "serde",
    serde(rename_all = "lowercase", deny_unknown_fields)
)]
pub enum FileType {
    Unknown,
    GRO,
    PDB,
    XTC,
    NDX,
    TRR,
    PQR,
    TPR,
    YAML,
}

impl Display for FileType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let string = match self {
            Self::Unknown => "unknown",
            Self::GRO => "gro",
            Self::PDB => "pdb",
            Self::XTC => "xtc",
            Self::NDX => "ndx",
            Self::TRR => "trr",
            Self::PQR => "pqr",
            Self::TPR => "tpr",
            Self::YAML => "yaml",
        };

        write!(f, "{}", string)
    }
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
            Some("tpr") => FileType::TPR,
            Some("yml") | Some("yaml") => FileType::YAML,
            Some(_) | None => FileType::Unknown,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::fmt::Write;

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
    fn identify_tpr() {
        assert_eq!(FileType::from_name("file.tpr"), FileType::TPR);
    }

    #[test]
    fn identity_yaml() {
        assert_eq!(FileType::from_name("file.yaml"), FileType::YAML);
        assert_eq!(FileType::from_name("file.yml"), FileType::YAML);
    }

    #[test]
    fn identify_unknown() {
        assert_eq!(FileType::from_name("file.txt"), FileType::Unknown);
    }

    #[test]
    fn identify_noextension() {
        assert_eq!(FileType::from_name("file"), FileType::Unknown);
    }

    #[test]
    fn display() {
        let files = vec![
            FileType::Unknown,
            FileType::GRO,
            FileType::PDB,
            FileType::XTC,
            FileType::NDX,
            FileType::TRR,
            FileType::PQR,
            FileType::TPR,
            FileType::YAML,
        ];

        let mut string = String::new();
        for file in files {
            write!(&mut string, "{} ", file).unwrap();
        }

        assert_eq!(string, "unknown gro pdb xtc ndx trr pqr tpr yaml ");
    }
}

#[cfg(test)]
#[cfg(feature = "serde")]
mod serde_tests {
    use super::*;

    #[test]
    fn filetype_to_yaml() {
        let files = vec![
            FileType::Unknown,
            FileType::GRO,
            FileType::PDB,
            FileType::XTC,
            FileType::NDX,
            FileType::TRR,
            FileType::PQR,
            FileType::TPR,
            FileType::YAML,
        ];

        let string = serde_yaml::to_string(&files).unwrap();
        let expected = "- unknown
- gro
- pdb
- xtc
- ndx
- trr
- pqr
- tpr
- yaml
";
        assert_eq!(string, expected);
    }

    #[test]
    fn filetype_from_yaml() {
        let string = "- unknown
- gro
- pdb
- xtc
- ndx
- trr
- pqr
- tpr
- yaml
";

        let filetypes: Vec<FileType> = serde_yaml::from_str(string).unwrap();

        let expected = vec![
            FileType::Unknown,
            FileType::GRO,
            FileType::PDB,
            FileType::XTC,
            FileType::NDX,
            FileType::TRR,
            FileType::PQR,
            FileType::TPR,
            FileType::YAML,
        ];

        assert_eq!(filetypes, expected);
    }
}
