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
    NC,
    DCD,
    TNG,
    LAMMPSTRJ,
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
            Self::NC => "nc",
            Self::DCD => "dcd",
            Self::TNG => "tng",
            Self::LAMMPSTRJ => "lammpstrj",
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
            Some("nc") => FileType::NC,
            Some("dcd") => FileType::DCD,
            Some("tng") => FileType::TNG,
            Some("lammpstrj") => FileType::LAMMPSTRJ,
            Some(_) | None => FileType::Unknown,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::fmt::Write;

    use super::*;

    #[test]
    fn identify() {
        let files = [
            "file.gro",
            "file.pdb",
            "file.xtc",
            "file.ndx",
            "file.trr",
            "file.pqr",
            "file.tpr",
            "file.yaml",
            "file.yml",
            "file.nc",
            "file.dcd",
            "file.tng",
            "file.lammpstrj",
            "file.txt",
            "file",
        ];

        let formats = [
            FileType::GRO,
            FileType::PDB,
            FileType::XTC,
            FileType::NDX,
            FileType::TRR,
            FileType::PQR,
            FileType::TPR,
            FileType::YAML,
            FileType::YAML,
            FileType::NC,
            FileType::DCD,
            FileType::TNG,
            FileType::LAMMPSTRJ,
            FileType::Unknown,
            FileType::Unknown,
        ];

        for (file, format) in files.into_iter().zip(formats.into_iter()) {
            assert_eq!(FileType::from_name(file), format);
        }
    }

    #[test]
    fn display() {
        let files = [
            FileType::Unknown,
            FileType::GRO,
            FileType::PDB,
            FileType::XTC,
            FileType::NDX,
            FileType::TRR,
            FileType::PQR,
            FileType::TPR,
            FileType::YAML,
            FileType::NC,
            FileType::DCD,
            FileType::TNG,
            FileType::LAMMPSTRJ,
        ];

        let mut string = String::new();
        for file in files {
            write!(&mut string, "{} ", file).unwrap();
        }

        assert_eq!(
            string,
            "unknown gro pdb xtc ndx trr pqr tpr yaml nc dcd tng lammpstrj "
        );
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
