use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SoftwareFormat {
    FragPipe,
    Byonic,
    GlycoDecipher,
    GPQuest,
    OPair,
    PGlyco,
}

impl std::fmt::Display for SoftwareFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SoftwareFormat::FragPipe => write!(f, "FragPipe"),
            SoftwareFormat::Byonic => write!(f, "Byonic"),
            SoftwareFormat::GlycoDecipher => write!(f, "Glyco-Decipher"),
            SoftwareFormat::GPQuest => write!(f, "GPQuest"),
            SoftwareFormat::OPair => write!(f, "OPair"),
            SoftwareFormat::PGlyco => write!(f, "pGlyco"),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConversionStats {
    pub format: String,
    pub identifications_written: usize,
    pub modifications_written: usize,
    pub processing_time_ms: u64,
    pub peak_memory_usage: Option<usize>,
}

#[derive(Debug, Clone, Serialize)]
pub struct IdentificationRecord {
    pub scan: String,
    pub sequence: String,
    pub charge: i32,
    pub modification: String,
    pub spectral_file: String,
}

#[derive(Debug, Clone, Serialize, PartialEq, Eq, Hash)]
pub struct ModificationRecord {
    pub modification_name: String,
    pub modification_mass: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct PGlycoIdentificationRecord {
    pub scan: String,
    pub peptide: String,
    pub charge: i32,
    pub modification: String,
    pub spectral_file: String,
}

impl IdentificationRecord {
    pub fn new(scan: String, sequence: String, charge: i32, modification: String, spectral_file: String) -> Self {
        Self { scan, sequence, charge, modification, spectral_file }
    }

    pub fn to_csv_row(&self) -> Vec<String> {
        vec![
            self.scan.clone(),
            self.sequence.clone(),
            self.charge.to_string(),
            self.modification.clone(),
            self.spectral_file.clone(),
        ]
    }
}

impl PGlycoIdentificationRecord {
    pub fn new(scan: String, peptide: String, charge: i32, modification: String, spectral_file: String) -> Self {
        Self { scan, peptide, charge, modification, spectral_file }
    }

    pub fn to_csv_row(&self) -> Vec<String> {
        vec![
            self.scan.clone(),
            self.peptide.clone(),
            self.charge.to_string(),
            self.modification.clone(),
            self.spectral_file.clone(),
        ]
    }
}

impl ModificationRecord {
    pub fn new(modification_name: String, modification_mass: String) -> Self {
        Self { modification_name, modification_mass }
    }

    pub fn to_csv_row(&self) -> Vec<String> {
        vec![self.modification_name.clone(), self.modification_mass.clone()]
    }
}