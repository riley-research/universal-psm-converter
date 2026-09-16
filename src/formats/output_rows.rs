use serde::Serialize;
use std::collections::HashMap;

#[derive(Debug, Serialize)]
pub struct IdentificationRow {
    #[serde(rename = "Scan")]
    pub scan: String,
    #[serde(rename = "Sequence")]
    pub sequence: String,
    #[serde(rename = "Charge")]
    pub charge: i32,
    #[serde(rename = "Modification")]
    pub modification: String,
    #[serde(rename = "spectral_file")]
    pub spectral_file: String,
    #[serde(flatten)]
    pub extra_columns: Option<HashMap<String, String>>,
}

#[derive(Debug, Serialize, Clone, PartialEq, Eq, Hash)]
pub struct ModificationRow {
    #[serde(rename = "Modification Name")]
    pub modification_name: String,
    #[serde(rename = "Modification Mass")]
    pub modification_mass: String,
}