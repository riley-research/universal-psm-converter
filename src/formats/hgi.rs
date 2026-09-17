use anyhow::{Result, anyhow};
use serde::Serialize;
use std::fs::File; 
use std::path::Path;
use std::collections::{HashMap, BTreeSet};
use regex::Regex;

#[derive(Debug)]
struct HGIRow {
    spectrum: String,
    peptide: String,
    charge: i32,
    modification_for_identification: Option<String>,
    modification_for_modification: Option<String>,
    file: String,
    extra_columns: HashMap<String, String>,
}

#[derive(Debug, Serialize, Clone, PartialEq, Eq, Hash)]
pub struct ModificationRow {
    #[serde(rename = "Modification Name")]
    modification_name: String,
    #[serde(rename = "Modification Mass")]
    modification_mass: String,
}

#[derive(Debug, Serialize)]
pub struct IdentificationRow {
    #[serde(rename = "Scan")]
    pub scan: String,
    #[serde(rename = "Sequence")]
    pub sequence: String,
    #[serde(rename = "Charge")]
    pub charge: i32,
    #[serde(rename = "Modification")]
    pub modification : String,
    #[serde(rename = "spectral_file")]
    spectrum_file: String,
    #[serde(flatten)]
    extra_columns: HashMap<String, String>,
}

fn row_from_record(mut record: HashMap<String, String>) -> Result<HGIRow> {
    let spectrum = record.remove("Scan").unwrap();
    let peptide =  record.remove("Sequence").unwrap();
    let charge = record.remove("Charge").unwrap().parse().unwrap_or(0);
    let modification_for_identification: Option<String> = record.remove("modification_for_identification");
    let modification_for_modification: Option<String> = record.remove("modification_for_modification");
    let file = record.remove("File").unwrap();
    let extra_columns = record;

    Ok(HGIRow{
        spectrum,
        peptide,
        charge,
        modification_for_identification,
        modification_for_modification,
        file,
        extra_columns,
    })
}

pub fn get_filename(raw_filename: &str) -> Result<String> {
    if raw_filename.contains("ModuleA_GlycoPSMs_Sample_1") {
        Ok(String::from("241202_KEM_HGI_ModuleA1.raw"))
    } else if raw_filename.contains("ModuleA_GlycoPSMs_Sample_2") {
        Ok(String::from("241202_KEM_HGI_ModuleA2.raw"))
    } else if raw_filename.contains("ModuleA_GlycoPSMs_Sample_3") {
        Ok(String::from("241202_KEM_HGI_ModuleA3.raw"))
    
    } else if raw_filename.contains("ModuleB1") {
        Ok(String::from("241231_KEM_HGI_ModuleB1_EThcD.raw"))
    } else if raw_filename.contains("ModuleB2") {
        Ok(String::from("241231_KEM_HGI_ModuleB2_EThcD.raw"))
    } else if raw_filename.contains("ModuleB3") {
        Ok(String::from("241231_KEM_HGI_ModuleB3_EThcD.raw"))
    
    } else if raw_filename.contains("ModuleC1_EThcD") {
        Ok(String::from("250324_VC_HGI_ModuleC_1_EThcD.raw"))
    } else if raw_filename.contains("ModuleC1_HCD") {
        Ok(String::from("250324_VC_HGI_ModuleC_1_HCD.raw"))
    
    } else if raw_filename.contains("ModuleC2_EThcD") {
        Ok(String::from("250324_VC_HGI_ModuleC_2_EThcD.raw"))
    } else if raw_filename.contains("ModuleC2_HCD") {
        Ok(String::from("250324_VC_HGI_ModuleC_2_HCD.raw"))
   
    }else if raw_filename.contains("ModuleC3_EThcD") {
        Ok(String::from("250324_VC_HGI_ModuleC_3_EThcD.raw"))
    } else if raw_filename.contains("ModuleC4_HCD") {
        Ok(String::from("250324_VC_HGI_ModuleC_3_HCD.raw"))
    
    }else {
        Err(anyhow!("No matching module found in filename: {}", raw_filename))
    }
}

fn clean_mods(value: Option<&str>) -> Vec<ModificationRow> {
    value
        .unwrap_or("")
        .split(';')
        .filter_map(|modification| {
            let modification = modification.trim();

            if modification.is_empty() {
                return None;
            }

            let (name, mass) = modification.rsplit_once(char::is_whitespace)?;

            Some(ModificationRow {
                modification_name: name.trim().to_string(),
                modification_mass: mass.trim().parse().ok()?,
            })
        })
        .collect()
}

pub fn convert_HGI_to_periscope(input_path: &Path, output_dir: &Path) -> Result<()> {
    let file = File::open(input_path)?;
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file);

    let mut identifications = Vec::new();
    let mut modifications = Vec::new();

    for result in reader.deserialize::<HashMap<String, String>>() {
        let row = row_from_record(result?)?;

        //Get the modifications ready
        let mods = clean_mods(row.modification_for_modification.as_deref());
        for m in mods {
            if !modifications.contains(&m){
                modifications.push(m);
            }
        }

        //Get the identifications ready
        let scan_number: String = row.spectrum;
        let sequence = row.peptide;
        let charge = row.charge;
        let modification = row.modification_for_identification;
        let spectrum_file = get_filename(&row.file).unwrap();
        let extra_columns = row.extra_columns;

        identifications.push(IdentificationRow {
            scan: scan_number,
            sequence: sequence,
            charge: charge,
            modification: modification.unwrap(),
            spectrum_file: spectrum_file,
            extra_columns: extra_columns,
        });

    }

        write_identifications(output_dir, identifications)?;
        write_modification(output_dir, modifications)?;

        Ok(())
}

fn write_modification(od: &Path, mods: Vec<ModificationRow>) -> Result<()> {

    let modifications_path = od.join("Modifications.csv");

    let mut wtr = csv::Writer::from_path(&modifications_path)?;
    for row in mods {
        wtr.serialize(row)?;
    }
    wtr.flush()?;
    Ok(())
}

fn write_identifications(od: &Path, identifications: Vec<IdentificationRow>) -> Result<()> {
    let identifications_path = od.join("Identifications.csv");
    let mut wtr = csv::Writer::from_path(&identifications_path)?;

    // Collect all extra column names (sorted) so header and row order match
    let extra_keys: Vec<String> = identifications
        .iter()
        .flat_map(|r| r.extra_columns.keys().cloned())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();

    let header: Vec<String> = [
        "Scan",
        "Sequence",
        "Charge",
        "Modification",
        "spectral_file",
    ]
    .into_iter()
    .map(String::from)
    .chain(extra_keys.clone())
    .collect();
    wtr.write_record(&header)?;

    for row in &identifications {
        let mut record = vec![
            row.scan.to_string(),
            row.sequence.clone(),
            row.charge.to_string(),
            row.modification.clone(),
            row.spectrum_file.clone(),
        ];
        for k in &extra_keys {
            record.push(row.extra_columns.get(k).cloned().unwrap_or_default());
        }
        wtr.write_record(&record)?;
    }
    wtr.flush()?;
    
    Ok(())
}

/*
#[cfg(test)]
mod lib_tests {
    use super::*;
    use std::path::Path;
    use crate::formats::hgi::convert_HGI_to_periscope;
    
   #[test]
fn test_convert_hgi_to_periscope() {
    let input = Path::new(
        r"C:\Users\tim_v\Documents\PostDoc\HGI\ModuleA_GlycoPSMs_Sample_1.csv"
    );

    let output = Path::new(
        r"C:\Users\tim_v\Documents\PostDoc\HGI\"
    );

    let result = convert_HGI_to_periscope(input, output);

    match result {
        Ok(_) => println!("Success"),
        Err(e) => panic!("Conversion failed: {}", e),
    }
}
} 

//cargo test formats::hgi::lib_tests::test_convert_hgi_to_periscope -- --nocapture
*/