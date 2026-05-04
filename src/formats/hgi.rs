use anyhow::{Result, anyhow};
use serde::Serialize;
use std::fs::File; 
use std::path::Path;
use std::collections::{HashMap, BTreeSet};
use regex::Regex;

struct HGIRow {
    spectrum: String,
    peptide: String,
    charge: i32,
    total_glycan_composition: Option<String>,
    assigned_modifications: Option<String>,
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
    let spectrum = record.remove("cleaned_scan").unwrap_or_default();
    let peptide =  record.remove("peptide").unwrap_or_default();
    let charge = record.remove("Observed charge state (z)").unwrap_or_default().parse().unwrap_or(0);
    let total_glycan_composition: Option<String> = record.remove("Calculated total glycan mass ");
    let assigned_modifications: Option<String> = record.remove("Identified peptide with modifications (as dalton shifts)");
    let extra_columns = record;

    Ok(HGIRow{
        spectrum,
        peptide,
        charge,
        total_glycan_composition,
        assigned_modifications,
        extra_columns, 
    })
}

pub fn get_filename(raw_filename: &str) -> Result<String> {
    if raw_filename.contains("ModuleA1") {
        Ok(String::from("241202_KEM_HGI_ModuleA1.raw"))
    } else if raw_filename.contains("ModuleA2") {
        Ok(String::from("241202_KEM_HGI_ModuleA2.raw"))
    } else if raw_filename.contains("ModuleA3") {
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

//Modification Name;Modification Mass
//N5H4; 1800.444
pub fn clean_mods(modified_peptide: Option<&str>) -> Vec<ModificationRow> {
    let mut results: Vec<ModificationRow> = Vec::new();

    let re = Regex::new(r"\[(.*?)\]").unwrap();

    for caps in re.captures_iter(modified_peptide.as_deref().unwrap_or("")){
        if let Some(cap) = caps.get(1){
            let resulti = cap.as_str().to_string().replace("+","").split(';').next().unwrap_or("").trim().to_string();

            let mod_row = ModificationRow{
                modification_name: resulti.clone(),
                modification_mass: resulti.clone(),
            };

            results.push(mod_row);
        }
    }

    results
}

//Scan;Sequence;Charge;Modification
//1053;TNATKAAGK;7;N6H4 T3;N2H1 T2
//spectrum, Identified peptide base sequence, charge, Identified peptide with modifications (as dalton shifts)

pub fn get_modification(modified_peptide: Option<&str>) -> Result<String>{
    //println!("{}", modified_peptide);
    let mut results: String = String::new();
    let re = Regex::new(r"(\[[^\]]+\])").unwrap();
    
    for caps in re.captures_iter(modified_peptide.as_deref().unwrap_or("")){
        //let capsi = caps.get(1);

        let mut mod_name: String = String::new();
        let mut amino_acid_number: String = String::new();

        for (i, group) in caps.iter().enumerate().skip(1) {
        if let Some(m) = group {
            //Get amino acids
            let start = m.start();

            let amino_acid = modified_peptide
            .as_deref()
            .unwrap_or("")
            .as_bytes()[start - 1] as char;

            //Get the modification name
            mod_name = m.as_str().to_string().replace("+","").replace("[","").replace("]","").split(';').next().unwrap_or("").trim().to_string();
             
            //Get the amino acid number
            let s = modified_peptide.as_deref().unwrap_or("");
            let before = &s[..start];

            let cleaned: String = before
                .chars()
                .filter(|c| c.is_ascii_uppercase()) // keep only A-Z
                .collect();

            amino_acid_number = cleaned.len().to_string();

            if results.len() == 0{
                results = mod_name + " " + &amino_acid.to_string() + &amino_acid_number
            }else{
                results = results + &mod_name + " " + &amino_acid.to_string() + &amino_acid_number
            }
        }
    }

}
    if results.len() == 0{
        return Ok(results)
    }else{
        return Ok(results)
    }
    
}

pub fn convert_HGI_to_periscope(input_path: &Path, output_dir: &Path) -> Result<()> {
    let file = File::open(input_path)?;
    let mut reader = csv::Reader::from_reader(file);

    let mut identifications = Vec::new();
    let mut modifications = Vec::new();

    for result in reader.deserialize::<HashMap<String, String>>() {
        let row = row_from_record(result?)?;

        //Get the modifications ready
        let mods = clean_mods(row.assigned_modifications.as_deref());
        for m in mods {
            if !modifications.contains(&m){
                modifications.push(m);
            }
        }

        //Get the identifications ready
        let scan_number: String = row.spectrum;
        let sequence = row.extra_columns
            .get("Identified peptide base sequence")
            .cloned()
            .unwrap_or_default();
        let charge = row.charge;
        let modification = match get_modification(row.assigned_modifications.as_deref()) {
            Ok(m) => m,
            Err(e) => {
                eprintln!(
                    "Failed to get modification for spectrum {}: {}",
                    &scan_number,
                    e
                );
                continue; // skip this row or provide a default
            }
        };
        
        let spectrum_file = get_filename(row.extra_columns.get("Data File").map(String::as_str).unwrap_or(""))?;
        let extra_columns = row.extra_columns;

        //Push the modifications
        identifications.push(IdentificationRow {
            scan: scan_number,
            sequence: sequence,
            charge: charge,
            modification: modification,
            spectrum_file: spectrum_file,
            extra_columns: extra_columns,
        });

        let mods = match get_modification(row.assigned_modifications.as_deref()) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("Failed to get modification: {}", e);
                "".to_string() // fallback empty string
            }
        };

        //println!("scan_number is {}, sequence is {}, charge is {}", scan_number, sequence, charge);
    }

    //Now writing
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
