// Copy the FragPipe implementation from the previous working converter
use anyhow::{anyhow, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fs::File;
use std::path::Path;
use std::collections::{HashMap, BTreeSet};

#[derive(Debug)]
struct FragPipeRow {
    spectrum: String,
    peptide: String,
    charge: i32,
    total_glycan_composition: Option<String>,
    assigned_modifications: Option<String>,
    extra_columns: HashMap<String, String>,
}

/// Build a FragPipeRow from a CSV record map. Known columns are extracted; the rest go into extra_columns.
fn row_from_record(mut record: HashMap<String, String>) -> Result<FragPipeRow> {
    let spectrum = record.remove("Spectrum").unwrap_or_default();
    let peptide = record.remove("Peptide").unwrap_or_default();
    let charge = record
        .remove("Charge")
        .unwrap_or_default()
        .parse::<i32>()
        .unwrap_or(0);
    let total_glycan_composition = record.remove("Total Glycan Composition").filter(|s| !s.is_empty());
    let assigned_modifications = record.remove("Assigned Modifications").filter(|s| !s.is_empty());
    // Whatever is left is "extra" columns
    let extra_columns = record;
    Ok(FragPipeRow {
        spectrum,
        peptide,
        charge,
        total_glycan_composition,
        assigned_modifications,
        extra_columns,
    })
}

#[derive(Debug, Serialize)]
struct IdentificationRow {
    #[serde(rename = "Scan")]
    scan: i32,
    #[serde(rename = "Sequence")]
    sequence: String,
    #[serde(rename = "Charge")]
    charge: i32,
    #[serde(rename = "Modification")]
    modification: String,
    #[serde(rename = "spectral_file")]
    spectrum_file: String,
    #[serde(flatten)]
    extra_columns: HashMap<String, String>,
}

#[derive(Debug, Serialize, Clone, PartialEq, Eq, Hash)]
pub struct ModificationRow {
    #[serde(rename = "Modification Name")]
    modification_name: String,
    #[serde(rename = "Modification Mass")]
    modification_mass: String,
}

/// Extract scan number from spectrum field
/// Equivalent to getSpecNr in R script
pub fn get_spec_nr(spectrum: &str) -> Result<i32> {
    let parts: Vec<&str> = spectrum.split('.').collect();
    if parts.len() < 2 {
        return Err(anyhow!("Invalid spectrum format"));
    }

    // Get the second-to-last element
    let scan_str = parts[parts.len() - 2];

    // Remove leading zeros and convert to number
    let scan_number = scan_str.trim_start_matches('0')
        .parse::<i32>()
        .or_else(|_| scan_str.parse::<i32>())
        .unwrap_or(0);

    Ok(scan_number)
}

/// Format modifications for identifications.csv
/// Equivalent to getMods in R script
pub fn get_mods(glycan_mod: Option<&str>, input_mod: Option<&str>) -> Option<String> {
    let input_mod = input_mod?;
    if input_mod.is_empty() {
        return None;
    }

    let mut formatted_mods = Vec::new();
    let mods: Vec<&str> = input_mod.split(',').collect();

    for modi in mods {
        let mut modi_processed = modi.trim().to_string();

        // Handle N-term
        if modi_processed.contains("N-term") {
            modi_processed = modi_processed.replace("N-term", "1N");
        }

        // Extract residue number (localization)
        let localization = modi_processed.chars()
            .take_while(|c| c.is_numeric())
            .collect::<String>();

        // Extract amino acid (residue)
        let residue_end = modi_processed.find('(').unwrap_or(modi_processed.len());
        let residue_str = &modi_processed[localization.len()..residue_end];
        let residue = residue_str.chars()
            .filter(|c| c.is_alphabetic())
            .collect::<String>();

        // Extract mass from parentheses
        let mass = if let (Some(start), Some(end)) = (modi_processed.find('('), modi_processed.find(')')) {
            modi_processed[start + 1..end].parse::<f64>()
                .map(|m| format!("{:.2}", m))
                .unwrap_or_default()
        } else {
            String::new()
        };

        // Check if it's a glycan modification
        let modification = if let Some(glycan) = glycan_mod {
            if !glycan.is_empty() {
                let parts: Vec<&str> = glycan.split(" % ").collect();
                if parts.len() == 2 {
                    let glycan_name = parts[0];
                    let glycan_mass = parts[1].parse::<f64>()
                        .map(|m| format!("{:.2}", m))
                        .unwrap_or_default();

                    if glycan_mass == mass && !glycan_mass.is_empty() {
                        glycan_name.to_string()
                    } else {
                        mass.clone()
                    }
                } else {
                    mass.clone()
                }
            } else {
                mass.clone()
            }
        } else {
            mass.clone()
        };

        // Format as "mod residue:localization"
        let mod_line = format!("{} {}:{}", modification, residue, localization);
        formatted_mods.push(mod_line);
    }

    if formatted_mods.is_empty() {
        None
    } else {
        Some(formatted_mods.join(";"))
    }
}

/// Format modifications for modifications.csv
/// Equivalent to cleanMods in R script
pub fn clean_mods(glycan_mod: Option<&str>, input_mod: Option<&str>) -> Vec<ModificationRow> {
    let mut results = Vec::new();

    let input_mod = match input_mod {
        Some(m) if !m.is_empty() => m,
        _ => return results,
    };

    let mods: Vec<&str> = input_mod.split(',').collect();

    for modi in mods {
        let mut modi_processed = modi.trim().to_string();

        // Handle N-term
        if modi_processed.contains("N-term") {
            modi_processed = modi_processed.replace("N-term", "1N");
        }

        // Extract amino acid
        let residue_end = modi_processed.find('(').unwrap_or(modi_processed.len());
        let residue_str = &modi_processed[..residue_end];
        let residue = residue_str.chars()
            .filter(|c| c.is_alphabetic())
            .collect::<String>();

        // Extract mass from parentheses
        let mass = if let (Some(start), Some(end)) = (modi_processed.find('('), modi_processed.find(')')) {
            modi_processed[start + 1..end].to_string()
        } else {
            String::new()
        };

        let mass_rounded = mass.parse::<f64>()
            .map(|m| format!("{:.2}", m))
            .unwrap_or(mass.clone());

        // Determine modification name
        let modification_name = if let Some(glycan) = glycan_mod {
            if !glycan.is_empty() {
                let parts: Vec<&str> = glycan.split(" % ").collect();
                if parts.len() == 2 {
                    let glycan_name = parts[0];
                    let glycan_mass = parts[1].parse::<f64>()
                        .map(|m| format!("{:.2}", m))
                        .unwrap_or_default();

                    if glycan_mass == mass_rounded && !glycan_mass.is_empty() {
                        format!("{} {}", glycan_name, residue)
                    } else {
                        format!("{} {}", mass_rounded, residue)
                    }
                } else {
                    format!("{} {}", mass_rounded, residue)
                }
            } else {
                format!("{} {}", mass_rounded, residue)
            }
        } else {
            format!("{} {}", mass_rounded, residue)
        };

        results.push(ModificationRow {
            modification_name,
            modification_mass: mass,
        });
    }

    results
}

pub fn convert_fragpipe_to_ipsa(input_path: &Path, output_dir: &Path) -> Result<()> {
    // Read the TSV file
    let file = File::open(input_path)?;
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file);

    let mut identifications = Vec::new();
    let mut modifications_set = HashSet::new();

    //for result in reader.deserialize::<FragPipeRow>() {
    //    let row = result?;
    for result in reader.deserialize::<HashMap<String, String>>() {
        let row = row_from_record(result?)?;
        //println!("{:#?}", row.extra_columns);

        // Extract scan number
        let scan = get_spec_nr(&row.spectrum)?;

        // Extract spectral file (first part before first dot)
        let spectral_file = row.spectrum.split('.').next()
            .unwrap_or("")
            .to_string();

        // Get modifications for identifications
        let modification = get_mods(row.total_glycan_composition.as_deref(), row.assigned_modifications.as_deref())
            .unwrap_or_default();

        // Create identification row
        identifications.push(IdentificationRow {
            scan,
            sequence: row.peptide.clone(),
            charge: row.charge,
            modification,
            spectrum_file: spectral_file,
            extra_columns: row.extra_columns,
        });

        // Collect modifications for modifications.csv
        if row.total_glycan_composition.is_some() && row.assigned_modifications.is_some() {
            let clean_mods = clean_mods(row.total_glycan_composition.as_deref(), row.assigned_modifications.as_deref());
            for mod_row in clean_mods {
                modifications_set.insert(mod_row);
            }
        }
    }

    // Write Identifications.csv
    let identifications_path = output_dir.join("Identifications.csv");
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
    
    // Write Modifications.csv
    let modifications_path = output_dir.join("Modifications.csv");
    let mut wtr = csv::Writer::from_path(&modifications_path)?;

    // Sort modifications by name for consistent output
    let mut modifications_vec: Vec<_> = modifications_set.into_iter().collect();
    modifications_vec.sort_by(|a, b| a.modification_name.cmp(&b.modification_name));

    for row in modifications_vec {
        wtr.serialize(row)?;
    }
    wtr.flush()?;

    Ok(())
}