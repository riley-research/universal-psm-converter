// Copy the FragPipe implementation from the previous working converter
use anyhow::{anyhow, Result};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fs::File;
use std::path::Path;
use std::collections::{HashMap, BTreeSet};
use regex::Regex;

#[derive(Debug)]
enum ColumnState {
    Missing,
    Empty,
    Value(String),
}

#[derive(Debug)]
struct FragPipeRow {
    spectrum: String,
    peptide: String,
    charge: i32,
    total_glycan_composition: Option<String>,
    assigned_modifications: Option<String>,
    glycan_site_compositions: ColumnState,
    extra_columns: HashMap<String, String>,
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

pub struct GlycanMassCalculator {
    pub library: HashMap<String, f64>,
}

impl GlycanMassCalculator {
    pub fn new() -> Self {
        let library = HashMap::from([
            ("N".to_string(), 203.07937),
            ("H".to_string(), 162.05282),
            ("A".to_string(), 291.09542),
            ("F".to_string(), 146.05791),
            ("G".to_string(), 307.09033),
            ("Pen".to_string(), 132.0423),
            ("Kdn".to_string(), 250.06889),
            ("HexA".to_string(), 176.03209),
            ("p".to_string(), 232.10592),
            ("P".to_string(), 79.96633),
            ("S".to_string(), 79.95682),
            ("C".to_string(), 42.01056),
            ("NH4".to_string(), 17.02655),
            ("Na".to_string(), 21.98194),
            ("Fe".to_string(), 52.91146),
            ("Ca".to_string(), 37.94694),
            ("Al".to_string(), 23.95806),
            ("K".to_string(), 38.96371),
            ("M".to_string(), 27.99491),
            ("U".to_string(), 100.01608),
        ]);

        Self { library }
    }
    pub fn calculate_mass(&self, composition: &str) -> Option<f64> {
        let re = Regex::new(r"([A-Za-z]+)(\d+)").unwrap();
        
        let mut total_mass = 0.0; 
        let mut found_any = false;

        for cap in re.captures_iter(composition) {
            let key = &cap[1];
            let count = cap[2].parse::<f64>().ok()?;

            match self.library.get(key) {
                Some(mass) => {
                    total_mass += mass * count;
                    found_any = true;
                }
                None => return Some(-10.0),
            }
        }

        if found_any { Some(total_mass) } else { None }
    }
}

pub struct GlycanNameConverter {
     pub library: HashMap<String, String>,
}

impl GlycanNameConverter {
    pub fn new() -> Self {
        let library = HashMap::from([
        ("N".to_string(), "HexNAc".to_string()),
        ("H".to_string(), "Hex".to_string()),
        ("A".to_string(), "NeuAc".to_string()),
        ("F".to_string(), "Fuc".to_string()),
        ("G".to_string(), "NeuGc".to_string()),
        ("Pen".to_string(), "Pent".to_string()),
        ("Kdn".to_string(), "KDN".to_string()),
        ("HexA".to_string(), "HexA".to_string()),
        ("p".to_string(), "pseudaminic".to_string()),
        ("P".to_string(), "Phospho".to_string()),
        ("S".to_string(), "Sulfo".to_string()),
        ("C".to_string(), "Acetyl".to_string()),
        ("NH4".to_string(), "NH4".to_string()),
        ("Na".to_string(), "Na".to_string()),
        ("Fe".to_string(), "Fe".to_string()),
        ("Ca".to_string(), "Ca".to_string()),
        ("Al".to_string(), "Al".to_string()),
        ("K".to_string(), "K".to_string()),
        ("M".to_string(), "Formyl".to_string()),
        ("U".to_string(), "Succinyl".to_string()),
    ]);

        Self { library }
    }

    pub fn convert_composition(&self, short_comp: &str) -> String {
            let re = Regex::new(r"([A-Za-z]+)(\d+)").unwrap();
            let mut result = Vec::new();

            for cap in re.captures_iter(short_comp) {
                let symbol = &cap[1]; // e.g., "NH4" or "N"
                let count = &cap[2];  // e.g., "1" or "2"

                // Lookup the full name. If not found, keep the symbol.
                let full_name = self.library.get(symbol)
                    .map(|s| s.as_str())
                    .unwrap_or(symbol);

                result.push(format!("{}({})", full_name, count));
            }

            if result.is_empty() { 
                short_comp.to_string() 
            } else { 
                result.join("") 
            }
        }
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
    let glycan_site_compositions = match record.get("Glycan Site Composition(s)") {
        None => ColumnState::Missing,
        Some(s) if s.is_empty() => ColumnState::Empty,
        Some(s) =>ColumnState::Value(s.clone()),
        };
    
    let extra_columns = record;
    Ok(FragPipeRow {
        spectrum,
        peptide,
        charge,
        total_glycan_composition,
        assigned_modifications,
        glycan_site_compositions,
        extra_columns,
    })
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

fn check_glycan_logic(glycan_mod: Option<&str>, mass: &str) -> String {
    if let Some(glycan) = glycan_mod {
        let parts: Vec<&str> = glycan.split(" % ").collect();
        if parts.len() == 2 {
            let glycan_mass = parts[1].parse::<f64>()
                .map(|m| format!("{:.2}", m))
                .unwrap_or_default();

            if glycan_mass == mass && !glycan_mass.is_empty() {
                return parts[0].to_string();
            }
        }
    }
    mass.to_string()
}

pub fn get_mods(
    glycan_mod: Option<&str>, 
    input_mod: Option<&str>, 
    glycan_site_comp: &ColumnState, 
    calc: &GlycanMassCalculator, 
    name_converter: &GlycanNameConverter,
) -> Option<String> {
    let input_mod = input_mod?;
    if input_mod.is_empty() { return None; }

    let mut formatted_mods = Vec::new();
    let mods: Vec<&str> = input_mod.split(',').collect();

    for modi in mods {
        let mut modi_processed = modi.trim().to_string();

        if modi_processed.contains("N-term") {
            modi_processed = modi_processed.replace("N-term", "1N");
        }

        let localization = modi_processed.chars()
            .take_while(|c| c.is_numeric())
            .collect::<String>();

        let residue_end = modi_processed.find('(').unwrap_or(modi_processed.len());
        let residue = modi_processed[localization.len()..residue_end]
            .chars()
            .filter(|c| c.is_alphabetic())
            .collect::<String>();

        let mass = if let (Some(start), Some(end)) = (modi_processed.find('('), modi_processed.find(')')) {
            modi_processed[start + 1..end].parse::<f64>()
                .map(|m| format!("{:.2}", m))
                .unwrap_or_default()
        } else {
            String::new()
        };

        let modification = match glycan_site_comp {
            ColumnState::Missing | ColumnState::Empty => {
                check_glycan_logic(glycan_mod, &mass) 
            }
            ColumnState::Value(s) => {
                let sites: Vec<&str> = s.split(',').collect();
                let mut found_match = None;

                for site in sites {
                    let site_trimmed = site.trim();
                    if let Some(calculated_mass) = calc.calculate_mass(site_trimmed) {
                        let formatted_calc_mass = format!("{:.2}", calculated_mass);
                        //println!("Formatted_calc_mass {:?}, site {:?}, formatted {:?}", formatted_calc_mass, site_trimmed, name_converter.convert_composition(site_trimmed));
                        if formatted_calc_mass == mass {
                            found_match = Some(name_converter.convert_composition(site_trimmed));
                            break;
                        }
                    }
                }
                found_match.unwrap_or_else(|| check_glycan_logic(glycan_mod, &mass))
            }
        };

        let mod_line = format!("{} {}:{}", modification, residue, localization);
        formatted_mods.push(mod_line);
    }

    if formatted_mods.is_empty() { None } else { Some(formatted_mods.join(";")) }
}

/// Format modifications for modifications.csv
/// Equivalent to cleanMods in R script
pub fn clean_mods(
    glycan_mod: Option<&str>, 
    input_mod: Option<&str>,
    glycan_site_comp: &ColumnState, 
    calc: &GlycanMassCalculator, 
    name_converter: &GlycanNameConverter,
) -> Vec<ModificationRow> {
    let mut results = Vec::new();

    let input_mod_str = match input_mod {
        Some(m) if !m.is_empty() => m,
        _ => return results,
    };

    let mods: Vec<&str> = input_mod_str.split(',').collect();

    for modi in mods {
        let mut modi_processed = modi.trim().to_string();

        // Handle N-term
        if modi_processed.contains("N-term") {
            modi_processed = modi_processed.replace("N-term", "1N");
        }

        // Extract residue (amino acid)
        let residue_end = modi_processed.find('(').unwrap_or(modi_processed.len());
        let residue = modi_processed[..residue_end]
            .chars()
            .filter(|c| c.is_alphabetic())
            .collect::<String>();

        // Extract and format mass
        let mass_raw = if let (Some(start), Some(end)) = (modi_processed.find('('), modi_processed.find(')')) {
            modi_processed[start + 1..end].to_string()
        } else {
            String::new()
        };

        let mass_rounded = mass_raw.parse::<f64>()
            .map(|m| format!("{:.2}", m))
            .unwrap_or_else(|_| mass_raw.clone());

        // Determine modification name using site-specific logic
        let modification_name = match glycan_site_comp {
            ColumnState::Missing | ColumnState::Empty => {
                check_glycan_logic(glycan_mod, &mass_rounded) 
            }
            ColumnState::Value(s) => {
                let sites: Vec<&str> = s.split(',').collect();
                let mut found_match = None;

                for site in sites {
                    let site_trimmed = site.trim();
                    if let Some(calculated_mass) = calc.calculate_mass(site_trimmed) {
                        let formatted_calc_mass = format!("{:.2}", calculated_mass);
                        if formatted_calc_mass == mass_rounded {
                            found_match = Some(name_converter.convert_composition(site_trimmed));
                            break;
                        }
                    }
                }
                found_match.unwrap_or_else(|| check_glycan_logic(glycan_mod, &mass_rounded))
            }
        };

        results.push(ModificationRow {
            modification_name: format!("{} {}", modification_name, residue),
            modification_mass: mass_raw,
        });
    }

    results
}

pub fn convert_fragpipe_to_periscope(input_path: &Path, output_dir: &Path) -> Result<()> {
    // Read the TSV file
    let file = File::open(input_path)?;
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(file);

    let mut identifications = Vec::new();
    let mut modifications_set = HashSet::new();
    let calc = GlycanMassCalculator::new();
    let name_converter = GlycanNameConverter::new();

    //for result in reader.deserialize::<FragPipeRow>() {
    //    let row = result?;
    for result in reader.deserialize::<HashMap<String, String>>() {
        let row = row_from_record(result?)?;
        //println!("{:#?}", row.glycan_site_compositions);

        let scan = get_spec_nr(&row.spectrum)?;

        let spectral_file = row.spectrum.split('.').next()
            .unwrap_or("")
            .to_string();

        let modification = get_mods(row.total_glycan_composition.as_deref(), row.assigned_modifications.as_deref(),
                                    &row.glycan_site_compositions, &calc, &name_converter)
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
            let clean_mods = clean_mods(row.total_glycan_composition.as_deref(), row.assigned_modifications.as_deref(),
                                        &row.glycan_site_compositions, &calc, &name_converter);
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

#[cfg(test)]
mod lib_tests {
    use super::*;
    use std::path::Path;
    use crate::formats::fragpipe::convert_fragpipe_to_periscope;
    
    #[test]
    fn test_convert_fragpipe_to_periscope() {
        let input = Path::new("Z:\\Tim\\Periscope_test_files\\psm.tsv");
        let output = Path::new("Z:\\Tim\\Periscope_test_files\\");

        let result = convert_fragpipe_to_periscope(input, output);

        assert!(result.is_ok());
    }
}