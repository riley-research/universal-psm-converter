use anyhow::Result;
use regex::Regex;
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::Path;

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
}

#[derive(Debug, Serialize, Clone, PartialEq, Eq, Hash)]
pub struct ModificationRow {
    #[serde(rename = "Modification Name")]
    pub modification_name: String,
    #[serde(rename = "Modification Mass")]
    pub modification_mass: String,
}

pub struct GlycoDecipherConverter;

impl GlycoDecipherConverter {
    pub fn convert_to_ipsa(input_path: &Path, output_dir: &Path) -> Result<()> {
        let file = File::open(input_path)?;
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(file);

        let mut identifications = Vec::new();
        let mut modifications_set = HashSet::new();

        let headers = reader.headers()?.clone();
        let is_database_glycan = headers.iter().any(|h| h == "GlycanComposition");

        for record in reader.deserialize::<HashMap<String, String>>() {
            let row = record?;

            let scan = row.get("Scan").unwrap_or(&String::new()).clone();
            let peptide = row.get("Peptide").unwrap_or(&String::new()).clone();
            let charge = row.get("Charge").map_or("0", |v| v).parse::<i32>()?;
            let modification = row.get("Modification").unwrap_or(&String::new()).clone();
            let file_name = row.get("File").unwrap_or(&String::new()).clone();

            let variable_mods = Self::get_variable_mods(&modification, &peptide);
            let glycosite = Self::find_n_glyco_residue(&peptide);

            let (final_modification, glycan_mass, glycan_composition) = if is_database_glycan {
                let glycan_composition = row.get("GlycanComposition").unwrap_or(&String::new()).clone();
                let glycan_mass = row.get("GlycanMass").map_or("0", |v| v).parse::<f64>().unwrap_or(0.0);
                let final_mod = Self::add_glycan_to_modifications(&variable_mods, &glycan_composition, glycosite);
                (final_mod, glycan_mass, glycan_composition)
            } else {
                let glycan_exp_mass = row.get("GlycanExpMass").map_or("0", |v| v).parse::<f64>().unwrap_or(0.0);
                let final_mod = Self::add_glycan_to_modifications_psm(&variable_mods, glycan_exp_mass, glycosite);
                (final_mod, glycan_exp_mass, glycan_exp_mass.to_string())
            };

            identifications.push(IdentificationRow {
                scan,
                sequence: peptide,
                charge,
                modification: final_modification.clone(),
                spectral_file: file_name,
            });

            if !final_modification.is_empty() {
                Self::process_modifications_for_csv(&final_modification, glycan_mass, &glycan_composition, &mut modifications_set);
            }
        }

        let mut modifications: Vec<_> = modifications_set.into_iter().collect();
        modifications.sort_by(|a, b| a.modification_name.cmp(&b.modification_name));

        Self::write_outputs(identifications, modifications, output_dir)
    }

    fn get_variable_mods(mods: &str, sequence: &str) -> String {
        if mods.is_empty() {
            return String::new();
        }

        let mut formatted_mods = Vec::new();

        for mod_str in mods.split(';') {
            let mod_str = mod_str.trim();
            if mod_str.is_empty() {
                continue;
            }

            let parts: Vec<&str> = mod_str.split(',').collect();
            if parts.len() < 2 {
                continue;
            }

            let residue_number = parts[0].parse::<usize>().unwrap_or(0);

            let mod_name_part = parts[1];

            let mass_re = Regex::new(r"\(\+?([0-9.]+)\)").unwrap();
            let mass = if let Some(cap) = mass_re.captures(mod_name_part) {
                cap[1].to_string()
            } else {
                continue;
            };

            let modified_residue = if residue_number > 0 && residue_number <= sequence.len() {
                sequence.chars().nth(residue_number - 1).unwrap_or('X')
            } else {
                'X'
            };

            formatted_mods.push(format!("{} {}:{}", mass, modified_residue, residue_number));
        }

        formatted_mods.join(";")
    }

    fn find_n_glyco_residue(sequence: &str) -> i32 {
        let n_sequon_re = Regex::new(r"N[^P][TS]").unwrap();
        if let Some(mat) = n_sequon_re.find(sequence) {
            return (mat.start() + 1) as i32;
        }

        if let Some(pos) = sequence.rfind('N') {
            return (pos + 1) as i32;
        }

        -1
    }

    fn add_glycan_to_modifications(mods: &str, glycan_composition: &str, glycosite: i32) -> String {
        let glycan_mod = format!("{} N:{}", glycan_composition, glycosite);

        if mods.is_empty() {
            glycan_mod
        } else {
            format!("{};{}", mods, glycan_mod)
        }
    }

    fn add_glycan_to_modifications_psm(mods: &str, glycan_mass: f64, glycosite: i32) -> String {
        let glycan_mod = format!("{} N:{}", glycan_mass, glycosite);

        if mods.is_empty() {
            glycan_mod
        } else {
            format!("{};{}", mods, glycan_mod)
        }
    }

    fn process_modifications_for_csv(
        modification: &str,
        glycan_mass: f64,
        glycan_composition: &str,
        modifications_set: &mut HashSet<ModificationRow>
    ) {
        for mod_part in modification.split(';') {
            let mod_part = mod_part.trim();
            if mod_part.is_empty() {
                continue;
            }

            let parts: Vec<&str> = mod_part.split(':').collect();
            if parts.is_empty() {
                continue;
            }

            let mod_name_mass = parts[0].trim();
            let modification_name = Self::get_modification_name(mod_name_mass, glycan_mass, glycan_composition);
            let modification_mass = Self::get_modification_mass(mod_name_mass, glycan_mass, glycan_composition);

            modifications_set.insert(ModificationRow {
                modification_name,
                modification_mass,
            });
        }
    }

    fn get_modification_name(mod_str: &str, glycan_mass: f64, glycan_composition: &str) -> String {
        let parts: Vec<&str> = mod_str.split(' ').collect();
        if parts.is_empty() {
            return mod_str.to_string();
        }

        if let Ok(listed_mass) = parts[0].parse::<f64>() {
            if (listed_mass - glycan_mass).abs() < 0.1 {
                return mod_str.replace(parts[0], glycan_composition);
            }
        }

        mod_str.to_string()
    }

    fn get_modification_mass(mod_str: &str, glycan_mass: f64, glycan_composition: &str) -> String {
        if mod_str.contains(glycan_composition) {
            format!("{:.4}", glycan_mass)
        } else {
            let parts: Vec<&str> = mod_str.split(' ').collect();
            if let Some(mass_str) = parts.first() {
                if let Ok(mass) = mass_str.parse::<f64>() {
                    format!("{:.4}", mass)
                } else {
                    mass_str.to_string()
                }
            } else {
                "0".to_string()
            }
        }
    }

    fn write_outputs(identifications: Vec<IdentificationRow>, modifications: Vec<ModificationRow>, output_dir: &Path) -> Result<()> {
        let identifications_path = output_dir.join("Identifications.csv");
        let mut wtr = csv::Writer::from_path(&identifications_path)?;
        for row in identifications {
            wtr.serialize(row)?;
        }
        wtr.flush()?;

        let modifications_path = output_dir.join("Modifications.csv");
        let mut wtr = csv::Writer::from_path(&modifications_path)?;
        for row in modifications {
            wtr.serialize(row)?;
        }
        wtr.flush()?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::TempDir;

    #[test]
    fn test_glyco_decipher_conversion_against_reference() -> Result<()> {
        let input_path = Path::new(r"C:\Users\vishn\Documents\sample_files\Software\Glyco-Decipher\N\250722_TV_R00011_MG_c00012_SA_glyco_K7144_NT1_GPSM_DatabaseGlycan.txt");
        let reference_dir = Path::new(r"C:\Users\vishn\Documents\sample_files\Software\Glyco-Decipher\N");

        if !input_path.exists() || !reference_dir.exists() {
            eprintln!("Skipping test - reference files not found");
            return Ok(());
        }

        let temp_dir = TempDir::new()?;
        let output_dir = temp_dir.path();

        GlycoDecipherConverter::convert_to_ipsa(input_path, output_dir)?;

        let our_identifications = fs::read_to_string(output_dir.join("Identifications.csv"))?;
        let ref_identifications = fs::read_to_string(reference_dir.join("Identifications.csv"))?;

        let our_lines: Vec<&str> = our_identifications.lines().collect();
        let ref_lines: Vec<&str> = ref_identifications.lines().collect();

        println!("Our output has {} lines, reference has {} lines", our_lines.len(), ref_lines.len());

        for (i, (our_line, ref_line)) in our_lines.iter().zip(ref_lines.iter()).enumerate().take(10) {
            let our_clean = our_line.replace("\"", "");
            let ref_clean = ref_line.replace("\"", "");
            if our_clean != ref_clean {
                println!("Data difference at line {}:", i + 1);
                println!("Our:  {}", our_clean);
                println!("Ref:  {}", ref_clean);
            }
        }

        assert_eq!(our_lines.len(), ref_lines.len(), "Line count mismatch");

        Ok(())
    }
}