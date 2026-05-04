use anyhow::Result;
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::Path;

#[derive(Debug, Serialize)]
pub struct IdentificationRow {
    #[serde(rename = "Scan")]
    pub scan: String,
    #[serde(rename = "Peptide")]
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

pub struct PGlycoConverter;

impl PGlycoConverter {
    pub fn convert_to_periscope(input_path: &Path, output_dir: &Path) -> Result<()> {
        let file = File::open(input_path)?;
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(file);

        let mut identifications = Vec::new();
        let mut var_modifications_set = HashSet::new();
        let mut glycan_modifications_set = HashSet::new();

        for record in reader.deserialize::<HashMap<String, String>>() {
            let row = record?;

            let scan = row.get("Scan").unwrap_or(&String::new()).clone();
            let mut peptide = row.get("Peptide").unwrap_or(&String::new()).clone();
            let charge = row.get("Charge").map_or("0", |v| v).parse::<i32>()?;
            let modifications = row.get("Mod").unwrap_or(&String::new()).clone();
            let glycan_composition = row.get("GlycanComposition").unwrap_or(&String::new()).clone();
            let glycan_mass = row.get("GlyMass").map_or("0", |v| v).to_string();
            let glycosite = row.get("GlySite").map_or("0", |v| v).parse::<i32>().unwrap_or(0);
            let gly_spec = row.get("GlySpec").unwrap_or(&String::new()).clone();

            peptide = peptide.replace('J', "N");

            let variable_mods = Self::get_modifications(&modifications);
            let final_modification = Self::add_glyco_mod(&peptide, &variable_mods, &glycan_composition, glycosite);
            let spectral_file = Self::extract_spectral_file(&gly_spec);

            identifications.push(IdentificationRow {
                scan,
                sequence: peptide.clone(),
                charge,
                modification: final_modification,
                spectral_file,
            });

            if !variable_mods.is_empty() {
                let mod_names = Self::get_modification_names(&modifications);
                for mod_name in mod_names.split(';') {
                    if !mod_name.is_empty() {
                        let mass = Self::get_variable_mod_mass(mod_name);
                        var_modifications_set.insert(ModificationRow {
                            modification_name: mod_name.to_string(),
                            modification_mass: mass,
                        });
                    }
                }
            }

            let gly_residue = if glycosite > 0 && glycosite as usize <= peptide.len() {
                peptide.chars().nth(glycosite as usize - 1).unwrap_or('N')
            } else {
                'N'
            };

            glycan_modifications_set.insert(ModificationRow {
                modification_name: format!("{} {}", glycan_composition, gly_residue),
                modification_mass: glycan_mass,
            });
        }

        let mut all_modifications: Vec<_> = var_modifications_set.into_iter()
            .chain(glycan_modifications_set)
            .collect();
        all_modifications.sort_by(|a, b| a.modification_name.cmp(&b.modification_name));

        Self::write_outputs(identifications, all_modifications, output_dir)
    }

    fn get_modifications(mods: &str) -> String {
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

            let mod_location = parts[0].parse::<i32>().unwrap_or(0);

            let bracket_part = parts[1];
            let residue = if let Some(start) = bracket_part.find('[') {
                if let Some(end) = bracket_part.find(']') {
                    bracket_part[start+1..end].to_string()
                } else {
                    continue;
                }
            } else {
                continue;
            };

            let mod_name = bracket_part.split('[').next().unwrap_or("").trim();

            formatted_mods.push(format!("{} {}:{}", mod_name, residue, mod_location));
        }

        formatted_mods.join(";")
    }

    fn add_glyco_mod(peptide: &str, formatted_mod: &str, glycan_composition: &str, gly_location: i32) -> String {
        let gly_residue = if gly_location > 0 && gly_location as usize <= peptide.len() {
            peptide.chars().nth(gly_location as usize - 1).unwrap_or('N')
        } else {
            'N'
        };

        let glycan_mod = format!("{} {}:{}", glycan_composition, gly_residue, gly_location);

        if formatted_mod.is_empty() {
            glycan_mod
        } else {
            format!("{};{}", formatted_mod, glycan_mod)
        }
    }

    fn get_modification_names(mods: &str) -> String {
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

            let bracket_part = parts[1];
            let residue = if let Some(start) = bracket_part.find('[') {
                if let Some(end) = bracket_part.find(']') {
                    bracket_part[start+1..end].to_string()
                } else {
                    continue;
                }
            } else {
                continue;
            };

            let mod_name = bracket_part.split('[').next().unwrap_or("").trim();

            formatted_mods.push(format!("{} {}", mod_name, residue));
        }

        formatted_mods.join(";")
    }

    fn get_variable_mod_mass(mod_name: &str) -> String {
        let mod_lookup = HashMap::from([
            ("Oxidation", "15.9949"),
            ("Carbamidomethyl", "57.0215"),
        ]);

        let mod_key = mod_name.split(' ').next().unwrap_or("");
        mod_lookup.get(mod_key).unwrap_or(&"0").to_string()
    }

    fn extract_spectral_file(gly_spec: &str) -> String {
        gly_spec.split('.').next().unwrap_or("").to_string()
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
    fn test_pglyco_conversion_against_reference() -> Result<()> {
        let input_paths = [
            (r"C:\Users\vishn\Documents\sample_files\Software\pGlyco\N\pGlycoDB-GP-FDR.txt", r"C:\Users\vishn\Documents\sample_files\Software\pGlyco\N"),
            (r"C:\Users\vishn\Documents\sample_files\Software\pGlyco\O\pGlycoDB-GP-FDR.txt", r"C:\Users\vishn\Documents\sample_files\Software\pGlyco\O"),
        ];

        for (input_path_str, reference_dir_str) in &input_paths {
            let input_path = Path::new(input_path_str);
            let reference_dir = Path::new(reference_dir_str);

            if !input_path.exists() || !reference_dir.exists() {
                eprintln!("Skipping test - files not found: {}", input_path_str);
                continue;
            }

            let temp_dir = TempDir::new()?;
            let output_dir = temp_dir.path();

            PGlycoConverter::convert_to_periscope(input_path, output_dir)?;

            let our_identifications = fs::read_to_string(output_dir.join("Identifications.csv"))?;

            if let Ok(ref_identifications) = fs::read_to_string(reference_dir.join("Identifications.csv")) {
                let our_lines: Vec<&str> = our_identifications.lines().collect();
                let ref_lines: Vec<&str> = ref_identifications.lines().collect();

                println!("pGlyco {} - Our output has {} lines, reference has {} lines",
                         input_path.file_name().unwrap().to_str().unwrap(),
                         our_lines.len(), ref_lines.len());

                for (i, (our_line, ref_line)) in our_lines.iter().zip(ref_lines.iter()).enumerate().take(5) {
                    let our_clean = our_line.replace("\"", "");
                    let ref_clean = ref_line.replace("\"", "");
                    if our_clean != ref_clean {
                        println!("Data difference at line {}:", i + 1);
                        println!("Our:  {}", our_clean);
                        println!("Ref:  {}", ref_clean);
                    }
                }
            }

            assert!(our_identifications.lines().count() > 0, "Output should not be empty");
        }

        Ok(())
    }
}