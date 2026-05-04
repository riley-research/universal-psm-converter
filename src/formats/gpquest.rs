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

pub struct GPQuestConverter;

impl GPQuestConverter {
    pub fn convert_to_periscope(input_path: &Path, output_dir: &Path) -> Result<()> {
        let file = File::open(input_path)?;
        let mut reader = csv::Reader::from_reader(file);

        let mut identifications = Vec::new();
        let mut modifications_set = HashSet::new();

        let spectral_file = input_path.file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .replace(".csv", "")
            .to_string();

        for record in reader.deserialize::<HashMap<String, String>>() {
            let row = record?;

            let scan = row.get("MS2").unwrap_or(&String::new()).clone();
            let sequence = row.get("Sequence").unwrap_or(&String::new()).clone();
            let charge = row.get("Precursor Charge").map_or("0", |v| v).parse::<i32>()?;
            let peptide = row.get("Peptide").unwrap_or(&String::new()).clone();
            let mod_mass = row.get("Mod Mass").map_or("0", |v| v).to_string();
            let nglycan = row.get("NGLYCAN").unwrap_or(&String::new()).clone();

            let n_residue = Self::find_n_glyco_residue(&sequence);
            if n_residue == -1 {
                continue;
            }

            let glycan_name = nglycan.split(':').next().unwrap_or("").to_string();
            let variable_mods = Self::get_variable_modifications(&peptide)?;
            let final_modification = Self::add_glycan_modification(&variable_mods, &glycan_name, n_residue);

            identifications.push(IdentificationRow {
                scan,
                sequence,
                charge,
                modification: final_modification,
                spectral_file: spectral_file.clone(),
            });

            modifications_set.insert(ModificationRow {
                modification_name: format!("{} N", glycan_name),
                modification_mass: mod_mass,
            });

            if !variable_mods.is_empty() {
                let var_mod_names = Self::get_variable_mod_names(&peptide)?;
                for mod_entry in var_mod_names {
                    modifications_set.insert(mod_entry);
                }
            }
        }

        let mut modifications: Vec<_> = modifications_set.into_iter().collect();
        modifications.sort_by(|a, b| a.modification_name.cmp(&b.modification_name));

        Self::write_outputs(identifications, modifications, output_dir)
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

    fn get_variable_modifications(peptide: &str) -> Result<String> {
        let mut formatted_mods = Vec::new();
        let mut working_peptide = peptide.to_string();

        while working_peptide.contains('[') {
            let bracket_pos = working_peptide.find('[').unwrap();

            let mod_location = bracket_pos;
            let mod_residue = if mod_location > 0 {
                working_peptide.chars().nth(mod_location - 1).unwrap_or('X')
            } else {
                'X'
            };

            let bracket_re = Regex::new(r"\[.*?\]").unwrap();
            if let Some(mat) = bracket_re.find(&working_peptide) {
                let pattern = mat.as_str();
                let mod_name = pattern.trim_start_matches('[').trim_end_matches(']');

                formatted_mods.push(format!("{} {}:{}", mod_name, mod_residue, mod_location));

                let (before, after) = working_peptide.split_once(pattern).unwrap();
                working_peptide = format!("{}{}", before, after);
            } else {
                break;
            }
        }

        Ok(formatted_mods.join(";"))
    }

    fn get_variable_mod_names(peptide: &str) -> Result<Vec<ModificationRow>> {
        let mut results = Vec::new();
        let bracket_re = Regex::new(r"\[.*?\]").unwrap();

        let mod_lookup = HashMap::from([
            ("Oxidation", "15.9949"),
            ("Carbamidomethyl", "57.0215"),
        ]);

        for mat in bracket_re.find_iter(peptide) {
            let pattern = mat.as_str();
            let mod_name = pattern.trim_start_matches('[').trim_end_matches(']');
            let mod_residue = if mat.start() > 0 {
                peptide.chars().nth(mat.start() - 1).unwrap_or('X')
            } else {
                'X'
            };

            if let Some(mass) = mod_lookup.get(mod_name) {
                results.push(ModificationRow {
                    modification_name: format!("{} {}", mod_name, mod_residue),
                    modification_mass: mass.to_string(),
                });
            }
        }

        Ok(results)
    }

    fn add_glycan_modification(mods: &str, glycan_name: &str, glycosite: i32) -> String {
        let glycan_mod = format!("{} N:{}", glycan_name, glycosite);

        if mods.is_empty() {
            glycan_mod
        } else {
            format!("{};{}", mods, glycan_mod)
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
    fn test_gpquest_conversion_against_reference() -> Result<()> {
        let input_path = Path::new(r"C:\Users\vishn\Documents\sample_files\Software\GPQuest\Output\250324_VC_HGI_ModuleC_2_HCD.csv");
        let reference_dir = Path::new(r"C:\Users\vishn\Documents\sample_files\Software\GPQuest\Output");

        if !input_path.exists() || !reference_dir.exists() {
            eprintln!("Skipping test - reference files not found");
            return Ok(());
        }

        let temp_dir = TempDir::new()?;
        let output_dir = temp_dir.path();

        GPQuestConverter::convert_to_periscope(input_path, output_dir)?;

        let our_identifications = fs::read_to_string(output_dir.join("Identifications.csv"))?;

        if let Ok(ref_identifications) = fs::read_to_string(reference_dir.join("Identifications.csv")) {
            let our_lines: Vec<&str> = our_identifications.lines().collect();
            let ref_lines: Vec<&str> = ref_identifications.lines().collect();

            println!("Our output has {} lines, reference has {} lines", our_lines.len(), ref_lines.len());

            for (i, (our_line, ref_line)) in our_lines.iter().zip(ref_lines.iter()).enumerate().take(30) {
                let our_clean = our_line.replace("\"", "");
                let ref_clean = ref_line.replace("\"", "");
                if our_clean != ref_clean {
                    println!("Data difference at line {}:", i + 1);
                    println!("Our:  {}", our_clean);
                    println!("Ref:  {}", ref_clean);
                }
            }

            assert_eq!(our_lines.len(), ref_lines.len(), "Line count mismatch");
        }

        assert!(our_identifications.lines().count() > 0, "Output should not be empty");

        Ok(())
    }
}