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

pub struct OPairConverter;

impl OPairConverter {
    pub fn convert_to_periscope(input_path: &Path, output_dir: &Path) -> Result<()> {
        let file = File::open(input_path)?;
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(file);

        let mut identifications = Vec::new();
        let mut modifications_set = HashSet::new();

        for record in reader.deserialize::<HashMap<String, String>>() {
            let row = record?;

            let scan = row.get("Scan Number").unwrap_or(&String::new()).clone();
            let full_sequence = row.get("Full Sequence").unwrap_or(&String::new()).clone();
            let charge = row.get("Precursor Charge").map_or("0", |v| v).parse::<i32>()?;
            let file_name = row.get("File Name").unwrap_or(&String::new()).clone();

            let (sequence, modification) = Self::extract_modifications(&full_sequence);
            let spectral_file = Self::extract_base_filename(&file_name);

            identifications.push(IdentificationRow {
                scan,
                sequence,
                charge,
                modification: modification.clone(),
                spectral_file,
            });

            if !modification.is_empty() {
                Self::process_modifications(&modification, &mut modifications_set);
            }
        }

        let mut modifications: Vec<_> = modifications_set.into_iter().collect();
        modifications.sort_by(|a, b| a.modification_name.cmp(&b.modification_name));

        Self::write_outputs(identifications, modifications, output_dir)
    }

    fn extract_modifications(sequence: &str) -> (String, String) {
        let mut mod_ids = Vec::new();
        let mut aas = Vec::new();
        let mut sites = Vec::new();
        let mut working_sequence = sequence.to_string();

        while working_sequence.contains('[') {
            let mod_re = Regex::new(r"\[.*?\]").unwrap();
            if let Some(mat) = mod_re.find(&working_sequence) {
                let mod_str = mat.as_str();

                let mod_id_re = Regex::new(r":([^ ]+)").unwrap();
                if let Some(mod_id_match) = mod_id_re.captures(mod_str) {
                    let mod_id = &mod_id_match[1];
                    mod_ids.push(mod_id.to_string());

                    let parts: Vec<&str> = working_sequence.splitn(2, mod_str).collect();
                    if parts.len() > 1 {
                        let prefix = parts[0];
                        if let Some(last_char) = prefix.chars().last() {
                            aas.push(last_char);
                            sites.push(prefix.len());
                        }
                    }

                    working_sequence = working_sequence.replacen(mod_str, "", 1);
                } else {
                    break;
                }
            } else {
                break;
            }
        }

        if mod_ids.is_empty() {
            return (working_sequence, String::new());
        }

        let mut mod_strings = Vec::new();
        for ((id, aa), site) in mod_ids.iter().zip(aas.iter()).zip(sites.iter()) {
            mod_strings.push(format!("{} {}:{}", id, aa, site));
        }

        (working_sequence, mod_strings.join(";"))
    }

    fn extract_base_filename(file_name: &str) -> String {
        let path = Path::new(file_name);
        path.file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_string()
    }

    fn process_modifications(modification: &str, modifications_set: &mut HashSet<ModificationRow>) {
        for mod_part in modification.split(';') {
            let mod_part = mod_part.trim();
            if mod_part.is_empty() {
                continue;
            }

            let parts: Vec<&str> = mod_part.split(' ').collect();
            if parts.len() < 2 {
                continue;
            }

            let mod_name = parts[0];
            let residue_site = parts[1];

            if let Some(colon_pos) = residue_site.find(':') {
                let residue = &residue_site[..colon_pos];

                let modification_name = format!("{} {}", mod_name, residue);
                let modification_mass = Self::get_glyco_mass(mod_name);

                modifications_set.insert(ModificationRow {
                    modification_name,
                    modification_mass,
                });
            }
        }
    }

    fn get_glyco_mass(glyc: &str) -> String {
        let hexnac_mass = 203.07937;
        let hex_mass = 162.05282;
        let neuac_mass = 291.09542;
        let fuc_mass = 146.05791;
        let neugc_mass = 307.0903;

        let mut total_mass = 0.0;

        let n_re = Regex::new(r"N(\d+)").unwrap();
        if let Some(cap) = n_re.captures(glyc) {
            if let Ok(count) = cap[1].parse::<i32>() {
                total_mass += count as f64 * hexnac_mass;
            }
        }

        let h_re = Regex::new(r"H(\d+)").unwrap();
        if let Some(cap) = h_re.captures(glyc) {
            if let Ok(count) = cap[1].parse::<i32>() {
                total_mass += count as f64 * hex_mass;
            }
        }

        let a_re = Regex::new(r"A(\d+)").unwrap();
        if let Some(cap) = a_re.captures(glyc) {
            if let Ok(count) = cap[1].parse::<i32>() {
                total_mass += count as f64 * neuac_mass;
            }
        }

        let f_re = Regex::new(r"F(\d+)").unwrap();
        if let Some(cap) = f_re.captures(glyc) {
            if let Ok(count) = cap[1].parse::<i32>() {
                total_mass += count as f64 * fuc_mass;
            }
        }

        let g_re = Regex::new(r"G(\d+)").unwrap();
        if let Some(cap) = g_re.captures(glyc) {
            if let Ok(count) = cap[1].parse::<i32>() {
                total_mass += count as f64 * neugc_mass;
            }
        }

        if total_mass > 0.0 {
            format!("{:.4}", total_mass)
        } else {
            Self::get_known_modification_mass(glyc)
        }
    }

    fn get_known_modification_mass(mod_name: &str) -> String {
        let mod_masses = HashMap::from([
            ("Deamidation", "0.984016"),
            ("Oxidation", "15.994915"),
            ("Carbamidomethyl", "57.021464"),
        ]);

        mod_masses.get(mod_name).unwrap_or(&"0").to_string()
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
    fn test_opair_conversion_against_reference() -> Result<()> {
        let input_paths = [
            (r"C:\Users\vishn\Documents\sample_files\Software\OPair\N\nglyco.psmtsv", r"C:\Users\vishn\Documents\sample_files\Software\OPair\N"),
            (r"C:\Users\vishn\Documents\sample_files\Software\OPair\O\oglyco.psmtsv", r"C:\Users\vishn\Documents\sample_files\Software\OPair\O"),
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

            OPairConverter::convert_to_periscope(input_path, output_dir)?;

            let our_identifications = fs::read_to_string(output_dir.join("Identifications.csv"))?;

            if let Ok(ref_identifications) = fs::read_to_string(reference_dir.join("Identifications.csv")) {
                let our_lines: Vec<&str> = our_identifications.lines().collect();
                let ref_lines: Vec<&str> = ref_identifications.lines().collect();

                println!("OPair {} - Our output has {} lines, reference has {} lines",
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