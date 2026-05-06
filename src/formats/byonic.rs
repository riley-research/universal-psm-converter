use anyhow::{anyhow, Result};
use calamine::{open_workbook, Reader, Xlsx};
use regex::Regex;
use serde::Serialize;
use std::collections::HashSet;
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

#[derive(Debug)]
struct ByonicRow {
    scan: String,
    sequence: String,
    charge: i32,
    glycan: String,
    modifications: String,
    comment: String,
}

pub struct ByonicConverter;

impl ByonicConverter {
    pub fn convert_to_periscope(input_path: &Path, output_dir: &Path) -> Result<()> {
        let extension = input_path.extension()
            .and_then(|ext| ext.to_str())
            .unwrap_or("")
            .to_lowercase();

        let mut rows = match extension.as_str() {
            "xlsx" | "xls" => Self::parse_excel(input_path)?,
            "csv" => Self::parse_csv(input_path)?,
            _ => return Err(anyhow!("Unsupported Byonic file format: {}", extension)),
        };

        let mut identifications = Vec::new();
        let mut modifications_set = HashSet::new();

        for row in &mut rows {
            row.sequence = Self::correct_nterm(&row.sequence);
            row.modifications = Self::get_nterm_residue(&row.sequence, &row.modifications);

            let scan = Self::extract_scan_from_query(&row.scan)?;

            let cleaned_sequence = Self::extract_clean_sequence(&row.sequence);

            let modification = Self::get_mods(&row.glycan, &cleaned_sequence, &row.modifications)?;
            let spectral_file = Self::extract_spectral_file_name(&row.comment);

            let final_clean_sequence = Self::final_clean_sequence(&cleaned_sequence);

            identifications.push(IdentificationRow {
                scan,
                sequence: final_clean_sequence,
                charge: row.charge,
                modification: modification.unwrap_or_default(),
                spectral_file,
            });

            if let Some(clean_mods) = Self::clean_mods(&row.glycan, &cleaned_sequence, &row.modifications)? {
                for mod_str in clean_mods.split(';') {
                    if let Some((name, mass)) = Self::parse_modification_pair(mod_str) {
                        modifications_set.insert(ModificationRow {
                            modification_name: name,
                            modification_mass: mass,
                        });
                    }
                }
            }
        }

        let mut modifications: Vec<_> = modifications_set.into_iter().collect();
        modifications.sort_by(|a, b| a.modification_name.cmp(&b.modification_name));

        Self::write_outputs(identifications, modifications, output_dir)
    }

    fn parse_excel(file_path: &Path) -> Result<Vec<ByonicRow>> {
        let mut workbook: Xlsx<_> = open_workbook(file_path)?;
        let range = workbook.worksheet_range("Spectra")
            .map_err(|e| anyhow!("Could not find Spectra sheet: {}", e))?;

        let headers: Vec<String> = range.rows().next()
            .ok_or_else(|| anyhow!("Empty sheet"))?
            .iter()
            .map(|cell| cell.to_string())
            .collect();

        let scan_idx = Self::find_column_index(&headers, "Scan #")?;
        let peptide_idx = Self::find_column_index(&headers, "Peptide\n< ProteinMetrics Confidential >")?;
        let charge_idx = Self::find_column_index(&headers, "z")?;
        let glycan_idx = Self::find_column_index(&headers, "Glycans\nNHFAGNa").ok();
        let modification_idx = Self::find_column_index(&headers, "Modification Type(s)").ok();
        let comment_idx = Self::find_column_index(&headers, "Comment")?;

        let mut rows = Vec::new();

        for row in range.rows().skip(1) {
            if row.is_empty() {
                continue;
            }

            rows.push(ByonicRow {
                scan: row[scan_idx].to_string(),
                sequence: row[peptide_idx].to_string(),
                charge: row[charge_idx].to_string().parse::<i32>().unwrap_or(0),
                glycan: glycan_idx.map(|idx| row[idx].to_string()).unwrap_or_default(),
                modifications: modification_idx.map(|idx| row[idx].to_string()).unwrap_or_default(),
                comment: row[comment_idx].to_string(),
            });
        }

        Ok(rows)
    }

    fn parse_csv(file_path: &Path) -> Result<Vec<ByonicRow>> {
        let file = File::open(file_path)?;
        let mut reader = csv::Reader::from_reader(file);

        let headers = reader.headers()?.clone();
        let scan_idx = Self::find_column_index_str(&headers, "Scan..")?;
        let peptide_idx = Self::find_column_index_str(&headers, "Peptide...ProteinMetrics.Confidential..")?;
        let charge_idx = Self::find_column_index_str(&headers, "z")?;
        let glycan_idx = Self::find_column_index_str(&headers, "Glycans.NHFAGNa").ok();
        let modification_idx = Self::find_column_index_str(&headers, "Modification.Type.s.").ok();
        let comment_idx = Self::find_column_index_str(&headers, "Comment")?;

        let mut rows = Vec::new();

        for record in reader.records() {
            let record = record?;

            rows.push(ByonicRow {
                scan: record[scan_idx].to_string(),
                sequence: record[peptide_idx].to_string(),
                charge: record[charge_idx].parse::<i32>().unwrap_or(0),
                glycan: glycan_idx.map(|idx| record[idx].to_string()).unwrap_or_default(),
                modifications: modification_idx.map(|idx| record[idx].to_string()).unwrap_or_default(),
                comment: record[comment_idx].to_string(),
            });
        }

        Ok(rows)
    }

    fn correct_nterm(seq: &str) -> String {
        if seq.contains(".[") && seq.contains("].") {
            let re = Regex::new(r"\.(\[.*?\])\.").unwrap();
            if let Some(cap) = re.captures(seq) {
                let mod_part = &cap[1];
                let seq_without_nterm_mod = seq.replace(&format!(".{}.", mod_part), ".");
                if seq_without_nterm_mod.len() >= 4 {
                    let pt1 = &seq_without_nterm_mod[..3];
                    let pt2 = &seq_without_nterm_mod[3..];
                    return format!("{}{}{}", pt1, mod_part, pt2);
                }
            }
        }
        seq.to_string()
    }

    fn get_nterm_residue(seq: &str, mods: &str) -> String {
        if seq.len() >= 3 {
            let nterm_aa = seq.chars().nth(2).unwrap_or('X');
            mods.replace('.', &nterm_aa.to_string())
        } else {
            mods.to_string()
        }
    }

    fn extract_scan_from_query(scan_str: &str) -> Result<String> {
        if let Some(parts) = scan_str.split("=").last() {
            Ok(parts.to_string())
        } else {
            let re = Regex::new(r"(\d+)").unwrap();
            if let Some(captures) = re.captures(scan_str) {
                Ok(captures[1].to_string())
            } else {
                Err(anyhow!("Could not extract scan number from: {}", scan_str))
            }
        }
    }

    fn extract_clean_sequence(seq: &str) -> String {
        if seq.len() >= 4 && seq.chars().nth(1) == Some('.') && seq.chars().nth(seq.len()-2) == Some('.') {
            seq[2..seq.len()-2].to_string()
        } else {
            seq.to_string()
        }
    }

    fn final_clean_sequence(seq: &str) -> String {
        let re = Regex::new(r"\[.*?\]|\d+|\+|\.").unwrap();
        re.replace_all(seq, "").to_string()
    }

    fn extract_spectral_file_name(comment: &str) -> String {
        if let Some(parts) = comment.split(".ScanId").next() {
            parts.to_string()
        } else {
            comment.to_string()
        }
    }

    fn get_mods(glycan_mod: &str, cleaned_seq: &str, mods: &str) -> Result<Option<String>> {
        if mods.is_empty() {
            return Ok(None);
        }

        let mut formatted_mods = Vec::new();
        let mod_parts: Vec<&str> = cleaned_seq.split(']').collect();

        for (i, part) in mod_parts.iter().enumerate() {
            if !part.contains('[') {
                continue;
            }

            let bracket_start = part.find('[').ok_or_else(|| anyhow!("Invalid modification format"))?;
            let modi_substr = if bracket_start > 0 {
                &part[bracket_start-1..]
            } else {
                part
            };

            let residue = modi_substr.chars().next().unwrap_or('X').to_string();

            let full_pep: String = mod_parts[..=i].join("").chars().filter(|c| c.is_ascii_uppercase()).collect();
            let localization = full_pep.len();

            let mod_content_start = part.find('[').unwrap() + 1;
            let mod_mass_str = &part[mod_content_start..];
            let mod_mass_str = mod_mass_str.trim_start_matches('+');
            let mod_mass = mod_mass_str.parse::<f64>().unwrap_or(0.0);

            let modification = if !glycan_mod.is_empty() {
                Self::match_glycan_modification(glycan_mod, mod_mass).unwrap_or_else(|| mod_mass_str.to_string())
            } else {
                mod_mass_str.to_string()
            };

            formatted_mods.push(format!("{} {}:{}", modification, residue, localization));
        }

        if formatted_mods.is_empty() {
            Ok(None)
        } else {
            Ok(Some(formatted_mods.join(";")))
        }
    }

    fn clean_mods(glycan_mod: &str, cleaned_seq: &str, mods: &str) -> Result<Option<String>> {
        if mods.is_empty() {
            return Ok(None);
        }

        let mut formatted_mods = Vec::new();
        let mod_parts: Vec<&str> = cleaned_seq.split(']').collect();

        for part in mod_parts.iter() {
            if !part.contains('[') {
                continue;
            }

            let bracket_start = part.find('[').ok_or_else(|| anyhow!("Invalid modification format"))?;
            let modi_substr = if bracket_start > 0 {
                &part[bracket_start-1..]
            } else {
                part
            };

            let residue = modi_substr.chars().next().unwrap_or('X').to_string();

            let mod_content_start = part.find('[').unwrap() + 1;
            let mod_mass_str = &part[mod_content_start..];
            let mod_mass_str = mod_mass_str.trim_start_matches('+');
            let mod_mass = mod_mass_str.parse::<f64>().unwrap_or(0.0);

            let modification = if !glycan_mod.is_empty() {
                Self::match_glycan_modification(glycan_mod, mod_mass).unwrap_or_else(|| mod_mass.to_string())
            } else {
                mod_mass.to_string()
            };

            formatted_mods.push(format!("{} {},{}", modification, residue, mod_mass_str));
        }

        if formatted_mods.is_empty() {
            Ok(None)
        } else {
            Ok(Some(formatted_mods.join(";")))
        }
    }

    fn match_glycan_modification(glycan_mod: &str, mass: f64) -> Option<String> {
        let glycan_masses = Self::calculate_glycan_masses(glycan_mod);

        for (composition, calc_mass) in glycan_masses {
            if (calc_mass - mass).abs() < 0.1 {
                return Some(composition);
            }
        }
        None
    }

    fn calculate_glycan_masses(glycan_mod: &str) -> Vec<(String, f64)> {
        let mut results = Vec::new();

        let hexnac_mass = 203.07937;
        let hex_mass = 162.05282;
        let neuac_mass = 291.09542;
        let fuc_mass = 146.05791;
        let neugc_mass = 307.0903;

        for glycan in glycan_mod.split(',') {
            let glycan = glycan.trim();
            if glycan.is_empty() {
                continue;
            }

            let mut total_mass = 0.0;

            if let Some(cap) = Regex::new(r"HexNAc\((\d+)\)").unwrap().captures(glycan) {
                if let Ok(count) = cap[1].parse::<i32>() {
                    total_mass += count as f64 * hexnac_mass;
                }
            }

            if let Some(cap) = Regex::new(r"Hex\((\d+)\)").unwrap().captures(glycan) {
                if let Ok(count) = cap[1].parse::<i32>() {
                    total_mass += count as f64 * hex_mass;
                }
            }

            if let Some(cap) = Regex::new(r"NeuAc\((\d+)\)").unwrap().captures(glycan) {
                if let Ok(count) = cap[1].parse::<i32>() {
                    total_mass += count as f64 * neuac_mass;
                }
            }

            if let Some(cap) = Regex::new(r"Fuc\((\d+)\)").unwrap().captures(glycan) {
                if let Ok(count) = cap[1].parse::<i32>() {
                    total_mass += count as f64 * fuc_mass;
                }
            }

            if let Some(cap) = Regex::new(r"NeuGc\((\d+)\)").unwrap().captures(glycan) {
                if let Ok(count) = cap[1].parse::<i32>() {
                    total_mass += count as f64 * neugc_mass;
                }
            }

            results.push((glycan.to_string(), total_mass));
        }

        results
    }

    fn parse_modification_pair(mod_str: &str) -> Option<(String, String)> {
        if let Some(comma_pos) = mod_str.rfind(',') {
            let name = mod_str[..comma_pos].trim().to_string();
            let mass = mod_str[comma_pos+1..].trim().to_string();
            Some((name, mass))
        } else {
            None
        }
    }

    fn find_column_index(headers: &[String], target: &str) -> Result<usize> {
        headers.iter().position(|h| h == target)
            .ok_or_else(|| anyhow!("Column '{}' not found", target))
    }

    fn find_column_index_str(headers: &csv::StringRecord, target: &str) -> Result<usize> {
        headers.iter().position(|h| h == target)
            .ok_or_else(|| anyhow!("Column '{}' not found", target))
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
    fn test_byonic_conversion_against_reference() -> Result<()> {
        let input_path = Path::new(r"C:\Users\vishn\Documents\sample_files\Software\Byonic\N_Trypsin.xlsx");
        let reference_dir = Path::new(r"C:\Users\vishn\Documents\sample_files\Software\Byonic\original");

        if !input_path.exists() || !reference_dir.exists() {
            eprintln!("Skipping test - reference files not found");
            return Ok(());
        }

        let temp_dir = TempDir::new()?;
        let output_dir = temp_dir.path();

        ByonicConverter::convert_to_periscope(input_path, output_dir)?;

        let our_identifications = fs::read_to_string(output_dir.join("Identifications.csv"))?;
        let ref_identifications = fs::read_to_string(reference_dir.join("Identifications.csv"))?;

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

        Ok(())
    }
}