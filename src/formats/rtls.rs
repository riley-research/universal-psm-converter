use std::{collections::{HashMap, BTreeSet}, path::Path, fs::File};
use anyhow::{Context, Result};
use super::output_rows::{IdentificationRow, ModificationRow};


#[derive(Debug)]
struct RTLSRow {
    scan_number: u32,
    sequence: Option<String>,
    precursor_mz: f64,
    mz_diff: f64,
    charge: i32,
    mods: Option<String>,
    capital_y_ions: Option<String>,
    file_name: String,
    extra_columns: HashMap<String, String>,
}

pub fn convert_RTLS_to_periscope(input_path: &Path, output_dir: &Path) -> Result<()> {
    let file = File::open(input_path)?;
    let mut reader = csv::ReaderBuilder::new()
        .from_reader(file);

    let mut identifications = Vec::new();
    let mut modifications = Vec::new();

    for result in reader.deserialize::<HashMap<String, String>>() {
        let row = row_from_record(result?)?;

        //Get the modifications ready
        let mods = clean_mods(row.mods.as_deref());
        for m in mods {
            if !modifications.contains(&m){
                modifications.push(m);
            }
        }

        let sequence = row.sequence.clone().unwrap_or_default();
        let modification = format_modifications(row.mods.as_deref(), &sequence);

        //Get the identifications ready
        let scan_number: String = row.scan_number.to_string();
        let charge = row.charge;
        let spectrum_file = row.file_name;
        let extra_columns = row.extra_columns;

        identifications.push(IdentificationRow {
            scan: scan_number,
            sequence: sequence,
            charge: charge,
            modification: modification,
            spectral_file: spectrum_file,
            extra_columns: Option::from(extra_columns),
        });

    }

        write_identifications(output_dir, identifications)?;
        write_modification(output_dir, modifications)?;

        Ok(())
}

fn get_required(record: &mut HashMap<String, String>, key: &str) -> Result<String> {
    record
        .remove(key)
        .with_context(|| {
            format!(
                "Missing required column '{}'. Available columns: {:?}",
                key,
                record.keys().collect::<Vec<_>>()
            )
        })
}

fn row_from_record(mut record: HashMap<String, String>) -> Result<RTLSRow> {
    let scan_number = get_required(&mut record, "ScanNumber")?
        .parse()
        .unwrap_or(0);
    let sequence = record.remove("Sequence");
    let precursor_mz =  get_required(&mut record, "PrecursorMZ")?
        .parse()
        .unwrap_or(0.0);
    let mz_diff = get_required(&mut record, "MZDiff(exp-lib)")?
        .parse()
        .unwrap_or(0.0);
    let charge: i32 = get_required(&mut record, "Charge")?
        .parse()
        .unwrap_or(0);

    let mods  = record.remove("Mods");
    let file_name = record.remove("FileName").unwrap_or_default();
    let capital_y_ions =record.remove("Yions");
    let extra_columns = record;

    Ok(RTLSRow{
        scan_number,
        sequence,
        precursor_mz,
        mz_diff,
        charge,
        mods,
        capital_y_ions,
        file_name,
        extra_columns,  
    })
}

fn format_modifications(mods: Option<&str>, sequence: &str) -> String {
    mods.unwrap_or("")
        .split(';')
        .filter_map(|entry| {
            let entry = entry.trim();
            if entry.is_empty() {
                return None;
            }

            let mut parts = entry.splitn(2, '@');
            let mass = parts.next()?.trim();
            let position: usize = parts.next()?.trim().parse().ok()?;

            let residue = sequence.chars().nth(position.checked_sub(1)?)?;

            Some(format!("{} {}:{}", mass, residue, position))
        })
        .collect::<Vec<_>>()
        .join(";")
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
            let mod_info = modification.split('@').collect::<Vec<&str>>()[0];

            let name:String;
            let mass: f64 ;
            if mod_info == "NGlycan" {
                name = "NGlycan".to_string();
                mass = 203.07937;
            }
            else{
                name =  mod_info.to_string();
                mass = mod_info.parse().unwrap_or(0.0);
            }
        
            Some(ModificationRow {
                modification_name: name.trim().to_string(),
                modification_mass: mass.to_string(),
            })
        })
        .collect()
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
        .flat_map(|r| r.extra_columns.as_ref().unwrap().keys().cloned())
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
        if row.sequence.is_empty(){
            continue; // Skip rows with empty or None sequence
        }
        let mut record = vec![
            row.scan.to_string(),
            row.sequence.clone(),
            row.charge.to_string(),
            row.modification.clone(),
            row.spectral_file.clone(),
        ];
        for k in &extra_keys {
            record.push(row.extra_columns.as_ref().unwrap().get(k).cloned().unwrap_or_default());
        }
        wtr.write_record(&record)?;
    }
    wtr.flush()?;
    
    Ok(())
}



#[cfg(test)]

mod lib_tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_convert_rtls_to_periscope() {

    let input = Path::new(r"Z:\R00016_RTLS\RawFiles\26-09-11_SerumRTLS_searches\GlycoEnriched\search.csv");
    let output = Path::new(r"Z:\R00016_RTLS\RawFiles\26-09-11_SerumRTLS_searches\GlycoEnriched\");

    let result = convert_RTLS_to_periscope(input, output);

    match result {

        Ok(_) => println!("Success"),

        Err(e) => panic!("Conversion failed: {}", e),

    }

}

}