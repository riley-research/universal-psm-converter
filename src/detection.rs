use crate::error::ConversionError;
use crate::types::SoftwareFormat;
use anyhow::Result;
use calamine::{open_workbook, Reader, Xlsx};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub fn detect_format_from_path(path: &Path) -> Result<SoftwareFormat> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    detect_format_from_reader(reader)
}

pub fn detect_format_from_reader<R: BufRead>(mut reader: R) -> Result<SoftwareFormat> {
    let mut first_line = String::new();
    reader.read_line(&mut first_line)?;

    let headers = first_line.to_lowercase();

    if headers.contains("spectrum") && headers.contains("total glycan composition") && headers.contains("assigned modifications") {
        return Ok(SoftwareFormat::FragPipe);
    }

    if headers.contains("glysite") && headers.contains("glyspec") && headers.contains("glycancomposition") {
        return Ok(SoftwareFormat::PGlyco);
    }

    if headers.contains("scan number") && headers.contains("full sequence") && headers.contains("plausible glycancomposition") {
        return Ok(SoftwareFormat::OPair);
    }

    if headers.contains("glycancomposition") || headers.contains("glycanexpmas") {
        return Ok(SoftwareFormat::GlycoDecipher);
    }

    if headers.contains("ms2") && headers.contains("nglycan") && headers.contains("mod mass") {
        return Ok(SoftwareFormat::GPQuest);
    }

    if headers.contains("Data File") && headers.contains("Identified peptide base sequence") && headers.contains("Observed charge state (z)") {
        return Ok(SoftwareFormat::hgi);
    }

    Err(ConversionError::UnsupportedFormat.into())
}

pub fn detect_format_from_excel(path: &Path) -> Result<SoftwareFormat> {
    let workbook: Xlsx<_> = open_workbook(path)?;
    let sheet_names = workbook.sheet_names().to_owned();

    if sheet_names.contains(&"Spectra".to_string()) {
        return Ok(SoftwareFormat::Byonic);
    }

    Err(ConversionError::UnsupportedFormat.into())
}

pub fn detect_format(path: &Path) -> Result<SoftwareFormat> {
    let extension = path.extension()
        .and_then(|ext| ext.to_str())
        .unwrap_or("")
        .to_lowercase();

    match extension.as_str() {
        "xlsx" | "xls" => detect_format_from_excel(path),
        "tsv" | "txt" | "psmtsv" => {
            let file = File::open(path)?;
            let reader = BufReader::new(file);
            detect_format_from_reader(reader)
        },
        "csv" => {
            let file = File::open(path)?;
            let mut reader = csv::Reader::from_reader(file);

            if let Ok(headers) = reader.headers() {
                let header_str = headers.iter().collect::<Vec<_>>().join(",").to_lowercase();

                if header_str.contains("scan..") && header_str.contains("peptide...proteinmetrics") {
                    return Ok(SoftwareFormat::Byonic);
                }

                if header_str.contains("ms2") && header_str.contains("nglycan") && header_str.contains("mod mass") {
                    return Ok(SoftwareFormat::GPQuest);
                }
            }

            Err(ConversionError::UnsupportedFormat.into())
        },
        _ => Err(ConversionError::UnsupportedFormat.into()),
    }
}