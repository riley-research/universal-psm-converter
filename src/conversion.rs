use crate::detection::detect_format;
use crate::types::{ConversionStats, IdentificationRecord, ModificationRecord, PGlycoIdentificationRecord, SoftwareFormat};
use crate::{ConverterConfig, formats};
use anyhow::Result;
use std::fs;
use std::path::Path;
use std::time::Instant;

pub fn convert_file_to_periscope(input: &Path, output_dir: &Path) -> Result<ConversionStats> {
    let config = ConverterConfig::default();
    convert_file_to_periscope_with_config(input, output_dir, config)
}

pub fn convert_file_to_periscope_with_config(
    input: &Path,
    output_dir: &Path,
    _config: ConverterConfig,
) -> Result<ConversionStats> {
    let format = detect_format(input)?;
    let start_time = Instant::now();

    match format {
        SoftwareFormat::FragPipe => {
            formats::fragpipe::convert_fragpipe_to_periscope(input, output_dir)?;
        },
        SoftwareFormat::Byonic => {
            use crate::formats::byonic::ByonicConverter;
            ByonicConverter::convert_to_periscope(input, output_dir)?;
        },
        SoftwareFormat::GlycoDecipher => {
            use crate::formats::glyco_decipher::GlycoDecipherConverter;
            GlycoDecipherConverter::convert_to_periscope(input, output_dir)?;
        },
        SoftwareFormat::GPQuest => {
            use crate::formats::gpquest::GPQuestConverter;
            GPQuestConverter::convert_to_periscope(input, output_dir)?;
        },
        SoftwareFormat::OPair => {
            use crate::formats::opair::OPairConverter;
            OPairConverter::convert_to_periscope(input, output_dir)?;
        },
        SoftwareFormat::PGlyco => {
            use crate::formats::pglyco::PGlycoConverter;
            PGlycoConverter::convert_to_periscope(input, output_dir)?;
        },
        SoftwareFormat::hgi => {
            crate::formats::hgi::convert_HGI_to_periscope(input, output_dir)?;
        },
    }

    let processing_time = start_time.elapsed();

    let (identifications_written, modifications_written) = count_actual_output(output_dir)?;

    Ok(ConversionStats {
        format: format.to_string(),
        identifications_written,
        modifications_written,
        processing_time_ms: processing_time.as_millis() as u64,
        peak_memory_usage: None,
    })
}

fn count_actual_output(output_dir: &Path) -> Result<(usize, usize)> {
    let ident_path = output_dir.join("Identifications.csv");
    let mods_path = output_dir.join("Modifications.csv");

    let ident_count = if ident_path.exists() {
        let content = fs::read_to_string(&ident_path)?;
        content.lines().count().saturating_sub(1)
    } else {
        0
    };

    let mods_count = if mods_path.exists() {
        let content = fs::read_to_string(&mods_path)?;
        content.lines().count().saturating_sub(1)
    } else {
        0
    };

    Ok((ident_count, mods_count))
}

pub fn convert_to_memory(_input: &Path) -> Result<(Vec<IdentificationRecord>, Vec<ModificationRecord>)> {
    Ok((vec![], vec![]))
}

pub fn convert_to_pglyco_memory(_input: &Path) -> Result<(Vec<PGlycoIdentificationRecord>, Vec<ModificationRecord>)> {
    Ok((vec![], vec![]))
}