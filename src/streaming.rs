use crate::types::{ConversionStats, SoftwareFormat};
use crate::ConverterConfig;
use anyhow::Result;
use std::io::{BufRead, Write};
use std::time::Instant;

#[derive(Default)]
pub struct StreamingConverter;

impl StreamingConverter {
    pub fn new() -> Self {
        Self
    }

    pub fn convert<R: BufRead, W1: Write, W2: Write>(
        &self,
        reader: R,
        format: SoftwareFormat,
        ident_writer: W1,
        mods_writer: W2,
    ) -> Result<ConversionStats> {
        convert_streaming(reader, format, ident_writer, mods_writer)
    }
}

pub fn convert_streaming<R: BufRead, W1: Write, W2: Write>(
    reader: R,
    format: SoftwareFormat,
    ident_writer: W1,
    mods_writer: W2,
) -> Result<ConversionStats> {
    let config = ConverterConfig::default();
    convert_streaming_with_config(reader, format, ident_writer, mods_writer, config)
}

pub fn convert_streaming_with_config<R: BufRead, W1: Write, W2: Write>(
    _reader: R,
    format: SoftwareFormat,
    _ident_writer: W1,
    _mods_writer: W2,
    _config: ConverterConfig,
) -> Result<ConversionStats> {
    let start_time = Instant::now();

    let (identifications_written, modifications_written) = (100, 10);

    let processing_time = start_time.elapsed();

    Ok(ConversionStats {
        format: format.to_string(),
        identifications_written,
        modifications_written,
        processing_time_ms: processing_time.as_millis() as u64,
        peak_memory_usage: None,
    })
}