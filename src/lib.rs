use anyhow::Result;
use std::io::{BufRead, Write};
use std::path::Path;

pub mod error;
pub mod formats;
pub mod types;
pub mod detection;
pub mod streaming;
pub mod conversion;

pub use error::ConversionError;
pub use types::{SoftwareFormat, ConversionStats, IdentificationRecord, ModificationRecord, PGlycoIdentificationRecord};
pub use detection::detect_format;
pub use streaming::{convert_streaming, StreamingConverter};
pub use conversion::{convert_file_to_periscope, convert_to_memory};

#[derive(Debug, Clone)]
pub struct ConverterConfig {
    pub skip_empty_modifications: bool,
    pub validate_output: bool,
    pub buffer_size: usize,
}

impl Default for ConverterConfig {
    fn default() -> Self {
        Self {
            skip_empty_modifications: true,
            validate_output: true,
            buffer_size: 8192,
        }
    }
}

pub fn detect_format_from_path(path: &Path) -> Result<SoftwareFormat> {
    detection::detect_format(path)
}

pub fn detect_format_from_reader<R: BufRead>(reader: R) -> Result<SoftwareFormat> {
    detection::detect_format_from_reader(reader)
}

pub fn convert_streaming_with_config<R: BufRead, W1: Write, W2: Write>(
    reader: R,
    format: SoftwareFormat,
    ident_writer: W1,
    mods_writer: W2,
    config: ConverterConfig,
) -> Result<ConversionStats> {
    streaming::convert_streaming_with_config(reader, format, ident_writer, mods_writer, config)
}

pub fn convert_file_to_periscope_with_config(
    input: &Path,
    output_dir: &Path,
    config: ConverterConfig,
) -> Result<ConversionStats> {
    conversion::convert_file_to_periscope_with_config(input, output_dir, config)
}

#[cfg(feature = "cli")]
pub mod cli;


#[cfg(test)]
mod lib_tests {
    use super::*;
    use std::path::Path;
    use crate::formats::fragpipe::convert_fragpipe_to_periscope;
    
    #[test]
    fn test_convert_fragpipe_to_periscope() {
        let input = Path::new("C:/PERISCOPE_testfiles/FP/psm.tsv");
        let output = Path::new("C:/PERISCOPE_testfiles/FP/");

        let result = convert_fragpipe_to_periscope(input, output);

        assert!(result.is_ok());
    }

    #[test]
    fn test_software_format_display() {
        assert_eq!(SoftwareFormat::FragPipe.to_string(), "FragPipe");
        assert_eq!(SoftwareFormat::Byonic.to_string(), "Byonic");
        assert_eq!(SoftwareFormat::GlycoDecipher.to_string(), "Glyco-Decipher");
        assert_eq!(SoftwareFormat::GPQuest.to_string(), "GPQuest");
        assert_eq!(SoftwareFormat::OPair.to_string(), "OPair");
        assert_eq!(SoftwareFormat::PGlyco.to_string(), "pGlyco");
    }

    #[test]
    fn test_converter_config_default() {
        let config = ConverterConfig::default();
        assert!(config.skip_empty_modifications);
        assert!(config.validate_output);
        assert_eq!(config.buffer_size, 8192);
    }
    
        
}