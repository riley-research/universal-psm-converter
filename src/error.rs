use thiserror::Error;

#[derive(Error, Debug)]
pub enum ConversionError {
    #[error("Unsupported file format")]
    UnsupportedFormat,

    #[error("File not found: {path}")]
    FileNotFound { path: String },

    #[error("Invalid file structure: {message}")]
    InvalidStructure { message: String },

    #[error("Missing required column: {column}")]
    MissingColumn { column: String },

    #[error("Parse error at row {row}: {message}")]
    ParseError { row: usize, message: String },

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("CSV error: {0}")]
    CsvError(#[from] csv::Error),

    #[error("Excel error: {0}")]
    ExcelError(#[from] calamine::Error),

    #[error("Regex error: {0}")]
    RegexError(#[from] regex::Error),

    #[error("Conversion failed: {message}")]
    ConversionFailed { message: String },
}