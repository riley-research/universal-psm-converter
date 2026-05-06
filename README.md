# Universal PSM Converter

Universal converter for glycoproteomics software outputs to PERISCOPE format.

## Supported Formats

- **FragPipe** (.tsv) - MSFragger glycoproteomics output
- **Byonic** (.xlsx, .csv) - Protein Metrics Byonic output
- **Glyco-Decipher** (.txt) - DatabaseGlycan and PSM formats
- **GPQuest** (.csv) - GPQuest glycoproteomics output
- **OPair** (.psmtsv) - OPair N-glycan and O-glycan output
- **pGlyco** (.txt) - pGlyco3.0 output files

## Usage

### As Library

```rust
use universal_psm_converter::{convert_file_to_periscope, detect_format_from_path};

// Detect format
let format = detect_format_from_path(&input_path)?;
println!("Detected: {}", format);

// Convert to PERISCOPE format
let stats = convert_file_to_periscope(&input_path, &output_dir)?;
println!("Converted {} identifications", stats.identifications_written);
```

### As CLI Tool

```bash
cargo run --features cli -- --input sample.tsv --output ./results/
```

## Output Format

All converters produce standardized PERISCOPE CSV files:

**Identifications.csv:**
```csv
Scan,Sequence,Charge,Modification,spectral_file
2124,KQLGPGKK,2,,EThcD_Glyco_1
2162,STNHEPSEMSNR,3,15.99 M:9;HexNAc(2)Hex(8) N:11,EThcD_Glyco_1
```

**Modifications.csv:**
```csv
Modification Name,Modification Mass
15.99 M,15.9949
HexNAc(2)Hex(8) N,1702.5813
```

## Testing

```bash
cargo test -- --nocapture
```