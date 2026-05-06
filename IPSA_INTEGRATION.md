# PERISCOPE Integration Guide

## Add Dependency

Add to PERISCOPE's `src-tauri/Cargo.toml`:

```toml
[dependencies]
universal-psm-converter = "0.1.0"
```

## Replace Current Parsers

### Remove Old Code

Remove from `src-tauri/src/lib.rs`:
- `check_fragpipe_file`
- `parse_fragpipe_file`
- `FragPipeResult` struct

### Add Universal Converter

```rust
#[tauri::command]
fn convert_universal_format(
    app: tauri::AppHandle,
    state: State<'_, AppState>,
    file_path: String,
) -> Result<UniversalConversionResult, String> {
    use universal_psm_converter::{detect_format_from_path, convert_file_to_periscope};
    use std::path::Path;

    let input_path = Path::new(&file_path);
    let format = detect_format_from_path(input_path)
        .map_err(|e| format!("Unsupported format: {}", e))?;

    let temp_dir = tempfile::tempdir().map_err(|e| e.to_string())?;
    let stats = convert_file_to_periscope(input_path, temp_dir.path())
        .map_err(|e| format!("Conversion failed: {}", e))?;

    let identifications_path = temp_dir.path().join("Identifications.csv");
    let modifications_path = temp_dir.path().join("Modifications.csv");

    let ident_result = upload_identification_csv(
        app.clone(), state.clone(),
        identifications_path.to_string_lossy().to_string(),
        None
    )?;

    let mods_result = upload_modification_csv(
        state.clone(),
        modifications_path.to_string_lossy().to_string()
    )?;

    if ident_result.file_id > 0 && mods_result.file_id > 0 {
        link_modification_files(state, ModificationLinkRequest {
            modification_file_ids: vec![mods_result.file_id],
            identification_file_id: ident_result.file_id,
            priorities: vec![],
        })?;
    }

    Ok(UniversalConversionResult {
        format: format.to_string(),
        identification_file_id: ident_result.file_id,
        modification_file_id: mods_result.file_id,
        stats: stats.identifications_written,
    })
}

#[derive(Debug, Serialize)]
pub struct UniversalConversionResult {
    pub format: String,
    pub identification_file_id: i64,
    pub modification_file_id: i64,
    pub stats: usize,
}
```

## Update UI

Replace format-specific logic in `FileBrowser.vue`:

```typescript
// Replace parseFragPipeFile and format-specific checks
const handleUniversalFile = async (filePath: string) => {
  const extension = filePath.toLowerCase().split('.').pop();

  if (['tsv', 'xlsx', 'txt', 'psmtsv'].includes(extension)) {
    const toastId = toast.loading(`Converting ${extension.toUpperCase()} file...`);

    try {
      const result = await invoke('convert_universal_format', { filePath });
      toast.success(`Converted ${result.format}: ${result.stats} identifications`, { id: toastId });
      await loadDatabaseFiles();
    } catch (error) {
      toast.error(`Conversion failed: ${error}`, { id: toastId });
    }
  } else {
    // Existing CSV/raw file handling
    await handleExistingFileTypes(filePath);
  }
};
```

## Benefits

- **6 formats supported** vs current 1 (FragPipe only)
- **Cleaner codebase** - format logic separated
- **Better performance** - optimized parsers
- **Community reusable** - separate crate