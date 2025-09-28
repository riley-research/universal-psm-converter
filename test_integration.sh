#!/bin/bash

# Integration testing script for universal-psm-converter
# Run this separately from the main crate to test against local sample files

echo "Running integration tests..."

cargo test --test integration_tests test_all_formats_detect_correctly -- --nocapture

cargo test --test integration_tests test_all_formats_convert_successfully -- --nocapture

cargo test --test integration_tests test_output_format_compatibility -- --nocapture

cargo test --test integration_tests test_error_handling -- --nocapture

echo "Integration testing complete."