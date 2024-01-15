#!/bin/bash
# Execute run_multi_signal.sh, run_multi_efficiency.sh, and run_multi_spectra.sh

# Signal
echo "Running run_multi_signal.sh"
scripts/run_multi_signal.sh
echo "Finished running run_multi_signal.sh"

# Efficiency
echo "Running run_multi_efficiency.sh"
scripts/run_multi_efficiency.sh
echo "Finished running run_multi_efficiency.sh"

# Spectra
echo "Running run_multi_spectra.sh"
scripts/run_multi_spectra.sh
echo "Finished running run_multi_spectra.sh"