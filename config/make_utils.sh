#!/bin/bash

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Configuration Variables ---
# This is the single, fixed output file your C++ code will read from.
OUTPUT_CONFIG="guiutils_current.txt"

# Define the paths for your campaign-specific input files containing the GUIID
GUIID_NOM_2020="guiutils_nominal2020.txt"
GUIID_NOM_2025="guiutils_nominal2025.txt"
GUIID_EXT_2020="guiutils_extracted2020.txt"
GUIID_EXT_2025="guiutils_extracted2025.txt"

# --- Script Logic ---

# Check if a campaign argument was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 [Nominal2020|Nominal2025|Extracted2020|Extracted2025]"
    echo "Example: $0 Nominal2025"
    exit 1
fi 

CAMPAIGN=$1

# 1. Ensure the output directory exists
mkdir -p $(dirname "$OUTPUT_CONFIG")

# 2. Determine which source file to use based on the campaign argument
SOURCE_FILE=""
if [ "$CAMPAIGN" == "Nominal2020" ]; then
    SOURCE_FILE="$GUIID_NOM_2020"
    echo "Selected configuration for MDC2020 Nominal"
elif [ "$CAMPAIGN" == "Nominal2025" ]; then
    SOURCE_FILE="$GUIID_NOM_2025"
    echo "Selected configuration for MDC2025 Nominal"
elif [ "$CAMPAIGN" == "Extracted2020" ]; then
    SOURCE_FILE="$GUIID_EXT_2020"
    echo "Selected configuration for Extracted MDC2020"
elif [ "$CAMPAIGN" == "Extracted2025" ]; then
    SOURCE_FILE="$GUIID_EXT_2025"
    echo "Selected configuration for Extracted MDC2025"
else
    echo "Error: Invalid campaign specified. Use 'Nominal2020', 'Nominal2025', 'Extracted2020', or 'Extracted2025'."
    exit 1
fi

# 3. Copy the selected source file to the fixed output file name
# This effectively makes 'guiutils_current.txt' an alias for the selected config file
cp "$SOURCE_FILE" "$OUTPUT_CONFIG"

echo "Successfully prepared $OUTPUT_CONFIG from $SOURCE_FILE"

