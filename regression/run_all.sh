#!/bin/bash

# Script: run_all_tests.sh
# Purpose: Enter each folder in current directory, run CTurboBFM with input.ini,
#          report success/failure.

# Store current directory
BASE_DIR=$(pwd)

# Initialize counters
SUCCESS=0
FAILURE=0

# Loop over all directories in current folder
for dir in */ ; do
    # Remove trailing slash
    dir="${dir%/}"
    
    echo " "
    echo "=== Running test in $dir ==="
    
    # Enter directory
    cd "$dir" || { echo "  Failed to enter $dir"; continue; }
    
    # Run solver
    python main.py > log_test.txt
    EXIT_CODE=$?
    
    # Check exit code
    if [[ $EXIT_CODE -eq 0 ]]; then
        echo "  SUCCESS"
        ((SUCCESS++))
    else
        echo "  FAILURE (exit code $EXIT_CODE)"
        ((FAILURE++))
    fi

    # delete generated files 
    rm -rf Results
    
    # Return to base directory
    cd "$BASE_DIR"
done

# Summary
echo "==============================="
echo "Total directories: $((SUCCESS + FAILURE))"
echo "Succeeded: $SUCCESS"
echo "Failed: $FAILURE"
echo "==============================="

if [[ $FAILURE -eq 0 ]]; then
    exit 0
else
    exit 1
fi
