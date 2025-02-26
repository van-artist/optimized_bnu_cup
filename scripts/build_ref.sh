#!/bin/bash

# Get the input path parameters
HISAT3N_BUILD=$1
GENOME_FA=$2
GENOME_INDEX=$3
NCRNA_FA=$4
NCRNA_INDEX=$5

# Check if all parameters are provided
if [ -z "$HISAT3N_BUILD" ] || [ -z "$GENOME_FA" ] || [ -z "$GENOME_INDEX" ] || [ -z "$NCRNA_FA" ] || [ -z "$NCRNA_INDEX" ]; then
    echo "=> Error: Missing required parameters! Please provide the following paths:"
    echo "1. HISAT3N_BUILD path"
    echo "2. GENOME_FA path"
    echo "3. GENOME_INDEX path"
    echo "4. NCRNA_FA path"
    echo "5. NCRNA_INDEX path"
    exit 1
fi

# Ensure the index output directories exist
mkdir -p $(dirname $GENOME_INDEX)
mkdir -p $(dirname $NCRNA_INDEX)

echo "Starting HISAT3N index building at $(date)"
start_time=$(date +%s)

# Build genome index
echo "Building HISAT3N genome index..."
$HISAT3N_BUILD --base-change C,T -p 32 $GENOME_FA $GENOME_INDEX

# Build ncRNA index
echo "Building HISAT3N ncRNA index..."
$HISAT3N_BUILD --base-change C,T -p 32 $NCRNA_FA $NCRNA_INDEX

# Record end time
end_time=$(date +%s)
echo "HISAT3N index building completed at $(date)"
echo "Total runtime: $((end_time - start_time)) seconds"

# Display the generated index files
ls -lh ${GENOME_INDEX}*
ls -lh ${NCRNA_INDEX}*

echo "HISAT3N genome and ncRNA indexes built successfully!"