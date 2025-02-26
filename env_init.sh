#!/bin/bash

WORKSPACE_DIR="$(pwd)/workspace"
PACKAGES_DIR="$(pwd)/packages"
ROOT_DIR="$(pwd)"
SCRIPTS_DIR="./scripts"
DATA_DIR="./data"
REW_DIR="$DATA_DIR/raw"
REF_DIR="$DATA_DIR/ref"

MINICONDA_DIR="./miniconda"

HISAT_3N_DIR="$PACKAGES_DIR/hisat3n"


mkdir -p "$WORKSPACE_DIR"
mkdir -p "$MINICONDA_DIR"
mkdir -p "$PACKAGES_DIR"
mkdir -p "$DATA_DIR" "$REW_DIR" 

chmod +x "$SCRIPTS_DIR/conda_init.sh"

bash "$SCRIPTS_DIR/conda_init.sh"

cd "$HISAT_3N_DIR"
make





