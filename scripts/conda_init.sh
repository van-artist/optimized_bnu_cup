#!/bin/bash

WORKSPACE_DIR="$PWD"
MINICONDA_DIR="$WORKSPACE_DIR/miniconda"
ENV_NAME="bio"
YAML_FILE="$WORKSPACE_DIR/bio.yml"

if command -v conda &>/dev/null; then
    echo "=> Conda installation detected."
    CONDA_INSTALLED=true
    CONDA_BASE_DIR="$(conda info --base)"
else
    echo "=> Conda not detected, downloading and installing Miniconda to workspace..."
    curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -o "$WORKSPACE_DIR/Miniconda3-latest-Linux.sh"
    bash "$WORKSPACE_DIR/Miniconda3-latest-Linux.sh" -b -p "$MINICONDA_DIR"
    echo "=> Miniconda installation completed."
    export PATH="$MINICONDA_DIR/bin:$PATH"
    export CONDA_INSTALLED=false
    CONDA_BASE_DIR="$MINICONDA_DIR"
fi

CONDA_EXE="$CONDA_BASE_DIR/bin/conda"

ENV_EXISTS=$($CONDA_EXE info --envs | awk '{print $1}' | grep -x "$ENV_NAME" || echo "")

if [ -n "$ENV_EXISTS" ]; then
    echo "=> Conda environment '$ENV_NAME' already exists, using it directly."
else
    echo "=> Creating Conda environment '$ENV_NAME'..."
    if [ -f "$YAML_FILE" ]; then
        $CONDA_EXE env create -f "$YAML_FILE"
        echo "=> Environment '$ENV_NAME' created with packages from the YAML file."
    else
        echo "=> Environment configuration file '$YAML_FILE' not found, please check the path."
        exit 1
    fi
fi

echo "=> Activating Conda environment '$ENV_NAME'..."
source "$CONDA_BASE_DIR/bin/activate" "$ENV_NAME"

if [[ "$CONDA_DEFAULT_ENV" == "$ENV_NAME" ]]; then
    echo "=> Conda environment '$ENV_NAME' activated successfully!"
else
    echo "=> Failed to activate Conda environment '$ENV_NAME', please check the error message."
    exit 1
fi

echo "=> Environment setup completed!"