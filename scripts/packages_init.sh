#!/bin/bash

if [ -z "$1" ]; then
    echo "=> Please provide the target directory path as an argument."
    exit 1
fi

TARGET_DIR="$1"
JAVA_URL="https://download.oracle.com/java/17/archive/jdk-17.0.12_linux-x64_bin.tar.gz"
HISAT3N_URL="https://api.bigcousin1.cn/resources/hisat-3n.zip"
M5C_UBSSEQ_URL="https://api.bigcousin1.cn/resources/m5C-UBSseq.zip"
UMICOLLAPSE_URL="https://api.bigcousin1.cn/resources/umicollapse.jar"
SNAPPY_URL="https://api.bigcousin1.cn/resources/snappy-java-1.1.9.1.jar"
HTSJDK_URL="https://api.bigcousin1.cn/resources/htsjdk-4.1.3-9-gc31bc92-SNAPSHOT.jar"

JAVA_DIR="$TARGET_DIR/java"
HISAT3N_DIR="$TARGET_DIR/hisat3n"
M5C_UBSSEQ_DIR="$TARGET_DIR/m5C-UBSseq"

# Ensure the target directory exists
mkdir -p "$TARGET_DIR"

# Check if necessary tools are installed
for cmd in curl unzip make xz; do
    if ! command -v $cmd &>/dev/null; then
        echo "=> Error: $cmd is missing in the system, cannot proceed."
        exit 1
    fi
done

# Download and extract Java 17 JDK
echo "=> Downloading Java 17 JDK to $JAVA_DIR ..."
mkdir -p "$JAVA_DIR"
curl -fsSL "$JAVA_URL" -o "$JAVA_DIR/jdk-17.tar.gz"
if [ $? -ne 0 ]; then
    echo "=> Java JDK download failed, exiting script."
    exit 1
fi

tar -xzf "$JAVA_DIR/jdk-17.tar.gz" -C "$JAVA_DIR"
rm "$JAVA_DIR/jdk-17.tar.gz"
echo "=> Java 17 JDK installation completed!"

# Set JAVA_HOME environment variable
export JAVA_HOME="$JAVA_DIR/jdk-17"
export PATH="$JAVA_HOME/bin:$PATH"
echo "=> Set JAVA_HOME=$JAVA_HOME"

# Download and extract HISAT3N package
echo "=> Downloading HISAT3N package to $HISAT3N_DIR ..."
mkdir -p "$HISAT3N_DIR"
curl -fsSL "$HISAT3N_URL" -o "$TARGET_DIR/hisat3n.zip"
if [ $? -ne 0 ]; then
    echo "=> HISAT3N download failed, exiting script."
    exit 1
fi

# Extract to target directory using unzip
echo "=> Extracting HISAT3N package ..."
unzip -q "$TARGET_DIR/hisat3n.zip" -d "$HISAT3N_DIR"
rm "$TARGET_DIR/hisat3n.zip"

# Check and handle nested directories
if [ -d "$HISAT3N_DIR/hisat-3n" ]; then
    mv "$HISAT3N_DIR/hisat-3n"/* "$HISAT3N_DIR/"
    rmdir "$HISAT3N_DIR/hisat-3n" # Remove empty nested directory
fi

# Ensure cd succeeds, otherwise exit
cd "$HISAT3N_DIR" || {
    echo "=> Failed to enter HISAT3N directory, exiting script."
    exit 1
}

# Execute compilation
echo "=> Compiling HISAT3N ..."
make || {
    echo "=> HISAT3N compilation failed, exiting script."
    exit 1
}

# Download and extract m5C-UBSseq package
echo "=> Downloading m5C-UBSseq package to $M5C_UBSSEQ_DIR ..."
mkdir -p "$M5C_UBSSEQ_DIR"
curl -fsSL "$M5C_UBSSEQ_URL" -o "$M5C_UBSSEQ_DIR/m5C-UBSseq.zip"
if [ $? -ne 0 ]; then
    echo "=> m5C-UBSseq download failed, exiting script."
    exit 1
fi

# Extract m5C-UBSseq package
echo "=> Extracting m5C-UBSseq package ..."
unzip -q "$M5C_UBSSEQ_DIR/m5C-UBSseq.zip" -d "$M5C_UBSSEQ_DIR"
rm "$M5C_UBSSEQ_DIR/m5C-UBSseq.zip"

echo "=> m5C-UBSseq installation completed!"