# RNA Methylation Analysis Pipeline

This pipeline is designed for processing RNA methylation sequencing data, including environment initialization, sequence alignment, UMI deduplication, and result analysis.

---

## Directory Structure Requirements

**Before running**, ensure the following directory structure is set up:

## Quick Start

1. Set Execution Permissions

```shell
chmod +x ./env_init.sh ./run.sh ./run_pipline.sh
```

2 . Prepare Input Files

- Place raw sequencing files (`*.fastq`) in the `./data/raw/` directory.
- Place reference genome files in the `./data/ref/` directory:
  - Homo_sapiens.GRCh38.dna.primary_assembly.fa (genome index)
  - Homo_sapiens.GRCh38.ncrna.fa (ncRNA index)

3. Run the Pipeline

```shell
./run.sh
```

The final results will be generated in the `./result/` directory, including:

## Notes

- Screenshots of the runtime results and execution times are stored in the `shots/` directory.
- run_pipline.sh is the workflow description file.
- The sites_file/ directory contains the m5C site files.
- The folder itself is a packaged, fully executable pipeline. Simply follow the instructions provided earlier to run it. The `packages/` directory contains the necessary toolkits. The `conda_init.sh` script in the `scripts/ `directory automatically sets up a virtual Python environment, while `packages_init.sh` can be used to re-download and compile various toolkits in case of file corruption.
