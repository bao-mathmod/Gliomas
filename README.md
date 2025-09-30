# Gliomas Pipeline

Reproducible pipelines (Nextflow) for downloading, QC, and processing sc/sn data.

## Quickstart
```bash
# create env, for example
# conda env create -f env.yml && conda activate glio

nextflow run main.nf -profile slurm
