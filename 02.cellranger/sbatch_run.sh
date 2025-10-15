#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=NF
#SBATCH --output=logs/nf_%j.out
#SBATCH --error=logs/nf_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=LocalQ

nextflow run main.nf \
   -profile slurm,singularity \
   --input /mnt/18T/chibao/gliomas/data/fastq/official/all_projects.csv \
   --cellranger_index /mnt/12T/chibao/env_tool/cellranger/ref_genome/refdata-gex-GRCh38-2024-A \
   --aligner cellranger \
   --outdir /mnt/18T/chibao/gliomas/data/output_cell/nextflow \
   --skip_fastqc \
   --skip_cellbender \
   -resume \
   2>&1 | tee /mnt/18T/chibao/gliomas/data/output_cell/nextflow/nextflow_run.log


nextflow run main.nf \
  -profile slurm,singularity \
  -w /mnt/10T/chibao/gliomas/code/scrnaseq-modified/work \
  --input /mnt/10T/chibao/gliomas/code/scrnaseq-modified/run_cmd/PRJNA869964.csv \
  --cellranger_index /mnt/12T/chibao/env_tool/cellranger/ref_genome/refdata-gex-GRCh38-2020-A \
  --aligner cellranger \
  --outdir /mnt/10T/chibao/gliomas/data/PRJNA869964 \
  --skip_fastqc \
  --skip_cellbender \
  -resume 30d49178-a7c6-4e79-84bc-c9397ce34505 \
  2>&1 | tee /mnt/10T/chibao/gliomas/data/PRJNA869964/nextflow_run.log

nextflow run main.nf \
  -profile slurm,singularity \
  --input /mnt/10T/chibao/gliomas/data/PRJNA961045/PRJNA961045.csv \
  --cellranger_index /mnt/12T/chibao/env_tool/cellranger/ref_genome/refdata-gex-GRCh38-2020-A \
  --aligner cellranger \
  --protocol ARC-v1 \
  --outdir /mnt/10T/chibao/gliomas/data/PRJNA961045 \
  --skip_fastqc \
  -resume \
  2>&1 | tee /mnt/10T/chibao/gliomas/data/PRJNA961045/nextflow_run.log

nextflow run main.nf \
  -profile slurm,singularity \
  --aligner cellrangerarc \
  --input /mnt/10T/chibao/gliomas/data/PRJNA961045/PRJNA961045_arc.csv \
  --cellrangerarc_reference /mnt/12T/chibao/env_tool/cellranger_arc/refgenome_arc/refdata-cellranger-arc-GRCh38-2024-A \
  --outdir /mnt/10T/chibao/gliomas/data/PRJNA961045 \
  -resume \
  2>&1 | tee /mnt/10T/chibao/gliomas/data/PRJNA961045/nextflow_arc_run.log

# Replace your run command with this (note the param name change):
nextflow run main.nf \
  -profile slurm,singularity \
  --aligner cellrangerarc \
  --input /mnt/10T/chibao/gliomas/data/PRJNA961045/PRJNA961045_arc.csv \
  --cellranger_index /mnt/12T/chibao/env_tool/cellranger_arc/refgenome_arc/refdata-cellranger-arc-GRCh38-2024-A \
  --save_align_intermeds true \
  --outdir /mnt/10T/chibao/gliomas/data/PRJNA961045 \
  -resume \
  2>&1 | tee /mnt/10T/chibao/gliomas/data/PRJNA961045/nextflow_arc_run.log
