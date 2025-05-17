
# Bioinformatics Module

Automated, genome-resolved metagenomics preprocessing pipeline for urban microbiome biosurveillance.

This folder contains SLURM-compatible shell scripts and Perl utilities for executing end-to-end preprocessing, quality control, MAG recovery, and output aggregation from SRA-derived shotgun metagenomic datasets.

---

## Pipeline Overview

This module enables clean, high-quality input for genome assembly, binning, annotation, and downstream analysis. It is optimized for high-throughput studies across hospital sewage, ambulances, public transit, and clinical environments.

---

## Scripts and Functionalities

### `01_run.sh` — Automated Preprocessing & SRA Integration

#### Overview:
- Automates data acquisition, QC, and preprocessing from NCBI SRA.
- Optimized for high-contact environmental metagenomes.

#### Key Steps:
- Downloads SRA data using `pysradb`
- Converts `.sra` to paired-end FASTQ (`fasterq-dump`)
- Quality assessment via `FastQC`
- Adapter and quality trimming using **BBTools**:
  - `reformat.sh` — ensure correct PE formatting
  - `bbmerge.sh` — merge overlapping reads
  - `clumpify.sh` — deduplicate
  - `bbduk.sh` — trim and filter adapters

---

### `02_fastq_screen.sh` — MAG Recovery & Functional Annotation

#### Overview:
- Performs MAG assembly and profiling from cleaned reads.

#### Key Steps:
- Assembly: `metaSPAdes`
- Binning: `MetaBAT2` with depth calculation
- MAG QC: `CheckM` (≥90% completeness, ≤10% contamination)
- Annotation:
  - `DIAMOND` BLASTX vs **VFDB** (virulence)
  - `DIAMOND` BLASTX vs **CARD** (AMR genes)
- Taxonomy: `Kraken2`, `Bracken`
- Outputs: Merged, structured annotation files per sample

---

### `03_collect_files.sh` — Output Aggregation & Cleanup

#### Overview:
- Collects all final outputs across samples and organizes by:
  - FastQ Screen logs
  - MAG bins
  - QUAST reports
  - Taxonomic profiles (species/genus/family)
  - VFDB and CARD annotations

#### Output Structure:
```
collect_AMR_Env/<project_name>/
├── FQS_Before/
├── FQS_After/
├── Bin/
├── bracken_species/
├── VFDB/
└── CARD/
```

---

### `simplifyFastaHeaders.pl` — Header Standardization

- Cleans and shortens FASTA headers for tool compatibility
- Use post-assembly/pre-annotation

---

##  Dependencies (Conda Envs Expected)

- BBMap suite (bbduk, bbmerge, reformat, clumpify)
- FastQC, FastQ Screen
- metaSPAdes, MetaBAT2, CheckM
- Kraken2, Bracken
- DIAMOND
- Perl 5.x
- pysradb, fasterq-dump

---

## Example Usage

```bash
bash run.sh PRJNA123456
bash fastq_screen.sh
perl simplifyFastaHeaders.pl contigs.fasta > renamed_contigs.fasta
bash collect_files.sh hospital_project
```

---

## Scientific Context

Supports:
- WHO AMR GLASS objectives
- One Health biosurveillance
- Predictive modeling of ecological risk

This module prepares microbial data for use in R/Python environments or visualization software.

---

## Notes

- Designed for HPC with SLURM
- Temporary files automatically cleaned
- Final outputs are sample-labeled and stored safely
