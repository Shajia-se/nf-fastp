# nf-fastp

`nf-fastp` is a Nextflow DSL2 module for paired-end FASTQ preprocessing using `fastp`.

## fastp Tool Function

`fastp` is an all-in-one FASTQ QC + preprocessing tool. In this module it is mainly used to:
- trim adapters
- trim low-quality bases
- filter low-quality reads
- generate QC reports before/after filtering

Compared with older multi-tool workflows, `fastp` runs these steps in one command and outputs both HTML and JSON reports.

## What This Module Does

1. Reads paired FASTQ files (priority: `samples_master` > `fastp_raw_data + fastp_pattern`).
2. Skips samples that already have `${sample}.fastp.html` in output.
3. Runs `fastp` for each sample pair.
4. Publishes trimmed FASTQ and reports to output directory.

## Input

Preferred mode (`samples_master`):
- `params.samples_master`: CSV with at least `sample_id,fastq_r1,fastq_r2` (optional `enabled`)
- module runs only rows with `enabled=true` (or empty)

Fallback mode (pattern-based):
- Directory: `params.fastp_raw_data`
- Pattern: `params.fastp_pattern` (default: `*_R{1,2}_001.fastq.gz`)
- Expected input: paired-end reads (`R1`/`R2`)

## Output

Under `${project_folder}/${fastp_output}` (default `./fastp_output/`):
- `${sample}_R1.fastp.trimmed.fastq.gz`
- `${sample}_R2.fastp.trimmed.fastq.gz`
- `${sample}.fastp.html`
- `${sample}.fastp.json`

## Key Parameters

- `project_folder`: output base folder (default: `$PWD`)
- `samples_master`: preferred input table
- `fastp_raw_data`: fallback input FASTQ folder
- `fastp_pattern`: paired FASTQ matching pattern
- `fastp_output`: output folder name
- `cpus`, `memory`, `time`: process resources

## Run

```bash
nextflow run main.nf -profile local
```

```bash
nextflow run main.nf -profile hpc
```

Recommended run (`samples_master`):
```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --fastp_output fastp_output
```

Fallback run (pattern):

```bash
nextflow run main.nf -profile hpc \
  --fastp_raw_data /your/raw_fastq \
  --fastp_pattern "*_R{1,2}_001.fastq.gz" \
  --fastp_output fastp_output
```

Resume:

```bash
nextflow run main.nf -profile hpc -resume
```

## How To Check If Trimming Happened

Open `${sample}.fastp.html` and compare before/after sections:
- total reads / total bases
- Q20 / Q30 rates
- adapter trimmed reads/bases
- read length distribution after filtering

For command-line check, inspect `${sample}.fastp.json`.

## Project Structure

```text
main.nf
nextflow.config
configs/
  local.config
  slurm.config
```
