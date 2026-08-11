# nf-fastp

`nf-fastp` trims and filters paired-end FASTQ files with `fastp`.

This module is normally the second step in the ChIP-seq pipeline:

```text
nf-fastqc -> nf-fastp -> nf-bwa
```

`nf-fastqc` checks raw FASTQ quality. `nf-fastp` then cleans the reads before alignment.

## What It Does

1. Finds paired-end FASTQ files.
2. Skips samples that already have all expected fastp outputs.
3. Runs `fastp` for each sample pair.
4. Saves trimmed FASTQ files and fastp QC reports.

## Before You Run

You need:

- Nextflow installed
- Access to paired-end FASTQ files
- For local runs: Docker available
- For HPC runs: Slurm and Singularity available
- For HPC runs: the fastp Singularity image path in `configs/slurm.config` must exist

Default HPC notification email:

```text
molendo.hpc@gmail.com
```

## Input Options

There are two ways to provide FASTQ pairs.

### Option 1: Samples Master CSV Recommended

Use this when running the full ChIP-seq pipeline or when processing multiple samples.

The CSV should contain these columns:

```csv
sample_id,fastq_r1,fastq_r2,enabled
WT_rep1,/path/to/WT_rep1_R1.fastq.gz,/path/to/WT_rep1_R2.fastq.gz,true
WT_rep2,/path/to/WT_rep2_R1.fastq.gz,/path/to/WT_rep2_R2.fastq.gz,true
```

Notes:

- `sample_id` becomes the output filename prefix.
- `fastq_r1` and `fastq_r2` are required.
- This module expects paired-end data.
- Rows with `enabled=false` are skipped.
- Empty `enabled` values are treated as enabled.

Run with:

```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --project_folder /path/to/output_project \
  --fastp_output fastp_output
```

### Option 2: FASTQ Folder + Pair Pattern

Use this for quick module testing.

```bash
nextflow run main.nf -profile hpc \
  --fastp_raw_data /path/to/fastq_folder \
  --fastp_pattern "*_R{1,2}_001.fastq.gz" \
  --project_folder /path/to/output_project \
  --fastp_output fastp_output
```

The default pattern expects filenames like:

```text
sample_R1_001.fastq.gz
sample_R2_001.fastq.gz
```

If neither `--samples_master` nor `--fastp_raw_data` is provided, the module stops with an error.

## Output

Results are written to:

```text
${project_folder}/${fastp_output}/
```

Example:

```text
/path/to/output_project/fastp_output/
  WT_rep1_R1.fastp.trimmed.fastq.gz
  WT_rep1_R2.fastp.trimmed.fastq.gz
  WT_rep1.fastp.html
  WT_rep1.fastp.json
```

Output files:

- `*_R1.fastp.trimmed.fastq.gz`: trimmed read 1 FASTQ
- `*_R2.fastp.trimmed.fastq.gz`: trimmed read 2 FASTQ
- `*.fastp.html`: human-readable fastp report
- `*.fastp.json`: machine-readable report, useful for MultiQC or parsing

## Recommended HPC Run

From inside the `nf-fastp` folder:

```bash
cd /path/to/nf-fastp

nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --project_folder /path/to/output_project \
  --fastp_output fastp_output
```

Resume a previous run:

```bash
nextflow run main.nf -profile hpc -resume \
  --samples_master /path/to/samples_master.csv \
  --project_folder /path/to/output_project \
  --fastp_output fastp_output
```

Override the HPC notification email:

```bash
nextflow run main.nf -profile hpc \
  --samples_master /path/to/samples_master.csv \
  --project_folder /path/to/output_project \
  --mail_user molendo.hpc@gmail.com
```

## Local Test Run

Use local mode only for small test data:

```bash
nextflow run main.nf -profile local \
  --fastp_raw_data /path/to/test_fastq \
  --fastp_pattern "*_R{1,2}_001.fastq.gz" \
  --project_folder ./test_output
```

## Key Parameters

| Parameter | Meaning | Default |
| --- | --- | --- |
| `samples_master` | CSV containing paired FASTQ paths | `null` |
| `fastp_raw_data` | FASTQ input folder when not using CSV | `null` |
| `fastp_pattern` | Paired FASTQ filename pattern | `*_R{1,2}_001.fastq.gz` |
| `project_folder` | Base output folder | current folder |
| `fastp_output` | fastp output subfolder | `fastp_output` |
| `cpus` | CPUs per fastp task | `4` |
| `memory` | Memory per task | `8GB` |
| `time` | Runtime limit per task | `2h` |
| `mail_user` | HPC notification email | `molendo.hpc@gmail.com` |

## Existing Results Are Skipped

For each sample, the module checks whether all four expected output files already exist:

```text
sample_R1.fastp.trimmed.fastq.gz
sample_R2.fastp.trimmed.fastq.gz
sample.fastp.html
sample.fastp.json
```

If all four exist, that sample is skipped. If any file is missing, fastp runs again for that sample.

## How To Check Results

Open the HTML report:

```text
${project_folder}/${fastp_output}/sample.fastp.html
```

Useful sections to inspect:

- before filtering quality
- after filtering quality
- adapter trimming
- read length distribution
- total reads retained

The trimmed FASTQ files are the inputs for the next module, usually `nf-bwa`.

## Troubleshooting

If the run fails:

1. Check the main Nextflow log:

```bash
less .nextflow.log
```

2. Check the failed task error file:

```bash
less work/<hash>/.command.err
```

3. Common problems:

- `fastq_r1` or `fastq_r2` path is wrong in `samples_master`.
- The data are single-end, but this module expects paired-end reads.
- The filename pattern does not match the FASTQ names.
- `configs/slurm.config` points to a missing Singularity image.
- The HPC bind path in `extra_mounts` does not include the FASTQ or output location.
- Docker is not running for local mode.

## Project Structure

```text
main.nf
nextflow.config
configs/
  local.config
  slurm.config
```
