# ascidian-ONT-genome-assemblies
This repository contains the code and workflows used to assemble and annotate ascidian genomes from Oxford Nanopore Technologies long-read sequencing data generated using R10.4.1 flow cells. These assemblies were produced for comparative genomic analyses investigating associations between gene-family evolution and geographic range size in ascidians with contrasting range limits.

We assembled genomes for five species in the order Aplousobranchia:

Aplidium sp. (Antarctic endemic), Aplidium coronum (New Zealand range-limited), Aplidium phortax (New Zealand broad-range), Didemnum marineae (New Zealand range-limited), and  Didemnum jucundum (trans-Tasman broad-range).

## Table of contents
- [Step 0: Basecalling and demultiplexing](#step-0-basecalling-and-demultiplexing)
- [Step 1: Read QC and filtering](#step-1-read-qc-and-filtering)
- [Step 2: Genome assembly](#step-2-genome-assembly)
- [Step 3: Assembly polishing](#step-3-assembly-polishing)


## Step 0: SUP basecalling and demultiplexing
Raw Oxford Nanopore signal data were basecalled and barcode-classified by Dr Annabel Whibley (Bragato Research Institute) using Dorado v1.1.1 with the SUP model. Reads were demultiplexed using the sequencing sample sheets and converted from BAM to compressed FASTQ format. The downstream assembly workflow presented in this repository began with the resulting per-species FASTQ files.

### Dorado SUP basecalling and barcode classification

#### Sample-sheet structure

The following is an example from one of the sample sheets used during basecalling. The `alias` column associates each barcode with its corresponding species. For example, reads assigned to `barcode26` belong to *Didemnum jucundum* and are represented by the alias `D_jucundum`.

| protocol_run_id | position_id | flow_cell_id | sample_id | experiment_id | flow_cell_product_code | kit | barcode | alias | type |
|---|---|---|---|---|---|---|---|---|---|
| 92a57ab4-3f57-44d0-b955-da7af5f12ef3 | 2F | PBC89758 | LX0146 | R0521 | FLO-PRO114M | SQK-NBD114-96 | barcode26 | D_jucundum | test_sample |
| 92a57ab4-3f57-44d0-b955-da7af5f12ef3 | 2F | PBC89758 | LX0146 | R0521 | FLO-PRO114M | SQK-NBD114-96 | barcode28 | A_coronum | test_sample |

The `dorado basecaller` command produces one major BAM file for the entire sequencing run. That BAM contains basecalled reads from every barcode/sample included in that run. In this project, we used three PromethION flow cells and 13 sequencing runs across those flow cells. We basecalled each run independently, producing one combined BAM containing all barcode-classified reads from that run. The combined BAMs were subsequently demultiplexed into sample-specific files.  
```
module purge
module load Dorado/1.1.1

POD5_DIR="/path/to/pod5_directory"
SAMPLE_SHEET="/path/to/sample_sheet.csv"

# PX071: project identifier
# R0521: experiment_id recorded in the sample sheet
OUTPUT_PREFIX="PX071_R0521"

# SUP basecalling and barcode classification
dorado basecaller sup \
    --kit-name SQK-NBD114-96 \
    --sample-sheet "${SAMPLE_SHEET}" \
    --device cuda:all \
    --recursive "${POD5_DIR}" \
    > "${RUN_ID}.sup.bam"

# Separate reads using the barcode assignments generated during basecalling
dorado demux \
    --output-dir "${RUN_ID}.sup_demux" \
    --no-classify \
    "${RUN_ID}.sup.bam"
```
### Demultiplexing
```
module purge
module load SAMtools/1.22-GCC-12.3.0
module load pigz/2.7

DEMUX_DIR="example_run.sup_demux"
SAMPLE_ALIAS="D_jucundum"

INPUT_BAM="${DEMUX_DIR}/${SAMPLE_ALIAS}.bam"
OUTPUT_FASTQ="${SAMPLE_ALIAS}.raw.fq.gz"

samtools fastq -@ 16 "${INPUT_BAM}" |
    pigz -p 16 -9 > "${OUTPUT_FASTQ}"
```
### Combining per-species raw FASTQ
```
cat run1_sample.fq.gz run2_sample.fq.gz > species_raw_reads.fq.gz
```

## Step 1: Read QC and filtering
