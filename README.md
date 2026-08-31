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
Raw Oxford Nanopore signal data were basecalled and barcode-classified by Annabel Whibley (Bragato Research Institute) using Dorado v1.1.1 with the SUP model. Reads were demultiplexed using the sequencing sample sheets and converted from BAM to compressed FASTQ format. The downstream assembly workflow presented in this repository began with the resulting per-species FASTQ files.

# Dorado SUP basecalling and barcode classification
```
module purge
module load Dorado/1.1.1

POD5_DIR="/path/to/pod5_directory"
SAMPLE_SHEET="/path/to/sample_sheet.csv"
RUN_ID="example_run"

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
# Demultiplexing
```
module purge
module load SAMtools/1.22-GCC-12.3.0
module load pigz/2.7

DEMUX_DIR="example_run.sup_demux"

for bam in "${DEMUX_DIR}"/*.bam; do
    [[ "${bam}" == *_unclassified.bam ]] && continue

    fq="${bam%.bam}.fq.gz"

    samtools fastq -@ 16 "${bam}" |
        pigz -p 16 -9 > "${fq}"
done
```
# Combining per-species raw FASTQ
```
cat run1_sample.fq.gz run2_sample.fq.gz > species_raw_reads.fq.gz
```

## Step 1: Read QC and filtering
