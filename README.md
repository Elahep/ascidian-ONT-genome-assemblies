# ascidian-ONT-genome-assemblies
This repository contains the code and workflows used to assemble and annotate ascidian genomes from Oxford Nanopore Technologies long-read sequencing data generated using R10.4.1 flow cells. These assemblies were produced for comparative genomic analyses investigating associations between gene-family evolution and geographic range size in ascidians with contrasting range limits.

We assembled genomes for five species in the order Aplousobranchia:

_Aplidium sp_. (Antarctic endemic), _Aplidium coronum_ (New Zealand range-limited), _Aplidium phortax_ (New Zealand broad-range), _Didemnum marineae_ (New Zealand range-limited), and  _Didemnum jucundum_ (trans-Tasman broad-range).

## Table of contents
- [Step 0: Basecalling and demultiplexing](#step-0-basecalling-and-demultiplexing)
- [Step 1: Read QC and filtering](#step-1-read-qc-and-filtering)
- [Step 2: Genome assembly](#step-2-genome-assembly)
- [Step 3: Haplotig purging and assembly polishing](#step-3-assembly-polishing)
- [Step 4: Contamination assessment](#step-3-contamination-assessment)

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
We examined read quality and length distribution for each species with Nanoplot. Based on the Nanoplot results, we filtered raw reads with a minimum quality score of 9 (Q9) and a minimum read length of 1,000 bp. For _Aplidium phortax_ we used a minimum read length of 500 bp to retain more reads in this sample which had shorter read-length distribution. 
We ran Nanoplot again after filtering to evaluate the retained read-length and quality distributions and to record summary statistics, including total read count, total yield and read N50.

```
module load NanoPlot/1.43.0-foss-2023a-Python-3.11.6
export PYTHONNOUSERSITE=1

NanoPlot --fastq ./Chopper_filtered/Filtered_D_marineae_1.fq -t 16 -o marinae_chopper_nanoplot

gunzip Combined_D_marineae_1_final.fq.gz
gunzip Combined_A_coronum_final.fq.gz
gunzip Combined_A_phortax_1_final.fq.gz
gunzip Combined_D_jucundum_final.fq.gz
gunzip Combined_A_sp_final.fq.gz

module load chopper
for i in D_jucundum A_coronum A_sp D_marineae_1
do
chopper --threads 12 -q 9 -l 1000 < Combined_${i}_final.fq > Filtered_${i}.fq ;
done
chopper --threads 8 -q 9 -l 500 < Combined_A_phortax_1_final.fq > Filtered_Aphortax.fastq
```

## Step 2: Genome assembly
We used FLYE for genome assembly using the filtered fq files from the previous step for each species separately. FLYE requires an estimate of the genome size (the `-g` parameter). We set `-g` as 600 Mb for the _Aplidium_ species and 1.1 Gb for the _Didemnum_ species.

```
module load Flye
flye --nano-hq ../Chopper_filtered/Filtered_A_coronum.fq --genome-size 600m --asm-coverage 40 --out-dir Sadareanum_flye --threads 20
flye --nano-hq ../Chopper_filtered/Filtered_A_phortax.fq --genome-size 600m --asm-coverage 40 --out-dir Sadareanum_flye --threads 20
flye --nano-hq ../Chopper_filtered/Filtered_A_sp.fq --genome-size 600m --asm-coverage 40 --out-dir Sadareanum_flye --threads 20
flye --nano-hq ../Chopper_filtered/Filtered_D_jucundum.fq --genome-size 1.1g --asm-coverage 40 --out-dir Sadareanum_flye --threads 20
flye --nano-hq ../Chopper_filtered/Filtered_D_marineae_1.fq --genome-size 1.1g --asm-coverage 40 --out-dir Sadareanum_flye --threads 20
```

## Step 3: Haplotig purging and assembly polishing
We first ran `Compleasm` with the `metazoa_odb10` lineage on the FLYE-assembled genomes to obtain a baseline genome completeness for each assembly. Assemblies containing more than 2% duplicated BUSCOs were processed with `purge_dups`. We aimed to reduce duplication below 10% while, where possible, retaining missing BUSCOs below approximately 10% (keeping D acceptably low (<10%) without a big penalty in M).

We  initially estimated purge-dups coverage thresholds from the mapped-read depth distribution. Compleasm was run again after purging, and the thresholds were adjusted where necessary to balance the removal of duplicated regions against the loss of genome completeness.

The Flye assembly of _Aplidium sp_. contained less than 2% duplicated BUSCOs and, therefore, we did not process it with purge_dups.

Following haplotig purging, the retained primary assemblies were polished using one round of `Racon` followed by one round of `Medaka`. The same polishing procedure was applied directly to the Flye assembly of _Aplidium sp_.

### 3.1 Baseline completeness assessment

The following example uses _Aplidium coronum_:
```
module purge
module load compleasm/0.2.5-gimkl-2022a

SPECIES="A_coronum"
FLYE_ASSEMBLY="sup_basecalls/Flye/Acoronum_flye/assembly.fasta"

compleasm.py run -a "${FLYE_ASSEMBLY}" -o "compleasm/${SPECIES}_flye" -l metazoa_odb10 -t 12
```
Assemblies with duplicated BUSCOs above 2% were carried forward to haplotig purging.

### 3.2 Mapping reads for purge_dups

Chopper-filtered ONT reads were mapped back to the corresponding Flye assembly. The alignments were converted to PAF format for read-depth estimation with `pbcstat`.

```
module purge
module load minimap2/2.28-GCC-12.3.0
module load SAMtools/1.22-GCC-12.3.0

SPECIES="A_coronum"
THREADS=14

ASM="sup_basecalls/Flye/Acoronum_flye/assembly.fasta"
READS="reads/${SPECIES}.q9.min1kb.fq.gz"
OUTDIR="sup_basecalls/Flye/Acoronum_flye/purge_dups"

mkdir -p "${OUTDIR}"

minimap2 -a -x map-ont -t "${THREADS}" "${ASM}" "${READS}" | samtools sort -@ "${THREADS}" -o "${OUTDIR}/reads_vs_assembly.bam"

samtools index "${OUTDIR}/reads_vs_assembly.bam"

samtools view -@ "${THREADS}" -h "${OUTDIR}/reads_vs_assembly.bam" > "${OUTDIR}/reads_vs_assembly.sam"
```

### 3.3 Haplotig identification and removal
```
module purge
module load purge_dups/1.2.6-gimkl-2022a-Python-3.10.5
module load minimap2/2.28-GCC-12.3.0

SPECIES="A_coronum"
ASM="sup_basecalls/Flye/Acoronum_flye/assembly.fasta"
OUTDIR="sup_basecalls/Flye/Acoronum_flye/purge_dups"

cd "${OUTDIR}"

# Convert mapped-read alignments to PAF
paftools.js sam2paf reads_vs_assembly.sam > "${SPECIES}.reads_vs_assembly.paf"

# Calculate read-depth statistics
pbcstat "${SPECIES}.reads_vs_assembly.paf"

# Estimate the initial depth thresholds
calcuts PB.stat > cutoffs 2> calcuts.log

# Split the assembly and perform a self-alignment
split_fa "${ASM}" > assembly.split.fasta

minimap2 -x asm5 -DP assembly.split.fasta assembly.split.fasta | gzip -c > assembly.self.paf.gz

# Identify duplicated regions
purge_dups -2 -T cutoffs -c PB.base.cov assembly.self.paf.gz > dups.bed 2> purge_dups.log

# Generate the retained primary assembly and removed haplotigs
get_seqs -e dups.bed "${ASM}"
```

This produces:
```
purged.fa    retained primary assembly
hap.fa       removed haplotigs
```
### 3.4 Evaluating and adjusting the purge-dups thresholds

Compleasm was rerun on purged.fa:
```
module purge
module load compleasm/0.2.5-gimkl-2022a

SPECIES="A_coronum"
PURGED_ASSEMBLY="sup_basecalls/Flye/Acoronum_flye/purge_dups/purged.fa"

compleasm.py run -a "${PURGED_ASSEMBLY}" -o "compleasm/${SPECIES}_purged" -l metazoa_odb10 -t 12
```
If excessive duplication remained or missing BUSCOs increased substantially, the read-depth histogram was inspected and the purge-dups thresholds were adjusted:
```
calcuts -l LOWER_CUTOFF -m HAPLOID_DIPLOID_CUTOFF -u UPPER_CUTOFF PB.stat > cutoffs.manual
```

The most important parameter was `-m`, which defines the depth transition between haploid/primary and diploid/duplicated sequence. We repeated purging and Compleasm assessment until an acceptable balance between duplication and completeness was achieved.

The selected cutoff file was then used to produce the retained purged.fa. We used the following  `-m` cutoff for different species:
_Aplidium coronum_ - `m=33'
_Aplidium phortax_ - `m=39`
_Didemnum jucundum_ - `m=25`
_Didemnum marineae_ = `m=32`

### 3.5 One round of Racon polishing

For species processed with purge_dups, purged.fa was used as the Racon input. For _Aplidium sp._, the original Flye assembly was used instead.
```
module load minimap2
module load SAMtools
module load Racon

ASM="sup_basecalls/Flye/Acoronum_flye/purge_dups/purged.fa"
READS="sup_basecalls/Chopper_filtered/Filtered_A_coronum.fq"
OUTDIR="sup_basecalls/Flye/Acoronum_flye/purge_dups/polishing"

mkdir -p "${OUTDIR}"

minimap2 -a -x map-ont -t "${THREADS}" "${ASM}" "${READS}" | samtools sort -@ "${THREADS}" -o "${OUTDIR}/reads_vs_purged.bam"

samtools view -h "${OUTDIR}/reads_vs_purged.bam" > "${OUTDIR}/reads_vs_purged.sam"

racon -t "${THREADS}" "${READS}" "${OUTDIR}/reads_vs_purged.sam" "${ASM}" > "${OUTDIR}/${SPECIES}.racon1.fasta"
```
### 3.6 One round of Medaka polishing

The Racon-polished assembly was subsequently polished with Medaka v2.2.0:
```
module purge
module load medaka/2.2.0-Miniforge3-25.3.1-0

SPECIES="A_coronum"
THREADS=6

READS="sup_basecalls/Chopper_filtered/Filtered_A_coronum.fq"
RACON_ASSEMBLY="sup_basecalls/Flye/Acoronum_flye/purge_dups/polishing/${SPECIES}.racon1.fasta"
MEDAKA_OUTDIR="sup_basecalls/Flye/Acoronum_flye/purge_dups/polishing/medaka1"

medaka_consensus -i "${READS}" -d "${RACON_ASSEMBLY}" -o "${MEDAKA_OUTDIR}" -t "${THREADS}"

The polished assembly is written by Medaka as:

`sup_basecalls/Flye/Acoronum_flye/purge_dups/polishing/medaka1/consensus.fasta`

For consistent downstream naming, this can be renamed to:

`A_coronum.racon1.medaka1.fasta`

We again ran Compleasm on the 1x Racon and 1x medaka polished assembly. These assemblies were carried forward to contamination analyses.

