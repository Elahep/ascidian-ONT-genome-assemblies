# ascidian-ONT-genome-assemblies
This repository contains the code and workflows used to assemble and annotate ascidian genomes from Oxford Nanopore Technologies long-read sequencing data generated using R10.4.1 flow cells. These assemblies were produced for comparative genomic analyses investigating associations between gene-family evolution and geographic range size in ascidians with contrasting range limits.

We assembled genomes for five species in the order Aplousobranchia:

_Aplidium sp_. (Antarctic endemic), _Aplidium coronum_ (New Zealand range-limited), _Aplidium phortax_ (New Zealand broad-range), _Didemnum marineae_ (New Zealand range-limited), and  _Didemnum jucundum_ (trans-Tasman broad-range).

The initial design of this workflow and several command structures were adapted from Dr Meeran Hussain’s [ONT and Illumina genome assembly and annotation workflow](https://github.com/meeranhussain/Genome_assembly_AND_annotation) for *Microctonus aethiopoides* parasitoid wasps. The workflow presented here was substantially modified for ONT-only assembly and annotation of larger ascidian genomes, with parameters and processing decisions evaluated separately for each species.

## Table of contents
- [Step 0: Basecalling and demultiplexing](#step-0-basecalling-and-demultiplexing)
- [Step 1: Read QC and filtering](#step-1-read-qc-and-filtering)
- [Step 2: Genome assembly](#step-2-genome-assembly)
- [Step 3: Haplotig purging and assembly polishing](#step-3-assembly-polishing)
- [Step 4: Contamination assessment and removal](#step-4-contamination-assessment-and-removal)
- [Step 5: Scaffolding and gap filling](#step-5-scaffolding-and-gap-filling)
- [Step 6: Repeat identification and masking](#step-6-repeat-identification-and-masking)
- [Step 7: Structural annotation and gene-model refinement](#step-6-structural-annotation-and-gene-model-refinement)

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

#### Exploring k-mer profiling with Jellyfish and GenomeScope

Canonical 21-mers were counted from the combined per-species ONT reads using `Jellyfish`. The resulting k-mer frequency histograms were examined with the web-based `GenomeScope 2.0` as an exploratory assessment of genome size, heterozygosity and repeat content.

GenomeScope is designed primarily for low-error short-read data. Because the present histograms were generated from ONT reads, sequencing errors may inflate the number of low-frequency k-mers and bias the fitted genome characteristics. The GenomeScope estimates were therefore treated as exploratory and were not used to specify assembly parameters, select purge-dups thresholds or evaluate final genome sizes.

```
module load Jellyfish
for i in A_coronum A_phortax_1 D_jucundum A_sp D_marineae_1; do 
zcat Combined_${i}_final.fq.gz | jellyfish count -C -m 21 -s 1000000000 -t 12 -o reads_${i}.jf /dev/stdin
jellyfish histo -t 12 reads_${i}.jf > reads_${i}.histo
done
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
_Aplidium coronum_ - `m=33`; 
_Aplidium phortax_ - `m=39`; 
_Didemnum jucundum_ - `m=25`; 
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
```
The polished assembly is written by Medaka as:

`sup_basecalls/Flye/Acoronum_flye/purge_dups/polishing/medaka1/consensus.fasta`

For consistent downstream naming, this can be renamed to:

`A_coronum.racon1.medaka1.fasta`

We again ran Compleasm on the 1x Racon- and 1x medaka-polished assemblies. We also ran `QUAST` to evaluate the quality of the polished assemblies, including obtaining total assembly length, N50, number of contigs, GC content, etc. These assemblies were carried forward to contamination analyses.

```
ml QUAST

mkdir -p quast_compare_5
quast.py \
  -t 4 \
  --eukaryote \
  --large \
  -o quast_compare_5 \
  Acoronum_1xRacon_medaka.fasta \
  Aphortax_1xRacon_medaka.fasta \
  Djucundum_1xRacon_medaka.fasta \
  Dmarineae_1xRacon_medaka.fasta \
  Asp_1xRacon_medaka.fasta
```

## Step 4: Contamination assessment and removal

To assess contamination in our assemblies, we used `BlobToolKit v.1.1``. We needed three complementary sources of information for this step:

| Input                   | Information provided                           |
| ----------------------- | ---------------------------------------------- |
| Polished assembly FASTA | Contig length and GC content                   |
| ONT read-alignment BAM  | Mean read coverage for each contig             |
| MegaBLAST results       | Taxonomic similarity to sequences in NCBI `nt` |

### 4.1 Mapping ONT reads to the polished assembly

The Chopper-filtered ONT reads were mapped back to the corresponding polished assembly using minimap2. The resulting sorted BAM file was used by BlobTools to calculate mean read coverage for each contig.

The following example uses _Aplidium coronum_:
```
ml SAMtools
ml minimap2

minimap2 -t 10 -ax map-ont Acoronum_1xRacon_medaka.fasta ../PX071_tunicate_ont/sup_basecalls/Chopper_filtered/Filtered_A_coronum.fq | samtools view -b - | samtools sort -@ 8 -o Acoronum.bam -

samtools index -c Acoronum.bam
```
This produced:

`Acoronum.bam`
`Acoronum.bam.csi`

The BAM file was not used by MegaBLAST. It was generated independently as the coverage input for BlobTools.

### 4.2 Splitting the assembly into FASTA chunks

The polished assemblies contained thousands of contigs. Before running MegaBLAST, each assembly was divided into smaller FASTA files containing 50 sequences.

The following Python script was saved as `split_fasta_by_nseq.py`:
```
#!/usr/bin/env python3
import sys

fa = sys.argv[1]
prefix = sys.argv[2]
nseq = int(sys.argv[3])

out = None
count = 0
idx = 0

def open_new():
    global out, idx
    idx += 1
    fn = f"{prefix}_{idx:04d}.fasta"
    return open(fn, "w")

with open(fa) as f:
    for line in f:
        if line.startswith(">"):
            if count % nseq == 0:
                if out:
                    out.close()
                out = open_new()
            count += 1
        if out:
            out.write(line)

if out:
    out.close()
```

### 4.3 Running MegaBLAST on the assembly chunks

We ran MegaBLAST independently on every contig chunk.
```
set -euo pipefail

module purge
module load BLASTDB/2026-01 BLAST/2.16.0-GCC-12.3.0
module load Parallel/20220922 Python

FA="Acoronum_1xRacon_medaka.fasta"
CHUNK_PREFIX="Acoronum_chunk"
NSEQ_PER_CHUNK=50

THREADS_TOTAL=${SLURM_CPUS_PER_TASK}   # 16
THREADS_PER_BLAST=2
JOBS=$(( THREADS_TOTAL / THREADS_PER_BLAST ))  # 8

mkdir -p megablast_chunks megablast_out
cp "$FA" megablast_chunks/
cd megablast_chunks

# Split the assembly into chunks before running MegaBLAST
python ../split_fasta_by_nseq.py "$FA" "$CHUNK_PREFIX" "$NSEQ_PER_CHUNK"

# Run MegaBLAST independently on each assembly chunk
ls ${CHUNK_PREFIX}_*.fasta | parallel -j "$JOBS" --halt soon,fail=1 \
  "blastn -query {} \
     -task megablast \
     -db nt \
     -max_target_seqs 5 \
     -culling_limit 1 \
     -evalue 1e-10 \
     -num_threads $THREADS_PER_BLAST \
     -outfmt '6 qseqid staxids pident length evalue bitscore sscinames sskingdoms stitle' \
     -out ../megablast_out/{/.}.tsv"

cd ..
```

The output format contained:

| Column | Information         |
| -----: | ------------------- |
|      1 | Query contig ID     |
|      2 | NCBI TaxID          |
|      3 | Percentage identity |
|      4 | Alignment length    |
|      5 | E-value             |
|      6 | Bitscore            |
|      7 | Scientific name     |
|      8 | Superkingdom        |
|      9 | Hit description     |

### 4.4 Combining the MegaBLAST results
After all chunk-level MegaBLAST searches had finished, the individual result files were concatenated:

`cat megablast_out/*.tsv > Acoronum_megablast.tsv`

### 4.5 Creating a best-hit table

A separate table containing the highest-bitscore hit for each contig was generated for manual inspection:
```
awk -F'\t' '
{
  q=$1; bits=$6;
  if(!(q in best) || bits > best[q]){ best[q]=bits; line[q]=$0 }
}
END{
  for(q in line) print line[q]
}' Aphortax_megablast.tsv > Aphortax_besthit.tsv
```
The complete MegaBLAST output was retained for BlobTools. The best-hit table was used as an additional aid when manually investigating suspicious contigs.

### 4.6 BlobTools analysis

Meeran's workflow used the newer interactive BlobToolKit available through NeSI. However, this version did not run successfully in my NeSI environment. I therefore used command-line BlobTools v1.1, installed directly from GitHub.

The Conda installation did not work, so BlobTools was invoked using the Python executable from its Conda environment.

NCBI taxonomy files were first downloaded and used to construct the BlobTools taxonomy database:
```
mkdir taxdump
cd taxdump

wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz
~/.conda/envs/blobtools/bin/python ./blobtools nodesdb \
  --nodes data/nodes.dmp \
  --names data/names.dmp
```

The following example uses _Aplidium coronum_:
```
~/.conda/envs/blobtools/bin/python ./blobtools create \
  -i ../fq_purged_polished/Acoronum_1xRacon_medaka.fasta \
  -b ../fq_purged_polished/Acoronum.bam \
  -t ../fq_purged_polished/megablast_out_Acoronum/Acoronum_megablast.tsv \
  --nodes ./data/nodes.dmp \
  --names ./data/names.dmp \
  -o Acoronum_blob
```

The FASTA, BAM and MegaBLAST files represented the same polished assembly stage.

Blob plots and tabular summaries were then generated:
```
~/.conda/envs/blobtools/bin/python ./blobtools plot \
  -i Acoronum_blob.blobDB.json \
  -o Acoronum_blob
~/.conda/envs/blobtools/bin/python ./blobtools view \
  -i Acoronum_blob.blobDB.json \
  -o Acoronum_blob
```
The plots and tables were used to examine contig taxonomy, GC content and ONT read coverage. In the BlobTools table, `bam0_mean` represents the mean mapped-read depth for each contig.

### 4.7 Conservative contamination removal

Candidate contaminants were assessed using several lines of evidence:
- Taxonomic assignment
- Position relative to the main chordate GC–coverage distribution
- Mean ONT read coverage
- Contig length
- MegaBLAST percentage identity and alignment length
- MegaBLAST hit description
- Effect of removal on Compleasm completeness

An initial broad strategy removed zero-coverage contigs and most sequences assigned outside Chordata. This substantially reduced Compleasm completeness in several assemblies, indicating that many non-chordate best hits reflected limited taxonomic representation of ascidians rather than genuine contamination.

We therefore rejected the broad filtering strategy. Final removal was restricted to clear microbial contaminants and obvious assembly artefacts supported by the combined BlobTools, coverage and MegaBLAST evidence. Ambiguous metazoan and no-hit contigs were retained unless additional evidence supported their removal.

Species-specific removal lists were generated from the BlobTools results and manually reviewed. The selected contigs were removed using SeqKit:

```
ml SeqKit
seqkit grep -v -f Acoronum.removebacterial.ids ../fq_purged_polished/Acoronum_1xRacon_medaka.fasta > Acoronum_nobact.fasta
```
Sequences shorter than 200 bp were then removed:
```
seqkit seq -m 200 Acoronum_nobact.fasta > Acoronum_nobact_min200.fasta
```
### 4.8 Validation of the cleaned assemblies

The number of removed contigs and the sequence counts before and after filtering were recorded:

`wc -l Acoronum.removebacterial.ids`
`grep -c "^>" Acoronum_1xRacon_medaka.fasta`
`grep -c "^>" Acoronum_nobact_min200.fasta`

We also reran Compleasm and QUAST to confirm that contamination removal had not substantially reduced genome completeness, changed total assembly size or adversely affected contiguity. The cleaned assemblies were then carried forward to scaffolding.

## Step 5: Scaffolding and gap filling
The contamination-filtered assemblies were scaffolded with ntLink using the corresponding Chopper-filtered ONT reads. Different combinations of minimizer k-mer size (k) and window size (w) were tested for each species, as recommended by the ntLink developers:

```
k = 24, 32 or 40
w = 100, 250 or 500
z = 1000
```
Each combination was run using the following code:
```
module load LongStitch
ntLink scaffold target=<cleaned_assembly.fasta> reads=<filtered_ONT_reads.fq> prefix=k${k}_w${w} k=${k} w=${w} z=1000 t=8
```
Candidate assemblies were compared using the number of scaffolds, average and maximum scaffold lengths, total assembly size, number of Ns, QUAST statistics and Compleasm completeness. The selected parameters were:
| Species             | `k` | `w` |
| ------------------- | --: | --: |
| *Aplidium coronum*  |  40 | 500 |
| *Aplidium phortax*  |  24 | 500 |
| *Aplidium* sp.      |  40 | 500 |
| *Didemnum jucundum* |  40 | 500 |
| *Didemnum marinae*  |  32 | 500 |

The selected scaffolded assembly was then subjected to ntLink gap filling. For example:
```
module load LongStitch
ntLink scaffold gap_fill target=Acoronum_nobact_min200.k40.w500.fa reads=Filtered_A_coronum.fq t=8
```
This produced:

`Acoronum_nobact_min200.k40.w500.fa.k32.w100.z1000.ntLink.scaffolds.gap_fill.fa`

Here, k40.w500 identifies the selected input scaffold, while k32.w100.z1000 records the default parameters used during the subsequent gap-filling run. The gap-filled assemblies were carried forward to a final round of Medaka polishing, followed by another Compleasm and QUAST quality assessment.

## Step 6: Repeat identification and masking

[EDTA](https://github.com/oushujun/EDTA) was used to identify and annotate repetitive elements and produce soft-masked genomes for gene prediction. The same workflow was applied to our assemblies and the publicly available *Aplidium turbinatum* and *Didemnum vexillum* genomes to represent two additional wide-range species for each genus.

### 6.1 Preparing genome FASTA files

Before EDTA, contig names were shortened and standardized using species-specific identifiers. This was done for both our assemblies and the NCBI genomes to avoid problems with long or complex FASTA headers.

For the NCBI genomes, a mapping file was retained to connect the renamed sequences to their original identifiers. For example, the *D. vexillum* assembly was renamed as follows:

```bash
awk '
BEGIN {OFS="\t"}
/^>/ {
n++
old=substr($0,2)
split(old,a," ")
new=sprintf("Dv_seq%03d",n)
print new,a[1],old >> "Dv_sequence_name_map.tsv"
print ">"new
next
}
{print}
' GCA_965643705.1_kaDidVexi2_genomic.fna > Dvexillum_renamed.fasta
```

This produced the renamed genome `Dvexillum_renamed.fasta` and the identifier mapping file `Dv_sequence_name_map.tsv`. The same procedure, using appropriate species-specific prefixes, was applied to *A. turbinatum* and our assembled genomes.

For the downloaded NCBI assemblies, we also checked whether the sequence had already been soft-masked:

```bash
grep -v "^>" Dvexillum_renamed.fasta | tr -cd 'a-z' | wc -c
```

A result of `0` indicates that no lowercase, soft-masked sequence is present. Where necessary, lowercase bases were converted to uppercase before EDTA:

```bash
awk '/^>/ {print; next} {print toupper($0)}' Dvexillum_renamed.fasta > Dvexillum_renamed_unmasked.fasta
```

### 6.2 Running EDTA

The following example shows the EDTA analysis for *Aplidium coronum*:

```bash
ml EDTA
EDTA.pl --genome /nesi/nobackup/uow04282/EDTA/Acoronum/Acoronum_ntLink_gapfilling_1xmedaka_renamed.fasta --species others --threads 12 --sensitive 1 --anno 1 --overwrite 0
```

`--sensitive 1` enables sensitive repeat discovery, and `--anno 1` generates genome-wide transposable-element annotations.

The EDTA run creates many files. The sum file contains the information about TE statistics that we need to report in the paper.

The genome is now, by default, hard-masked. We want it to be soft-masked, so we need to do the following.

This is an alternative to RepeatMasker. This will soft-mask the genome, keeping minlen 1000 bp so the genes and near genes won't be masked, as if that happens can interfere with annotation.

### 6.3 Generating the soft-masked genome

First, git clone the EDTA GitHub repository and chmod +x for all the scripts in bin:
git clone https://github.com/oushujun/EDTA.git
chmod -R +x EDTA.pl bin/*.pl

The EDTA annotation was then used to generate a soft-masked genome:

```bash
perl /nesi/project/uow04282/EDTA/EDTA/bin/make_masked.pl \
  -genome Acoronum_ntLink_gapfilling_1xmedaka_renamed.fasta \
  -minlen 1000 \
  -hardmask 0 \
  -t 8 \
  -rmout Acoronum_ntLink_gapfilling_1xmedaka_renamed.fasta.mod.EDTA.anno/Acoronum_ntLink_gapfilling_1xmedaka_renamed.fasta.mod.EDTA.TEanno.out
```

Setting `-hardmask 0` represents annotated repeats as lowercase bases rather than replacing them with `N`s. The resulting soft-masked genomes were used in the subsequent gene-annotation workflow.

## Step 7: Structural annotation and gene-model refinement

Structural annotation was performed using two complementary methods. [GALBA](https://github.com/Gaius-Augustus/GALBA) used protein-to-genome alignments to train AUGUSTUS and predict protein-coding genes. [GeMoMa](https://www.jstacs.de/index.php/GeMoMa) transferred conserved gene structures from annotated ascidian reference genomes. The two prediction sets were integrated using [EVidenceModeler](https://github.com/EVidenceModeler/EVidenceModeler), followed by gene-model rescue and repeat filtering with [QuickProt](https://github.com/thecgs/quickprot).

The soft-masked genome generated by EDTA was used for all analyses. The mitochondrial contig was removed from the genome and annotation files before EVM integration.

### 7.1 Reference annotation data

For GeMoMa, published genome assemblies and transcriptome-supported annotations were retrieved for representatives of all three ascidian orders:

| Order           | Reference species                                                   |
| --------------- | ------------------------------------------------------------------- |
| Aplousobranchia | *Aplidium turbinatum*, *Clavelina lepadiformis*, *Diplosoma virens* |
| Phlebobranchia  | *Ascidiella aspersa*, *Ascidia mentula*, *Phallusia fumigata*       |
| Stolidobranchia | *Styela clava*, *Polycarpa mamillaris*, *Botrylloides leachii*      |

Genome FASTA and GFF files were obtained from public resources, including NCBI and Tunome.

For GALBA, the protein database contained Metazoa proteins from OrthoDB v12 supplemented with available ascidian proteins. Duplicate proteins and sequences shorter than 50 amino acids were removed.

### 7.2 GALBA

GALBA was run through an Apptainer container using a writable AUGUSTUS configuration directory:

```bash
module purge
ml Apptainer/1.3.1

export APPTAINER_CACHEDIR="/nesi/nobackup/uow04282/apptainer-cache"
export APPTAINER_TMPDIR=${APPTAINER_CACHEDIR}
mkdir -p $APPTAINER_CACHEDIR

export AUGUSTUS_CONFIG_PATH=/nesi/nobackup/uow04282/Augustus_config/config

apptainer exec \
  --bind /nesi/nobackup/uow04282:/nesi/nobackup/uow04282 \
  --env AUGUSTUS_CONFIG_PATH=/nesi/nobackup/uow04282/Augustus_config/config \
  /nesi/project/uoo04689/Softwares/galba.sif \
  galba.pl \
  --species=Djucundum_gal \
  --genome=Djucundum_gnm_mskd.fasta \
  --prot_seq=protein_DB.faa \
  --gff3 \
  --AUGUSTUS_ab_initio \
  --threads=24
```

GALBA aligned the supplied proteins to the genome and used the resulting evidence to train AUGUSTUS and predict gene structures.

### 7.3 GeMoMa

GeMoMa v1.9 was run in a dedicated Conda environment with Java 11:

```bash
module purge
ml Miniconda3/23.10.0-1
ml Java/11.0.4

source activate /nesi/project/uow04282/software/gemoma-1.9

export JAVA_HOME=$EBROOTJAVA
export PATH=$JAVA_HOME/bin:$PATH
```

The target genome and reference FASTA and GFF files were linked into the working directory before running GeMoMa:

```bash
GeMoMa GeMoMaPipeline -Xms30G -Xmx60G \
  t=Djucundum_gnm_mskd.fasta \
  s=own i=bput a=Atur_gnm.gff g=Atur_gnm.fasta \
  s=own i=bris a=Clep_gnm.gff g=Clep_gnm.fasta \
  s=own i=isca a=Dvir_gnm.gff g=Dvir_gnm.fasta \
  s=own i=asim a=Aasp_gnm.gff g=Aasp_gnm.fasta \
  s=own i=bras a=Amen_gnm.gff g=Amen_gnm.fasta \
  s=own i=dmel a=Pfum_gnm.gff g=Pfum_gnm.fasta \
  s=own i=iele a=Scla_gnm.gff g=Scla_gnm.fasta \
  s=own i=rmai a=Pmam_gnm.gff g=Pmam_gnm.fasta \
  s=own i=ztil a=Bleachii_gnm.gff g=Bleachii_gnm.fasta \
  outdir=Djucundum_out \
  threads=8 \
  AnnotationFinalizer.r=NO \
  GeMoMa.Score=ReAlign \
  GeMoMa.m=57000
```

GeMoMa used amino-acid and conserved intron-position information from the reference annotations to predict homologous genes in the target genome.

### 7.4 Integrating GALBA and GeMoMa predictions with EVM

The GALBA and GeMoMa GFF3 files were placed in an `01_evidence` directory and renamed to `GALBA.evm.gff3` and `GeMoMa.evm.gff3`.

Auxiliary GeMoMa `GAF` entries were removed before EVM:

```bash
awk '$0 ~ /^#/ || $2!="GAF"' GeMoMa.evm.gff3 > GeMoMa.cleaned.evm.gff3

mv GeMoMa.evm.gff3 GeMoMa.evm_old.gff3
```

EVM combined the two prediction sets using a weight of `1` for the GALBA/AUGUSTUS predictions and `2` for the GeMoMa predictions:

```bash
module purge
ml Perl/5.34.1-GCC-11.3.0

SPECIES="Djucundum"
MASKED_GENOME_FILE="Djucundum_ntLink_gapfilling_1xmedaka_renamed.fasta.new.masked.nomito"
ANNOTATION_DIR="02_evm_annotation"
BRAKER_WEIGHT=1
GEMOMA_WEIGHT=2
EVM_THREADS=22

EVM_DIR="/nesi/project/uow04282/EVM/EVidenceModeler-v2.1.0"
EVM_BIN="${EVM_DIR}/EVidenceModeler"

export PATH="${EVM_DIR}:$PATH"
export PERL5LIB="${EVM_DIR}/PerlLib:${PERL5LIB:-}"

WORKDIR="${ANNOTATION_DIR}/GAL${BRAKER_WEIGHT}_G${GEMOMA_WEIGHT}"
EVIDENCE_DIR="01_evidence"

mkdir -p "${WORKDIR}"

cat "${EVIDENCE_DIR}"/*evm.gff3 > "${WORKDIR}/combined.gff3"

cd "${WORKDIR}"

cat > weights.txt <<EOF
ABINITIO_PREDICTION	AUGUSTUS	${BRAKER_WEIGHT}
OTHER_PREDICTION	GeMoMa	${GEMOMA_WEIGHT}
EOF

perl "${EVM_BIN}" \
  --sample_id "${SPECIES}" \
  --genome "../../${MASKED_GENOME_FILE}" \
  --weights weights.txt \
  --gene_predictions "combined.gff3" \
  --segmentSize 1000000 \
  --overlapSize 100000 \
  --CPU "${EVM_THREADS}"

sed '/^$/d' "${SPECIES}.EVM.gff3" > "${SPECIES}.EVM.mod.gff3"
```

The variable name `BRAKER_WEIGHT` was retained from the original EVM script, but in this analysis it represented the weight assigned to the GALBA/AUGUSTUS predictions.

### 7.5 Gene-model improvement with QuickProt

QuickProt was used to add conserved gene models that were detected by Compleasm/miniprot but absent from the EVM annotation. The resulting annotation was sorted, assigned consistent gene identifiers, and filtered to remove repeat-associated protein predictions.

The `miniprot_output.gff` file was obtained by running Compleasm in genome mode on the same mitochondrial-free soft-masked genome used by EVM.

```bash
ml Python/3.11.3-gimkl-2022a
/nesi/project/uow04282/compleasm_kit/compleasm.py run \
  -a ./Djucundum_renamed_unmasked.fasta.new.masked \
  -o Djucundum_compleasm \
  -l metazoa \
  --odb odb12 \
  -t 12
```

```bash
GENOME=Djucundum_ntLink_gapfilling_1xmedaka_renamed.fasta.new.masked.nomito
EVM=Djucundum.EVM.mod.gff3
MINIPROT=Djucundum_miniprot_output.gff

module purge

unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
unset PERL_MB_OPT
unset PERL_MM_OPT

ml miniprot/0.11-GCC-11.3.0 \
   TransDecoder/5.7.1-GCC-11.3.0-Perl-5.34.1 \
   Python/3.11.6-foss-2023a
```

Missing conserved gene models were added to the EVM annotation:

```bash
python3 /nesi/project/uow04282/software/quickprot/script/update_gff3_from_minibusco.py \
  -r "$EVM" \
  -m "$MINIPROT" \
  -g "$GENOME" \
  -o improve_busco.gff3

cat improve_busco.gff3 "$EVM" > genome.longest.gff.tmp

python3 /nesi/project/uow04282/software/quickprot/script/sort_gff3.py \
  genome.longest.gff.tmp \
| python3 /nesi/project/uow04282/software/quickprot/script/rename_gff3.py \
  - \
  -o genome.longest.gff3 \
  -p Djucundum
```

Protein and CDS sequences were extracted from the improved annotation:

```bash
perl /nesi/project/uow04282/software/quickprot/bin/TransDecoder-5.7.1/util/gff3_file_to_proteins.pl \
  --gff3 genome.longest.gff3 \
  --fasta "$GENOME" \
  --seqType prot \
  > Djucundum_longest.pep.fasta

perl /nesi/project/uow04282/software/quickprot/bin/TransDecoder-5.7.1/util/gff3_file_to_proteins.pl \
  --gff3 genome.longest.gff3 \
  --fasta "$GENOME" \
  --seqType CDS \
  > Djucundum_longest.cds.fasta
```

Finally, repeat-associated predictions were removed:

```bash
module purge

ml DIAMOND/2.1.14-GCC-12.3.0 \
   Python/3.11.6-foss-2023a

mkdir -p 01_filter_repeat
cd 01_filter_repeat

python3 /nesi/project/uow04282/software/quickprot/script/filter_repeatPeps_from_gff3.py \
  -q ../Djucundum_longest.pep.fasta \
  -g ../genome.longest.gff3

sed '/^$/d' retain.gff3 > Djucundum_repfil.gff3
```

The final repeat-filtered peptide set was extracted from the retained gene models:

```bash
module purge

unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
unset PERL_MB_OPT
unset PERL_MM_OPT

ml TransDecoder/5.7.1-GCC-11.3.0-Perl-5.34.1

perl /nesi/project/uow04282/software/quickprot/bin/TransDecoder-5.7.1/util/gff3_file_to_proteins.pl \
  --gff3 Djucundum_repfil.gff3 \
  --fasta "$GENOME" \
  --seqType prot \
  > Djucundum_repfil_pep.fasta
```

The principal outputs used in downstream analyses were:

* `Djucundum_repfil.gff3`: final repeat-filtered structural annotation.
* `Djucundum_repfil_pep.fasta`: corresponding predicted protein sequences.

