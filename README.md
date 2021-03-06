# transcriptM

## Overview
* **Process several sets of metatranscriptomic paired-end reads (sequenced with Illumina)**
  - Quality trimming
  - PhiX reads removal
  - rRNA, tRNA and tmRNA removal
* **Complete metagenomics analysis**
  - Map processed metatranscriptomic reads against metagenomic contigs
  - Remove reads which mapped with low stringency
  - Compute coverage of annotated genes from several population genomes (i.e. bins, sets of metagenomic contigs that might represent individual genome) in different samples.

```sh
transcriptm --paired_end sample1-R1.fq.gz sample1-R2.fq.gz sample2-R1.fq.gz sample2-R2.fq.gz --metaG_contigs assembly.fa --dir_bins dir_gff
```
## What does it produce ?
* **FastQC_raw**
  - Fastqc reports of raw reads
* **FastQC_processed**
  - Fastqc reports of processed reads
* TranscriptM_output_COUNT.csv  
  - Raw count of mapped reads per gene (provided in gff files)
* TranscriptM_output_NORM_COVERAGE.csv 
  - Average coverage of mapped reads per gene (provided in gff files) normalized by the total of mapped reads 
* **log**
  - Log file of each step
* summary_reads
  - Distribution of reads after each step (IMPORTANT: the unit is paired-end reads)

## Dependencies
* ruffus      (v2.6.3)
* fastqc      (v0.10.1)
* trimmomatic (v0.32)
* bamm        (v1.5.0)
* fxtract     (v1.2)
* sortmerna   (v2.0)
* samtools    (v0.1.19)
* numpy       (v1.9.1)
* dirseq      (v0.0.2)
* bedtools    (v2.20.1)
* tempdir     (v0.6)

## Databases
In order to remove contaminant sequences, TranscriptM requires 3 databases containing: 

1. The sequences of adapters using during the sequencing 
2. The sequence of the PhiX genome (used as control in Illumina sequencing)
3. The sequences of ribosomal, transfer and transfer-messenger RNA  

## Usage
```sh
$ transcriptm -h

optional arguments:
  -h, --help            show this help message and exit
  --paired_end PAIRED_END [PAIRED_END ...]
                        Input files: paired sequences files of raw metatranscriptomic reads (fq.gz format) 
                        e.g. --paired_end sample1_1.fq.gz sample1_2.fq.gz sample2_1.fq.gz sample2_2.fq.gz
  --metaG_contigs METAG_CONTIGS
                        All contigs from the reference metagenome in a fasta file
  --dir_bins DIR_BINS   Directory which contains several annotated population genomes (bins) 
                        -> gff format, the others files would be ignored
  --threads THREADS     Number of threads to use (default=20)
  --db_path DB_PATH     Directory which contains the TranscriptM databases
  --output_dir OUTPUT_DIR
                        Output directory (default: ./TranscriptM_output)
  --working_dir WORKING_DIR
                        Working directory (default: /tmp)
                        
Trimmomatic options:
  --adapters {nextera,truseq}
                        Type of adapters to clip (default='truseq')
  --min_len MIN_LEN     Minimum required length of read (default = 30)
  --min_avg_qc MIN_AVG_QC
                        Minimum average quality score for 4 bp windows (default = 25)
  --phred {phred33,phred64}
                        Quality encoding (default='phred33')
  --min_qc MIN_QC       Minimum quality score for leading and trailing bases (default = 20)
  --crop CROP           Cut read to a specific length (default = 10000)
  --headcrop HEADCROP   Cut specified number of bases from start of read (default=0)
  
  
SortMeRNA options:
  --path_db_smr PATH_DB_SMR
                        Path to databases and index (created with sortmerna/2.0/bin/indexdb_rna) 
                        e.g. path_db1,path_index1:path_db2,path_index2 

Mapping options (BamM filter):
  --percentage_id PERCENTAGE_ID
                        Minimum allowable percentage base identity of a mapped read (default=0.97)
  --percentage_aln PERCENTAGE_ALN
                        Minimum allowable percentage read bases mapped (default=0.95)
  --no_mapping_filter   Do not adjust the mapping srtingency by filtering alignments (default=False)

```
