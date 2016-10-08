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
* bamm        (v1.5.0)
* bedtools    (v2.20.1)
* dirseq      (v0.0.2)
* extern      (v0.1.0)
* fastqc      (v0.10.1)
* fxtract     (v1.2)
* graphviz    (v2.38.0)
* numpy       (v1.9.1)
* python      (v2.7.4)
* ruffus      (v2.6.3)
* samtools    (v0.1.19)
* sortmerna   (v2.0)
* tempdir     (v0.6)
* trimmomatic (v0.32)

## Databases
In order to remove contaminant sequences, TranscriptM requires 3 databases containing: 

1. The sequences of adapters using during the sequencing 
2. The sequence of the PhiX genome (used as control in Illumina sequencing)
3. The sequences of ribosomal, transfer and transfer-messenger RNA  

## Usage
```sh
$ transcriptm -h

usage: transcriptm [-h] [--verbose [VERBOSE]] [--version] [-L FILE] [-n]
                   [--flowchart FILE] [--draw_graph_horizontally]
                   [--flowchart_format FORMAT] --paired_end PAIRED_END
                   [PAIRED_END ...] [--metaG_contigs METAG_CONTIGS]
                   [--dir_bins DIR_BINS] [--threads THREADS]
                   [--halt_after_stage HALT_AFTER_STAGE]
                   [--restart_from_stage RESTART_FROM_STAGE]
                   [--db_path DB_PATH] [--output_dir OUTPUT_DIR]
                   [--working_dir WORKING_DIR] [--adapters {nextera,truseq}]
                   [--min_len MIN_LEN] [--min_avg_qc MIN_AVG_QC]
                   [--phred {phred33,phred64}] [--min_qc MIN_QC] [--crop CROP]
                   [--headcrop HEADCROP] [--path_db_smr PATH_DB_SMR]
                   [--percentage_id PERCENTAGE_ID]
                   [--percentage_aln PERCENTAGE_ALN] [--no_mapping_filter]

transcriptm v0.4.0: Metatranscriptomic data processing and complete
metagenomics analysis. Example usage: transcriptm --paired_end sample1-1.fq.gz
sample1-2.fq.gz sample2-1.fq.gz sample2-2.fq.gz --metaG_contigs assembly.fa
--dir_bins dir_gff

General options:
  -h, --help            show this help message and exit
  --paired_end PAIRED_END [PAIRED_END ...]
                        required input files: paired sequences files of raw
                        metatranscriptomics reads (.fq.gz format) e.g.
                        --paired_end sample1_1.fq.gz sample1_2.fq.gz
                        sample2_1.fq.gz sample2_2.fq.gz
  --metaG_contigs METAG_CONTIGS
                        All contigs from the reference metagenome in a fasta
                        file
  --dir_bins DIR_BINS   Directory which contains several annotated population
                        genomes (bins) -> gff format, the others files would
                        be ignored
  --threads THREADS     Number of threads to use
  --halt_after_stage HALT_AFTER_STAGE
                        An intermediate stage. After this stage is completed,
                        the pipeline run will halt. Available options:
                        [view_raw_reads, trim_raw_reads, phiX_removal,
                        filter_rna, prep_for_mapping, map_to_reference,
                        bam_stats, summary_tables]. Important: The pipeline
                        may later be restarted from the next stage (or
                        previous stages), but note: (i) the same
                        `--paired_end` arguments must be supplied in exactly
                        the SAME ORDER as before, else mismatching will occur;
                        and (ii) the `--working_dir` argument must be
                        specified explicitly and the data must be retained in
                        the specified directory between runs.
  --restart_from_stage RESTART_FROM_STAGE
                        An intermediate stage to restart from, using
                        previously saved data in the `--working_dir`
                        directory. Available options are the same as per
                        option `--halt_after_stage`. NB. See important notes
                        in the help for option `--halt_after_stage`.
  --db_path DB_PATH     Directory which contains the TranscriptM databases
  --output_dir OUTPUT_DIR
                        Output directory
  --working_dir WORKING_DIR
                        Working directory (which will be created if it does
                        not exist). Specifying a working directory is
                        essential if `--halt_after_stage` or
                        `--restart_from_stage` are required in this run or a
                        future restarted run. If a working directory is not
                        specified, a temporary directory is dynamically
                        created and disposed of.

Common options:
  --verbose [VERBOSE], -v [VERBOSE]
                        Print more verbose messages for each additional
                        verbose level.
  --version             show program's version number and exit
  -L FILE, --log_file FILE
                        Name and path of log file

pipeline arguments:
  -n, --just_print      Don't actually run any commands; just print the
                        pipeline.
  --flowchart FILE      Don't run any commands; just print pipeline as a
                        flowchart.
  --draw_graph_horizontally
                        Draw horizontal dependency graph.
  --flowchart_format FORMAT
                        format of dependency graph file. Can be 'pdf', 'svg',
                        'svgz' (Structured Vector Graphics), 'pdf', 'png'
                        'jpg' (bitmap graphics) etc

Trimmomatic options:
  --adapters {nextera,truseq}
                        Type of adapters to clip
  --min_len MIN_LEN     Minimum required length of read
  --min_avg_qc MIN_AVG_QC
                        Minimum average quality score for 4 bp windows
  --phred {phred33,phred64}
                        Quality encoding
  --min_qc MIN_QC       Minimum quality score for leading and trailing bases
  --crop CROP           Cut read to a specific length
  --headcrop HEADCROP   Cut specified number of bases from start of read

SortMeRNA options:
  --path_db_smr PATH_DB_SMR
                        Path to databases and index e.g.
                        path_db1,path_index1:path_db2,path_index2 [default:
                        rRNA and tRNA db] NB: index must be created with the
                        script sortmerna/2.0/bin/indexdb_rna

Mapping options (BamM filter):
  --percentage_id PERCENTAGE_ID
                        Minimum allowable percentage base identity of a mapped
                        read
  --percentage_aln PERCENTAGE_ALN
                        Minimum allowable percentage read bases mapped
  --no_mapping_filter   Do not adjust the mapping stringency by filtering
                        alignments

```
