# Exome Utilities
This repository holds useful scripts to carry out bespoke checks in exome data. Each script is detailed in the below list:

## Fastq_NextSeq_Prep.sh
A bash script to take the four gunzipped fastq files from an Illumina Nextseq and unify these into a single fastq file.

## Gender.R
An R script which takes a `ped` file (based on Plink's format structure), along with coverage files of the SRY gene (+-1Kb). The idea is to test the average coverage of SRY, to infer sex. The inferred sex is tested against the encoded sex in the ped file, those that are suspect mismatches are visualised, and a spreadsheet is output.  

## Module_PedigreeCheck.sh
A module from the PID-WES-GATK3.4-SGE repository, which checks the pedigree structure of the `ped` file, vs the samples that are available for analysis. For example, if sample `x`, has a maternal ID of sample `y`, but sample `y` is not physically present for analysis, the maternal ID for sample `x` will be set to `0`.

## Relatedness.R
This script checks for potential sample swaps and systematic errors. Providing a sample population is big enough (>20 individuals), then a relatedness coefficient can be calculated via one of two methods; Yang et al, Nature Genetics 2010, or Manichaikul et al., BIOINFORMATICS 2010. From a `ped` file of expected sample structures, and the results of a relatedness test, lineages can be derived and tested against known structures. Mismatches and potential swaps are flagged.

## RmDup_Run.sh & RmDup.sh
These scripts are designed to be run on a compute cluster running SoGE, and they take a clean GATK bam file and remove all duplicate marked SAM flags in the alignments (0x0400)
