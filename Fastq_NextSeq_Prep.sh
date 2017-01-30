#!/bin/bash
#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : Bash                                                                       |
#  Study       : Exome Project                                                              |
#  Data Owner  : Newcastle University - Prof. Sophie Hambleton                              |
#  Type        : Sulaco Script                                                              |
#  Description : Prep raw data from Illumina NextSeq, and concatenate gunzipped fastq files |
#                to unify the four lane split. Current version takes directly from          |
#                a mounted Illumina BaseSpace project directory. Exclusively designed to    |
#                deal with paired end reads                                                 |
#  Version     : 1.2                                                                        |
#  Input       : BaseSpace Project Directory (hard coded variable - PROJ_BASE)              |
#-------------------------------------------------------------------------------------------#


# Set Input directory, and output directory
PROJ_BASE="/home/andrew/BaseSpace/Projects/2016_016/Samples/"
OUTDIR="/home/andrew/Raw_WES/2016/June/"

# Create the output directory
mkdir -p $OUTDIR

# Loop through each subdirectory of the main Project
# directory, in essence, loop through each sample
for i in ${PROJ_BASE}*
do
    # Get Sample ID
    SAMPLE_ID=$(basename "$i")
    echo "New Sample ${SAMPLE_ID}"

    # Create the sample's directory in output
    # data structure, and create empty, gunzipped
    # fastq files
    mkdir -p ${OUTDIR}Sample_${SAMPLE_ID}/Raw_Data/
    touch ${OUTDIR}Sample_${SAMPLE_ID}/Raw_Data/${SAMPLE_ID}_R1.fastq.gz
    touch ${OUTDIR}Sample_${SAMPLE_ID}/Raw_Data/${SAMPLE_ID}_R2.fastq.gz

    # Iterate over the 4 fastq files and write the contents
    # to the empty fastq files...
    echo "Forward Reads"
    for x in ${i}/Files/*_R1_*
    do
        echo $x
        cat ${x} >> ${OUTDIR}Sample_${SAMPLE_ID}/Raw_Data/${SAMPLE_ID}_R1.fastq.gz
    done

    # ...Repeate for the reverse reads
    echo "Reverse Reads"
    for y in ${i}/Files/*_R2_*
    do
        echo $y
        cat ${y} >> ${OUTDIR}Sample_${SAMPLE_ID}/Raw_Data/${SAMPLE_ID}_R2.fastq.gz
    done
done
