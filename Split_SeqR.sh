#!/bin/bash

#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : bash                                                                       |
#  Study       : WES PID                                                                    |
#  Data Owner  : Newcastle University - Prof. Sophie Hambleton                              |
#  Description : Split VCF file using GATK's SelectVariants function, based on input ped.   |
#-------------------------------------------------------------------------------------------#


##'Set variables
##'-----------------------------------------------------------------------------------------#
VEP_DIR="/home/andrew/.vep/"
VCF_IN="/home/andrew/2016Oct_SeqR_WES/data_prepped/Hambleton_Jan2017.vep.vcf"
VCF_OUT="/home/andrew/2016Oct_SeqR_WES/data_prepped/Hambleton_Jan2017_NG052.vep.vcf"
WRK_DIR="/home/andrew/2016Oct_SeqR_WES/data_prepped/"
VEP="/home/andrew/Software/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl"
FA_IN="/home/andrew/.vep/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly2.fa"
PED="/home/andrew/2016Oct_SeqR_WES/data_prepped/Samples.ped"
JAR_IN="/opt/databases/GenomeAnalysisTK-3.4-protected/GenomeAnalysisTK.jar"
BUNDLE="/opt/databases/hg19/"
##'-----------------------------------------------------------------------------------------#


##'Split based on Pedigree file
##'-----------------------------------------------------------------------------------------#
# grep "FAMD" ${PED} > ${WRK_DIR}/tmp.ped
list=`cat "${WRK_DIR}/Temp.ped" | cut -f 2 | sed 's/^/--sample_name /' -`
echo ${list}
java -Xmx4g -jar \
    ${JAR_IN} \
      -T SelectVariants \
      -R ${BUNDLE}/ucsc.hg19.fasta \
      --downsampling_type NONE \
      ${list} \
      --variant ${VCF_IN} \
      --out ${VCF_OUT}
# rm ${WRK_DIR}/tmp.ped
##'-----------------------------------------------------------------------------------------#


##'Zip and tabix index
##'-----------------------------------------------------------------------------------------#
bgzip ${VCF_OUT}
tabix -p vcf ${VCF_OUT}.gz
##'-----------------------------------------------------------------------------------------#
