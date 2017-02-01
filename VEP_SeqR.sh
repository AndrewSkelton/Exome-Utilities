#!/bin/bash

#-------------------------------------------------------------------------------------------#
#  Author      : Andrew J Skelton                                                           |
#  Language    : bash                                                                       |
#  Study       : WES PID                                                                    |
#  Data Owner  : Newcastle University - Prof. Sophie Hambleton                              |
#  Description : Run Ensembl's Variant Effect Predictor on filtered VCF file from Joint     |
#                Calling. Add in specific annotations for dbNSFP (v3.2a), CADD, LoF, used   |
#                by seqr.                                                                   |
#-------------------------------------------------------------------------------------------#


##'Set variables
##'-----------------------------------------------------------------------------------------#
VEP_DIR="/home/andrew/.vep/"
VCF_IN="/home/andrew/2016Oct_SeqR_WES/data_in/VQSR_Recalibrated.vcf"
VCF_OUT="/home/andrew/2016Oct_SeqR_WES/data_prepped/Hambleton_Jan2017.vep.vcf"
WRK_DIR="/home/andrew/2016Oct_SeqR_WES/"
VEP="/home/andrew/Software/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl"
FA_IN="/home/andrew/.vep/homo_sapiens/86_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly2.fa"
##'-----------------------------------------------------------------------------------------#


##'Run VEP
##'-----------------------------------------------------------------------------------------#
perl ${VEP} \
      --everything --vcf --allele_number --no_stats --cache \
      --fork 35 --dir ${VEP_DIR} --assembly GRCh37 --hgvs \
      --offline --force_overwrite --sift b --polyphen b --symbol --numbers \
      --biotype --total_length --canonical --ccds --buffer_size 80000 \
      --plugin dbNSFP,${VEP_DIR}Plugins/dbNSFP/dbNSFPv3_2a.gz,Polyphen2_HVAR_pred,CADD_phred,SIFT_pred,FATHMM_pred,MutationTaster_pred,MetaSVM_pred \
      --plugin LoF,human_ancestor_fa:${VEP_DIR}Plugins/loftee/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15 \
      -i ${VCF_IN} -o ${VCF_OUT}
##'-----------------------------------------------------------------------------------------#


##'Run VEP - recommendation by seqr devs
##'-----------------------------------------------------------------------------------------#
# VEP SeqR
# perl ${VEP} \
#       --everything --vcf --allele_number --no_stats --cache \
#       --offline --force_overwrite --fork 20 --dir ${VEP_DIR} \
#       --fasta ${FA_IN}/ --cache_version 86 --assembly GRCh37 \
#       --tabix --buffer_size 75000 --fork 20 \
#       --plugin dbNSFP,${VEP_DIR}Plugins/dbNSFP/dbNSFPv3_2a.gz,Polyphen2_HVAR_pred,CADD_phred,SIFT_pred,FATHMM_pred,MutationTaster_pred,MetaSVM_pred \
#       --plugin LoF,human_ancestor_fa:${VEP_DIR}Plugins/loftee/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15 \
#       -i ${VCF_IN} -o ${VCF_OUT}
##'-----------------------------------------------------------------------------------------#
