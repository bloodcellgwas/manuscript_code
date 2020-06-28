#!/bin/bash

# use QCTOOL to pull the table of conditionally significant SNPS in BGEN format

#Rscript ./condsig_var_incl.R
#
#export MODULEPATH=/rds/project/who1000-1/rds-who1000-cbrc/user/cbrcmod/modules/out/modulefiles:$MODULEPATH 
#module load qctool/beta
#qctool_v2.0-rc9 \
#-g  ${IMPUTE_BGEN_FOLDER}/ukb_imp_chr#_v3.bgen \
#-s ${SAMPLE_FILE} \
#-incl-positions ${OUT_DIR}/condout/ldclump/pull_table.tsv  \
#-og ${OUT_DIR}/condout/ldclump/condsig_oldid.gen
#
#if [ $? -ne 0 ]; then { echo "qctool failed, aborting." ; exit 1; } fi
#
## plink requires the gen files to be in the following format with third column empty
#awk '{gsub("^0*", "", $1); $2=$1":"$4"_"$5"_"$6; $3=""; print}' ${OUT_DIR}/condout/ldclump/condsig_oldid.gen #> ${OUT_DIR}/condout/ldclump/condsig.gen
#
#if [ $? -ne 0 ]; then { echo "awk failed, aborting." ; exit 1; } fi
#
#awk '{print $2}' ${OUT_DIR}/condout/ldclump/condsig.gen > ${OUT_DIR}/condout/ldclump/ids.tsv
#
#if [ $? -ne 0 ]; then { echo "awk failed, aborting." ; exit 1; } fi
#
## create r2 table
#    --maf 0.00001 \
#    --geno 0.999999999999 \

module load plink/1.90beta
#plink \
#    --gen ${OUT_DIR}/condout/ldclump/condsig.gen \
#    --sample ${SAMPLE_FILE} \
#    --extract ${OUT_DIR}/condout/ldclump/pull_table_ids.tsv \
#    --make-bed \
#    --out ${OUT_DIR}/condout/ldclump/plink_conv
#
#if [ $? -ne 0 ]; then { echo "first plink failed, aborting." ; exit 1; } fi

#plink \
#    --bfile ${OUT_DIR}/condout/ldclump/plink_conv \
#    --ld-snp-list ${OUT_DIR}/condout/ldclump/pull_table_ids.tsv \
#    --ld-window-r2 0 \
#    --ld-snps 3:50035015_T_C \
#    --r2 inter-chr yes-really
plink --bfile ${OUT_DIR}/condout/ldclump/plink_conv --ld 3:50035015_T_C
if [ $? -ne 0 ]; then { echo "plink failed, aborting." ; exit 1; } fi
# move the results to the relevant folder

#### OLD CLUMPING METHOD IGNORE ######
# this method clumps the variants using PLINK, but in our method we do LD
# using PLINK and then a bespoke script to clump them together. 
# PLINK method will probably led to a higher number of clumps
# ${OUT_DIR}/create_clump_assoc.R
#/nfs/team151/software/plink2_18_April_2015/plink \
#    --gen ${OUT_DIR}/condout/ldclump/condsig.gen \
#    --sample ${SAMPLE_FILE} \
#    --clump-r2 0.8 \
#    --clump-p1 0.01 \
#    --clump-kb 250000 \
#    --clump ${OUT_DIR}/condout/ldclump/condsig.assoc \
#####################################
mv plink.nosex ${OUT_DIR}/condout/ldclump
mv plink.log ${OUT_DIR}/condout/ldclump
mv plink.ld ${OUT_DIR}/condout/ldclump
Rscript ./novel_clumps.R

