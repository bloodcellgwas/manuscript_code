#!/bin/bash
export OUT_DIR=/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_4/sysmex

### HPC ###
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment
export MODULEPATH=/rds/project/who1000-1/rds-who1000-cbrc/user/cbrcmod/modules/out/modulefiles:$MODULEPATH 
module load qctool/beta


# use QCTOOL to pull the table of conditionally significant SNPS in BGEN format
##$SUPPORT_FILES/qctool \
##-g ${OUT_DIR}/condout/ldclump_compare/chr#_gwsig.gen \
##qctool_v2.0-rc9 \
##-g ${IMPUTE_BGEN_FOLDER}/ukb_imp_chr#_v3.bgen \
##-s ${SAMPLE_FILE} \
##-incl-positions ${OUT_DIR}/condout/ldclump_compare/pull_table.tsv \
##-og ${OUT_DIR}/condout/ldclump_compare/all_oldid.gen

##if [ $? -ne 0 ]; then { echo "qctool failed, aborting." ; exit 1; } fi

awk '{$2=$1":"$4"_"$5"_"$6;$3=""; print $0}' ${OUT_DIR}/condout/ldclump_compare/all_oldid.gen > ${OUT_DIR}/condout/ldclump_compare/all.gen

# create r2 table
##/nfs/team151/software/plink2_18_April_2015/plink \
module load plink/1.90beta
plink \
    --gen ${OUT_DIR}/condout/ldclump_compare/all.gen \
    --sample ${SAMPLE_FILE} \
    --ld-window-r2 0 \
    --inter-chr \
    --r2

if [ $? -ne 0 ]; then { echo "awk failed, aborting." ; exit 1; } fi

# move the results to the relevant folder
mv plink.nosex ${OUT_DIR}/condout/ldclump_compare/plink.nosex.ld
mv plink.log ${OUT_DIR}/condout/ldclump_compare/plink.log.ld
mv plink.ld ${OUT_DIR}/condout/ldclump_compare

# find out how many of the clumps are novel (don't include any SNPs in the previous list)
${OUT_DIR}/novel_clumps_new.R

################ OLD METHOD #########################
# perform the clumping
# creates assoc file which is a list of SNPs to input to plink for clumping
#${OUT_DIR}/create_clump_compare_assoc.R
#/nfs/team151/software/plink2_18_April_2015/plink \
#    --gen ${OUT_DIR}/condout/ldclump_compare/all.gen \
#    --sample ${SAMPLE_FILE} \
#    --clump-r2 0.6 \
#    --clump-p1 0.01 \
#    --clump-kb 250000 \
#    --clump ${OUT_DIR}/condout/ldclump_compare/all.assoc \
# move the results to the relevant folder
#mv plink.nosex ${OUT_DIR}/condout/ldclump_compare
#mv plink.log ${OUT_DIR}/condout/ldclump_compare
#mv plink.clumped ${OUT_DIR}/condout/ldclump_compare
# format the plink output, set environment variable so this is saved in the right directory
#export COMPARE="true"
#${OUT_DIR}/format_plink_clumps.R
#####################################################
