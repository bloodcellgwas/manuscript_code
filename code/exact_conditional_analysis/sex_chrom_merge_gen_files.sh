#!/bin/bash
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
source /home/pa354/projects/2018_10_11_ukbb500k_2/envars
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment
export MODULEPATH=/rds/project/who1000-1/rds-who1000-cbrc/user/cbrcmod/modules/out/modulefiles:$MODULEPATH 
module load qctool/beta

## subset the first file
qctool_v2.0-rc9 \
    -g ${OUT_DIR}/genfiles/chrX_gwsig_clean_id.gen \
    -s /rds/project/who1000-1/rds-who1000-cbrc/user/wja24/shared/parsa_only/ukbb_pheno/ukbb_500k_final_adjusted_chrX.sample \
    -og ${OUT_DIR}/genfiles/chrX_gwsig_clean_id_sample_filt.gen \
    -threads 10 \
    -incl-samples ${OUT_DIR}/sex_chrom_x.sample
if [ $? -ne 0 ]; then { echo "first failed, aborting." ; exit 1; } fi 

### subset the second file
qctool_v2.0-rc9 \
    -g ${OUT_DIR}/genfiles/chrXY_gwsig_clean_id.gen \
    -s /rds/project/who1000-1/rds-who1000-cbrc/user/wja24/shared/parsa_only/ukbb_pheno/ukbb_500k_final_adjusted_chrXY.sample \
    -og ${OUT_DIR}/genfiles/chrXY_gwsig_clean_id_sample_filt.gen \
    -threads 10 \
    -incl-samples ${OUT_DIR}/sex_chrom_x.sample
if [ $? -ne 0 ]; then { echo "second failed, aborting." ; exit 1; } fi

## merge both files
qctool_v2.0-rc9 \
    -g ${OUT_DIR}/genfiles/chrX_gwsig_clean_id_sample_filt.gen \
    -s ${OUT_DIR}/sex_chrom_x.sample \
    -merge-in ${OUT_DIR}/genfiles/chrXY_gwsig_clean_id_sample_filt.gen ${OUT_DIR}/sex_chrom_x.sample \
    -threads 10 \
    -og ${OUT_DIR}/genfiles/chrXY_merged.gen \
    -os ${OUT_DIR}/genfiles/chrXY_merged.sample

if [ $? -ne 0 ]; then { echo "third failed, aborting." ; exit 1; } fi

awk '{print $2}' ${OUT_DIR}/genfiles/chrXY_merged.gen > ${OUT_DIR}/genfiles/ids_chrXY_merged.tsv
if [ $? -ne 0 ]; then { echo "awk failed, aborting." ; exit 1; } fi


