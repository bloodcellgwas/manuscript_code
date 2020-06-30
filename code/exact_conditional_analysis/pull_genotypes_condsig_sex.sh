#!/bin/bash

export CHR=23
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
source /home/pa354/projects/2018_07_06_ukbb500k_2/envars
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment
export MODULEPATH=/rds/project/who1000-1/rds-who1000-cbrc/user/cbrcmod/modules/out/modulefiles:$MODULEPATH 
module load qctool/beta

#-g ${OUT_DIR}/genfiles/chr${CHR}_gwsig_clean_id.gen \
#-s ${SAMPLE_FILE} \
qctool_v2.0-rc9 \
-g ${OUT_DIR}/genfiles/chrXY_merged.gen \
-s ${OUT_DIR}/genfiles/chrXY_merged.sample \
-incl-snpids ${OUT_DIR}/condout/pull_ids/pull_ids_chr${CHR}.tsv \
-og ${OUT_DIR}/condout/subgenhd5/chr${CHR}.gen

if [ $? -ne 0 ]; then { echo "qctool failed, aborting." ; exit 1; } fi

awk '{gsub("^0*","", $1); print $1":"$4"_"$5"_"$6}' ${OUT_DIR}/condout/subgenhd5/chr${CHR}.gen > ${OUT_DIR}/condout/subgenhd5/chr_ids${CHR}.tsv

if [ $? -ne 0 ]; then { echo "awk failed, aborting." ; exit 1; } fi

Rscript ./pull_genotypes_condsig_check.R ${CHR}

if [ $? -ne 0 ]; then { echo "extract CHECK failed, aborting." ; exit 1; } fi

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_2/gen2hd5_lib/

/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_2/gen2hd5.git/src/gen2hd5 \
-g ${OUT_DIR}/condout/subgenhd5/chr${CHR}.gen \
-o ${OUT_DIR}/condout/subgenhd5/chr${CHR}.h5

if [ $? -ne 0 ]; then { echo "hd5 failed, aborting." ; exit 1; } fi

echo completed
