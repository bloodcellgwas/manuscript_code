#!/bin/bash
export CHR=X
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
source /home/pa354/projects/2018_07_06_ukbb500k_2/envars
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment
export MODULEPATH=/rds/project/who1000-1/rds-who1000-cbrc/user/cbrcmod/modules/out/modulefiles:$MODULEPATH 
module load qctool/beta

qctool_v2.0-rc9 \
 -g ${IMPUTE_BGEN_FOLDER}/ukb_imp_chr${CHR}_v3.bgen \
-s /rds/project/who1000-1/rds-who1000-cbrc/user/wja24/shared/parsa_only/ukbb_pheno/ukbb_500k_final_adjusted_chr${CHR}.sample \
-incl-positions ${OUT_DIR}/genfiles/gwsig_snps_chr${CHR}_pos.tsv \
-bgen-bits 8 \
-bgen-compression zlib \
-threads 2 \
-og ${OUT_DIR}/genfiles/chr${CHR}_gwsig.bgen

if [ $? -ne 0 ]; then { echo "qctool bgen failed, aborting." ; exit 1; } fi


qctool_v2.0-rc9 \
-g ${IMPUTE_BGEN_FOLDER}/ukb_imp_chr${CHR}_v3.bgen \
-s /rds/project/who1000-1/rds-who1000-cbrc/user/wja24/shared/parsa_only/ukbb_pheno/ukbb_500k_final_adjusted_chr${CHR}.sample \
-incl-positions ${OUT_DIR}/genfiles/gwsig_snps_chr${CHR}_pos.tsv \
-threads 2 \
-og ${OUT_DIR}/genfiles/chr${CHR}_gwsig.gen

if [ $? -ne 0 ]; then { echo "qctool failed, aborting." ; exit 1; } fi

awk '{gsub("^0*", "", $1); $2=$1":"$4"_"$5"_"$6; $3=$1":"$4"_"$5"_"$6; print}' ${OUT_DIR}/genfiles/chr${CHR}_gwsig.gen > ${OUT_DIR}/genfiles/chr${CHR}_gwsig_clean_id.gen

if [ $? -ne 0 ]; then { echo "awk failed, aborting." ; exit 1; } fi

##rm ${OUT_DIR}/genfiles/chr${CHR}_gwsig.gen

if [ $? -ne 0 ]; then { echo "rm failed, aborting." ; exit 1; } fi

awk '{print $2}' ${OUT_DIR}/genfiles/chr${CHR}_gwsig_clean_id.gen > ${OUT_DIR}/genfiles/ids_${CHR}.tsv

if [ $? -ne 0 ]; then { echo "awk failed, aborting." ; exit 1; } fi

./check_genfiles.R ${CHR}
