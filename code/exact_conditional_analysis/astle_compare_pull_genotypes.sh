#!/bin/bash
export CHR=$1
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
source ./envars
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment
export MODULEPATH=/rds/project/who1000-1/rds-who1000-cbrc/user/cbrcmod/modules/out/modulefiles:$MODULEPATH 
module load qctool/beta

export SAMPLE_FILE_FH=${SAMPLE_FILE}

if [ $CHR = "X" ]
then
    SAMPLE_FILE_FH=/rds/project/who1000-1/rds-who1000-cbrc/user/wja24/shared/parsa_only/ukbb_pheno/ukbb_500k_final_adjusted_chr${CHR}.sample
fi
if [ $CHR = "XY" ]
then
    SAMPLE_FILE_FH=/rds/project/who1000-1/rds-who1000-cbrc/user/wja24/shared/parsa_only/ukbb_pheno/ukbb_500k_final_adjusted_chr${CHR}.sample
fi 
echo $SAMPLE_FILE_FH

if [ $? -ne 0 ]; then { echo "qctool bgen failed, aborting." ; exit 1; } fi

qctool_v2.0-rc9 \
-g ${IMPUTE_BGEN_FOLDER}/ukb_imp_chr${CHR}_v3.bgen \
-s ${SAMPLE_FILE_FH} \
-incl-positions ${OUT_DIR}/astle_compare/pull_tables/${CHR}.csv \
-threads 2 \
-og ${OUT_DIR}/astle_compare/genfiles/chr${CHR}.gen

if [ $? -ne 0 ]; then { echo "qctool gen failed, aborting." ; exit 1; } fi

awk '{gsub("^0*", "", $1); $2=$1":"$4"_"$5"_"$6; $3=$1":"$4"_"$5"_"$6; print}' ${OUT_DIR}/astle_compare/genfiles/chr${CHR}.gen > ${OUT_DIR}/astle_compare/genfiles/chr${CHR}_clean_id.gen

if [ $? -ne 0 ]; then { echo "awk failed, aborting." ; exit 1; } fi

rm ${OUT_DIR}/astle_compare/genfiles/chr${CHR}.gen

if [ $? -ne 0 ]; then { echo "rm failed, aborting." ; exit 1; } fi

awk '{print $2}' ${OUT_DIR}/astle_compare/genfiles/chr${CHR}_clean_id.gen > ${OUT_DIR}/astle_compare/genfiles/ids_${CHR}.tsv

if [ $? -ne 0 ]; then { echo "awk failed, aborting." ; exit 1; } fi

python ./astle_compare_check_genfiles.py --chrom ${CHR}

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_2/gen2hd5_lib/

/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_2/gen2hd5.git/src/gen2hd5 \
    -g ${OUT_DIR}/astle_compare/genfiles/chr${CHR}_clean_id.gen \
    -o ${OUT_DIR}/astle_compare/genfiles/chr${CHR}_clean_id.h5

if [ $? -ne 0 ]; then { echo "hd5 failed, aborting." ; exit 1; } fi

echo completed
