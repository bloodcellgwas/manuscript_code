#!/bin/bash

#export CHR=$SLURM_ARRAY_TASK_ID
export CHR=$1
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
source /home/pa354/projects/2018_10_11_ukbb500k_2/envars
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment
export MODULEPATH=/rds/project/who1000-1/rds-who1000-cbrc/user/cbrcmod/modules/out/modulefiles:$MODULEPATH 
module load qctool/beta


if [ $1 -gt 22 ]
then
    qctool_v2.0-rc9 \
	-g ${OUT_DIR}/genfiles/chrXY_merged.gen \
	-s ${OUT_DIR}/genfiles/chrXY_merged.sample \
	-incl-snpids ${OUT_DIR}/condout/pull_ids/pull_ids_chr${CHR}.tsv \
	-og ${OUT_DIR}/condout/subgenhd5/chr${CHR}_pos_clean_id.gen

    if [ $? -ne 0 ]; then { echo "qctool failed, aborting." ; exit 1; } fi

    awk '{gsub("^0*","", $1); print $1":"$4"_"$5"_"$6}' ${OUT_DIR}/condout/subgenhd5/chr${CHR}_pos_clean_id.gen > ${OUT_DIR}/condout/subgenhd5/ids_pos_${CHR}.tsv

    if [ $? -ne 0 ]; then { echo "awk failed, aborting." ; exit 1; } fi
    # we are not pulling by position so use a different checking script than below
    # where we are pulling by position
    Rscript ./pull_genotypes_condsig_check.R ${CHR}

    if [ $? -ne 0 ]; then { echo "extract CHECK failed, aborting." ; exit 1; } fi

fi
if [ $1 -lt 23 ]
then
    qctool_v2.0-rc9 \
	-g ${IMPUTE_BGEN_FOLDER}/ukb_imp_chr${CHR}_v3.bgen \
	-s ${SAMPLE_FILE} \
	-incl-positions ${OUT_DIR}/condout/pull_ids/pull_ids_chr${CHR}_pos.tsv \
	-og ${OUT_DIR}/condout/subgenhd5/chr${CHR}_pos.gen

    if [ $? -ne 0 ]; then { echo "qctool failed, aborting." ; exit 1; } fi

    awk '{gsub("^0*", "", $1); $2=$1":"$4"_"$5"_"$6; $3=$1":"$4"_"$5"_"$6; print}' ${OUT_DIR}/condout/subgenhd5/chr${CHR}_pos.gen > ${OUT_DIR}/condout/subgenhd5/chr${CHR}_pos_clean_id.gen

    if [ $? -ne 0 ]; then { echo "awk failed, aborting." ; exit 1; } fi

    #rm ${OUT_DIR}/condout/subgenhd5/chr${CHR}_pos.gen

    if [ $? -ne 0 ]; then { echo "rm failed, aborting." ; exit 1; } fi

    awk '{print $2}' ${OUT_DIR}/condout/subgenhd5/chr${CHR}_pos_clean_id.gen > ${OUT_DIR}/condout/subgenhd5/ids_pos_${CHR}.tsv
		
    if [ $? -ne 0 ]; then { echo "awk failed, aborting." ; exit 1; } fi
    
    Rscript ./check_genfiles_pos.R ${CHR}
fi

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_2/gen2hd5_lib/

/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_2/gen2hd5.git/src/gen2hd5 \
    -g ${OUT_DIR}/condout/subgenhd5/chr${CHR}_pos_clean_id.gen \
    -o ${OUT_DIR}/condout/subgenhd5/chr${CHR}_pos_clean_id.h5

if [ $? -ne 0 ]; then { echo "hd5 failed, aborting." ; exit 1; } fi

echo completed
