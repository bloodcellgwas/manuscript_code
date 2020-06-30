#!/bin/bash
export PATH="/home/pa354/software/miniconda3/bin:$PATH"
source envars
source activate itr
#snakemake data/all_pull_genotypes_gwsig.txt \
#snakemake data/all_genotype_pulls.txt --rerun-incomplete \
#snakemake data/ukbb_compare_previous_number_signals.csv \
#snakemake data/condout/vep/tabular_all_condsig.txt \
#snakemake  data/locuszoom/output/done.txt --use-conda --jobs 99 \
snakemake data/BCX_final_output_raw.csv \
    --jobs 20 #\
#    --cluster "\
#        sbatch -A BRIDGE-CORE-SL2-CPU --nodes=1 --ntasks=1 --mail-type=FAIL -p skylake-himem \
#        --cpus-per-task {threads} --mem {resources.mem_mb}mb -J {rule} -t 25:00:00 \
#        --output=./output/{params.section}/{rule}_%J.out \
#        --error=./output/{params.section}/{rule}_%J.err"

