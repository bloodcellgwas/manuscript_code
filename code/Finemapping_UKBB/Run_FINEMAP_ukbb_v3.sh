#!/usr/local/bin/bash

## RUN FINEMAP

#arguments
pheno=$1
b=${LSB_JOBINDEX}
root_dir=/lustre/scratch115/projects/ukbb500k_t151/final_fine_mapping/${pheno}/finemap_${pheno}

N=$(wc -l /lustre/scratch115/projects/ukbb500k_t151/final_fine_mapping/${pheno}/${pheno}_conditional_leads/${pheno}_finemap_block_${b}.txt | cut -d" " -f1)
let N=N-1
if [ $N -gt 10 ]
then
N=10
fi
#sed -i 's/n_ind/n-ind/g' ${root_dir}/masterfile.txt 

cd ${root_dir}

/software/team151/finemap_v1.3.1_x86_64/finemap_v1.3.1_x86_64 --sss --in-files masterfile.txt --n-causal-snps ${N} --prior-std 0.08 --dataset ${LSB_JOBINDEX} 