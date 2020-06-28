#!/usr/local/bin/bash

## compute LD with LD Store - run with interval range then filter only needed SNPs per block

#arguments
index=$1
chr=$2
start=$3
end=$4
out_file=$5
mem=${6:-0}
thr=${7:-15}
root_dir=${8:-/lustre/scratch115/projects/ukbb500k_t151/FINEMAP_500K/ld_all_regions}

bgen_dir=/lustre/scratch115/projects/ukbiobank/FullRelease/Imputed/EGAD00010001474



if [ $start -lt 0 ]
then
start=1
fi


let diff=$end-$start
echo $diff

#run LD store with 15 threads as default or otherwise as specified by user
if [ $mem -eq 0 ]
then
if [ $chr -eq 6 ] || [ $diff -gt 800000 ]
then 
### if we are on chr 6 or the region is big, run with extra memory 50GB
bsub -G ukbb500k_t151 -J "LD_${index}" -o ${root_dir}/log_files/LD_$index.log -R'select[mem>50000] rusage[mem=50000] span[hosts=1]' -M 50000 -n${thr} -- /nfs/team151/software/ldstore_v1.1_x86_64/ldstore --bgen ${bgen_dir}/ukb_imp_chr${chr}_v3.bgen --incl-range ${start}-${end} --bcor ${root_dir}/${out_file} --n-threads ${thr} --n-variants-chunk 50
else
### if not on chr 6 and region not large, run with standard memory 30GB
bsub -G ukbb500k_t151 -J "LD_${index}" -o ${root_dir}/log_files/LD_$index.log -R'select[mem>30000] rusage[mem=30000] span[hosts=1]' -M 30000 -n${thr} -- /nfs/team151/software/ldstore_v1.1_x86_64/ldstore --bgen ${bgen_dir}/ukb_imp_chr${chr}_v3.bgen --incl-range ${start}-${end} --bcor ${root_dir}/${out_file} --n-threads ${thr} --n-variants-chunk 50
fi
else
### custom run with memory specified by user
if [ $diff -gt 1000000 ]
then
bsub -G ukbiobank_t151 -J "LD_${index}" -q basement -o ${root_dir}/log_files/LD_$index.log -R"select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" -M ${mem} -n${thr} -- /nfs/team151/software/ldstore_v1.1_x86_64/ldstore --bgen ${bgen_dir}/ukb_imp_chr${chr}_v3.bgen --incl-range ${start}-${end} --bcor ${root_dir}/${out_file} --n-threads ${thr} --n-variants-chunk 50
else
bsub -G ukbiobank_t151 -J "LD_${index}" -q long -o ${root_dir}/log_files/LD_$index.log -R"select[mem>${mem}] rusage[mem=${mem}] span[hosts=1]" -M ${mem} -n${thr} -- /nfs/team151/software/ldstore_v1.1_x86_64/ldstore --bgen ${bgen_dir}/ukb_imp_chr${chr}_v3.bgen --incl-range ${start}-${end} --bcor ${root_dir}/${out_file} --n-threads ${thr} --n-variants-chunk 50
fi
fi

#merge output in single file 
bsub -G ukbb500k_t151 -w "done(LD_${index})" -J "mer_${index}" -o  ${root_dir}/log_files/mer_$index.log -R'select[mem>10000] rusage[mem=10000]' -M 10000 -- /nfs/team151/software/ldstore_v1.1_x86_64/ldstore --bcor ${root_dir}/${out_file} --merge ${thr}



