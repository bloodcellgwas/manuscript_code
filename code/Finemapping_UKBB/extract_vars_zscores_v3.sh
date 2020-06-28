#!/usr/local/bin/bash

# cat finemap_overview_blocks.txt | while read -a phenotypes
# do 
# phen_arg=${phenotypes[0]}
# Nblocks=${phenotypes[2]}
# echo "sh ~/Desktop/Scripts/BTRU-Theme1/FINEMAP_pipeline/extract_vars_zscores_v3.sh ${phen_arg} ${phen_arg}" | bsub -G ukbb500k_t151 -J "Ex_${phen_arg}[1-${Nblocks}]" -o /lustre/scratch115/projects/ukbb500k_t151/final_fine_mapping/${phen_arg}/log_files/extract_${phen_arg}.log -R'select[mem>2000] rusage[mem=2000]' -M2000
# done


#arguments
root_dir=/lustre/scratch115/projects/ukbb500k_t151/final_fine_mapping
res_dir=/lustre/scratch115/realdata/mdt3/projects/ukbb500k_t151/BOLT_analysis/output
pheno=$1
pheno_alt=$2
pheno_res_dir=${res_dir}

#Block range
b=${LSB_JOBINDEX}
chr=$(awk '{print $1}' ${root_dir}/${pheno}/${pheno}_block_ranges/range_str_end_${pheno}_finemap_block_${b}.txt)
start=$(awk '{print $2}' ${root_dir}/${pheno}/${pheno}_block_ranges/range_str_end_${pheno}_finemap_block_${b}.txt)
end=$(awk '{print $3}' ${root_dir}/${pheno}/${pheno}_block_ranges/range_str_end_${pheno}_finemap_block_${b}.txt)

#extract block from gwas results, filter for maf and info score
cd ${pheno_res_dir}
awk -v start="$start" -v end="$end" -v chr="$chr" '{if($2==chr && $3>=start-250000 && $3<=end+250000 && $7<=1-0.00005 && $7>=0.00005 && $8>=0.4) print $0 "\t" $2 ":" $3 "_" $5 "_" $6}' ${pheno}_gwas_normalised_imputed_full_panel.out > ${root_dir}/${pheno}/${pheno}_block_meta_snps/meta_snps_${pheno}_block_${b}.txt


#create snp list with variant id 
cd ${root_dir}/${pheno}/${pheno}_block_meta_snps/
awk '{print $2 ":" $3 "_" $5 "_" $6}' meta_snps_${pheno}_block_${b}.txt > ../finemap_snps_gen/${pheno}_block_${b}_snplist.txt

#create .z files in new format 
cut -f1-3,5-7,11,12 ${root_dir}/${pheno}/${pheno}_block_meta_snps/meta_snps_${pheno}_block_${b}.txt > ${root_dir}/${pheno}/finemap_${pheno}/${pheno}_block${b}.z

# change from tab to space delimited
sed 's/\t/ /g' ${root_dir}/${pheno}/finemap_${pheno}/${pheno}_block${b}.z > tmp_${pheno}_${b}
mv tmp_${pheno}_${b} ${root_dir}/${pheno}/finemap_${pheno}/${pheno}_block${b}.z
#compute maf
awk '{if($6>0.5) {print $1" "$2" "$3" "$4" "$5" "1-$6" "$7" "$8} else if($6<=0.5) print $0}' ${root_dir}/${pheno}/finemap_${pheno}/${pheno}_block${b}.z > tmp_${pheno}_${b}
mv tmp_${pheno}_${b} ${root_dir}/${pheno}/finemap_${pheno}/${pheno}_block${b}.z
#print 0 in front of chr 1-9
awk '{if($2<10) {print $1" "0$2" "$3" "$4" "$5" "$6" "$7" "$8} else if($2>=10) print $1" "$2" "$3" "$4" "$5" "$6" "$7" "$8}' ${root_dir}/${pheno}/finemap_${pheno}/${pheno}_block${b}.z  > tmp_${pheno}_${b}
cat ${root_dir}/tmp_hedd tmp_${pheno}_${b} > ${root_dir}/${pheno}/finemap_${pheno}/${pheno}_block${b}.z

rm tmp_${pheno}_${b}   







