## Files/Folder prep
#!/usr/local/bin/bash

# echo 'sh ~/Desktop/Scripts/BTRU-Theme1/FINEMAP_pipeline/folder_prep.sh' | bsub -G ukbiobank_t151 -J 'folders[1-43]' -o /lustre/scratch115/projects/ukbiobank_t151/fine-mapping/FINEMAP_sysmex_2/folder_prep.log -R'select[mem>100] rusage[mem=100]' -M100
 

traitnumber=${LSB_JOBINDEX}
awk 'NR=='$traitnumber /lustre/scratch115/projects/ukbiobank_t151/fine-mapping/FINEMAP_sysmex_2/final_finemap_overview_blocks.txt | while read -a line
do
block=${line[2]} ## number of blocks per phenotype
pheno=${line[0]} ## phenotype name
root_dir=/lustre/scratch115/projects/ukbiobank_t151/fine-mapping/FINEMAP_sysmex_2/

mkdir -p ${root_dir}/${pheno}/${pheno}_block_meta_snps
mkdir -p ${root_dir}/${pheno}/${pheno}_block_ranges
mkdir -p ${root_dir}/${pheno}/${pheno}_conditional_leads
mkdir -p ${root_dir}/${pheno}/finemap_${pheno}
mkdir -p ${root_dir}/${pheno}/finemap_snps_gen
mkdir -p ${root_dir}/${pheno}/log_files

mv ${root_dir}/${pheno}/range* ${root_dir}/${pheno}/${pheno}_block_ranges
mv ${root_dir}/${pheno}/${pheno}_finemap_block* ${root_dir}/${pheno}/${pheno}_conditional_leads/
done


