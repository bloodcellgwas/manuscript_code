#!/bin/bash

/home/pa354/software/samtools/bin/bgzip -c ${OUT_DIR}/tables/mapping_file.vcf > ${OUT_DIR}/tables/mapping_file.vcf.gz
/home/pa354/software/samtools/bin/tabix -p vcf ${OUT_DIR}/tables/mapping_file.vcf.gz

/home/pa354/software/bcftools/bin/bcftools annotate -a /rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_2/support_files/dbsnp/human_9606_b150_GRCh37p13/00-All.vcf.gz -c ID ${OUT_DIR}/tables/mapping_file.vcf.gz > ${OUT_DIR}/tables/mapping_file_rsids.vcf
