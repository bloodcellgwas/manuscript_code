#!/bin/bash
pheno=$1
DIR=/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/vep/VEP_84_GRCh37
LOFFILES=${DIR}/lof_files
PYTHON=/software/bin/python
export PATH=/home/pa354/software/samtools/bin/:$PATH
# Path VARIABLE
VEPSoft=/home/pa354/software/ensembl-tools-release-84/scripts/variant_effect_predictor/variant_effect_predictor.pl
PYTHON=python2.7

echo $DIR
echo $LOFFILES

$VEPSoft \
    --offline \
    --cache \
    --force_overwrite \
    --dir_cache ${DIR} \
    --dir ${DIR} \
    --i $OUT_DIR/vcf/${pheno}_gwsig.vcf \
    --allele_number \
    --plugin LoF,human_ancestor_fa:${LOFFILES}/human_ancestor.fa,filter_position:0.05,${LOFFILES}/phylocsf.sql \
    --vcf \
    --o ${OUT_DIR}/vep/${pheno}_gwsig_vep.vcf

if [ $? -ne 0 ]; then { echo "Failed, aborting." ; exit 1; } fi
echo Run VEP 84 ... built 37 with Loftee plugin :LoF,Human_ancestor_fa,filter_position-Transcript,phylocsf.sql... Done
echo "_____________"


### 2nd step:  Python to extract the CSQ and put into a table ###
${PYTHON} \
    ${LOFFILES}/loftee-master/src/tableize_vcf.py \
    --vcf ${OUT_DIR}/vep/${pheno}_gwsig_vep.vcf \
    --vep_info Allele,Consequence,IMPACT,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,ALLELE_NUM,DISTANCE,STRAND,SYMBOL_SOURCE,HGNC_ID,LoF_info,LoF_filter,LoF_flags,LoF \
    --functional_simplify \
    --include_id \
    --output ${OUT_DIR}/vep/tabular_${pheno}.txt
if [ $? -ne 0 ]; then { echo "Failed, aborting." ; exit 1; } fi
echo python extract from vcf VEP loftee... Done
echo "_____________"
