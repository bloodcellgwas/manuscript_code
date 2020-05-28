#!/bin/bash

use Samtools

ggsashimi="/broad/sankaranlab/tools/ggsashimi/sashimi-plot.py"
rsid="rs139178017"
chr=7
pos=100225847
start=$(($pos - 300))
end=$(($pos + 125))

# # rs8113779
# rsid="rs8113779"
# chr=19
# pos=45910003
# start=$(($pos - 150))
# end=$(($pos + 500))

rsid="rs12898397"
chr=15
pos=75130093
start=$(($pos - 425))
end=$(($pos + 80))

coord="${chr}:${start}-${end}"
controlbam="/broad/sankaranlab/clareau/extract164_merged_bams/${rsid}.merge.control.bam"
casebam="/broad/sankaranlab/clareau/extract164_merged_bams/${rsid}.merge.case.bam"
bamlist="/broad/sankaranlab/BCX_FINEMAP/ukb_500k_v2/splice_variants/${rsid}_bamlist.tsv"
outfile="/broad/sankaranlab/BCX_FINEMAP/ukb_500k_v2/splice_variants/splice_graphs/${rsid}_ggsashimi.pdf"

echo -e "carriers\t${casebam}" > $bamlist
echo -e "non-carriers\t${controlbam}" >> $bamlist

python $ggsashimi \
-b ${bamlist} \
--coordinates ${coord} \
-g /broad/sankaranlab/tools/gene_annotations/gencode.v34lift37.basic.annotation.nochr_prefix.gtf \
-M 100 \
--alpha 1 \
-C 1 \
-P /broad/sankaranlab/tools/ggsashimi/examples/palette.txt \
-F pdf \
--base-size=12 \
--height=1.3 --width=4.5 \
--ann-height 1 \
-o ${outfile}

