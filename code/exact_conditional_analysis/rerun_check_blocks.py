import os
import subprocess
OUT_DIR = os.environ['OUT_DIR']
lsbs = [5321, 5306, 5362, 5369, 5375, 5376, 5371, 5353, 5434, 5437, 5381, 5331, 5317, 5383, 5314, 5382, 5332, 5320, 5315, 5361, 5368, 5347, 5352, 5344, 5345, 5276, 5301, 5291, 5290, 5293, 5300, 5284, 5285, 5289, 5307, 5354, 5302, 5268, 5265, 5308, 5292, 5275, 5339]

f = open(OUT_DIR+"/blocks/block_chrs_order.tsv", 'r')
blocks_array = []
blocks_array_chr = []
blocks_array_chr_numeric = []
for line in f:
    row = line.split(" ")
    blocks_array.append(row[0])
    blocks_array_chr.append(row[1].replace("\n", ""))
    blocks_array_chr_numeric.append(row[2].replace("\n", ""))

def get_chromosome(lsb_jobindex):
    chromosome = str(blocks_array_chr[lsb_jobindex])
    return(chromosome)

for lsb_jobindex in lsbs:
    block_id = str(blocks_array[lsb_jobindex])
    chrom = get_chromosome(lsb_jobindex)
    if chrom != "XY" and chrom != "X": exit("this isn't sex chromosome "+str(lsb_jobindex)+" "+str(chrom))
    cmd_str = "Rscript ./check_block_subset.R "+block_id
    proc = subprocess.Popen(cmd_str, shell=True, executable='/bin/bash')
    proc.wait()
    if proc.returncode != 0: exit(block_id+" failed")
    
    cmd_str = "rm output/cond_analysis_extract_gen/XY/block_"+str(lsb_jobindex)+"_*"
    proc = subprocess.Popen(cmd_str, shell=True, executable='/bin/bash')
    proc.wait()

