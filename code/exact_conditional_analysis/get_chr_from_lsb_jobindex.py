#!/software/python-2.7.10/bin/python2.7
# give it a range of ids
import os
import sys
import subprocess
import re
import time
OUT_DIR = os.environ['OUT_DIR']
LOG_DIR = os.environ['LOG_DIR']

f = open(OUT_DIR+"/blocks/block_chrs_order.tsv", 'r')
blocks_array = []
blocks_array_chr = []
for line in f:
    row = line.split(" ")
    blocks_array.append(row[0])
    blocks_array_chr.append(row[1].replace("\n", ""))

def check_chromosome_exists(lsb_jobindex, blocks_array_chr):
    chromosome = blocks_array_chr[lsb_jobindex]
    if os.path.exists(OUT_DIR+"/genfiles/chr"+chromosome+"_gwsig.gen"):
        return(True)
    else:
        return(False)

def get_chromosome(lsb_jobindex):
    chromosome = str(blocks_array_chr[lsb_jobindex])
    return(chromosome)

print("Chromosome: "+get_chromosome(int(sys.argv[1])))
print("Block: "+str(blocks_array[int(sys.argv[1])]))
