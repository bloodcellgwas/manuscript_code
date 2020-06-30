#!/software/python-2.7.10/bin/python2.7
# give it a range of ids
import os
import sys
import subprocess
import re
import time
import pdb
import pandas as pd
OUT_DIR = os.environ['OUT_DIR']
LOG_DIR = os.environ['LOG_DIR']

f = open(OUT_DIR+"/blocks/block_chrs_order.tsv", 'r')
blocks_array = []
blocks_array_chr = []
blocks_array_chr_numeric = []
for line in f:
    row = line.split(" ")
    blocks_array.append(row[0])
    blocks_array_chr.append(row[1].replace("\n", ""))
    blocks_array_chr_numeric.append(row[2].replace("\n", ""))

def get_lsb_index_from_block(blocks):
    lsb_index_list = []
    for block in blocks:
        lsb_index_list.append(blocks_array.index(str(block)))
    return(lsb_index_list)

def check_chromosome_exists(lsb_jobindex, blocks_array_chr):
    chromosome = blocks_array_chr[lsb_jobindex]
    if os.path.exists(OUT_DIR+"/genfiles/chr"+chromosome+"_gwsig_clean_id.gen"):
        return(True)
    else:
        return(False)

def get_chromosome(lsb_jobindex):
    chromosome = str(blocks_array_chr[lsb_jobindex])
    return(chromosome)

def get_chromosome_numeric(lsb_jobindex):
    chromosome = str(blocks_array_chr_numeric[lsb_jobindex])
    return(chromosome)

def get_block_size(block_id):
    #### FINISH THIS
    block=pd.read_csv(OUT_DIR+"/blocks/pull_ids_block_"+str(block_id)+".tsv")
    return(block.shape[0])

# make dataframe of block sizes
df_lsb_jobindex = []
df_block_id = []
df_block_size = []
df_chr_num = []
for lsb_jobindex in range(1, len(blocks_array_chr_numeric)):
    df_lsb_jobindex.append(lsb_jobindex)
    df_block_id.append(blocks_array[lsb_jobindex])
    df_chr_num.append(get_chromosome_numeric(lsb_jobindex))
    df_block_size.append(get_block_size(blocks_array[lsb_jobindex]))

block_sizes = pd.DataFrame({'lsb_jobindex': df_lsb_jobindex, 'block_id': df_block_id, 'block_size': df_block_size, 'chr_num': df_chr_num})

block_sizes.to_csv(OUT_DIR+"/blocks/block_sizes_db.csv", index=False)
