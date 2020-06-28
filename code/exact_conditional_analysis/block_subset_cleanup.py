import pandas as pd
import os
import subprocess
OUT_DIR = os.environ['OUT_DIR']

block_sizes_db = pd.read_csv(OUT_DIR+"/blocks/block_sizes_db.csv")

f = open("./check_block_subset_failed.txt")
for line in f:
    lsb_jobindex = int(line.split("_")[1])
    print(lsb_jobindex)
    block_row = block_sizes_db.loc[block_sizes_db["lsb_jobindex"] == int(lsb_jobindex),:] 
    block_id =str(block_row['block_id'].values[0])
    assert(block_row.shape[0] == 1)
    cmd_str = "Rscript ./check_block_subset.R "+block_id
    try:
        subprocess.check_call(cmd_str)
        #proc = subprocess.Popen(cmd_str, shell=True, executable='/bin/bash')
    except:
        print("error")
        # if theres an error delete gen and id file
        if os.path.isfile(OUT_DIR+'/blocks/block_'+block_id+".gen") == True:
            print("deleting gen file")
        if os.path.isfile(OUT_DIR+'/blocks/ids_'+block_id+".tsvn") == True:
            print("deleting ids file")

        # check and delete if necessary
    

#            os.remove(OUT_DIR+'/blocks/block_'+block_id++"/.gen")
