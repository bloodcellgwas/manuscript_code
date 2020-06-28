#!python
import os
import h5py
import argparse
import numpy as np
import pandas as pd
import general_functions as gf

OUT_DIR = os.environ['OUT_DIR']
parser = argparse.ArgumentParser()
parser.add_argument('--chrom', type=int, help='chromosome', default=None)
args = parser.parse_args()
print(args)

print(args.chrom)
hd5_variants = pd.read_csv(OUT_DIR+"/astle_compare/genfiles/ids_"+str(args.chrom)+".tsv", header=None)[0]

# load all condsig on this chromosome
chrom_condsig = gf.get_all_condsig_variants(args.chrom)
# load all astle variants on this chromosome
chrom_astle = gf.get_all_astle_variants(args.chrom)
chrom_extract = set(chrom_condsig + chrom_astle)

missing = [a for a in chrom_extract if a not in hd5_variants.tolist()]

if len(missing) > 0:
    os.rename(gf.out_dir+"/astle_compare/genfiles/ids_"+str(args.chrom)+".tsv", gf.out_dir+"/astle_compare/genfiles/ids_"+str(args.chrom)+".tsv_WRONG")
    os.rename(gf.out_dir+"/astle_compare/genfiles/chr"+str(args.chrom)+"_clean_id.gen", gf.out_dir+"/astle_compare/genfiles/chr"+str(args.chrom)+"_clean_id.gen_WRONG")
    os.rename(gf.out_dir+"/astle_compare/genfiles/chr"+str(args.chrom)+"_clean_id.h5", gf.out_dir+"/astle_compare/genfiles/chr"+str(args.chrom)+"_clean_id.h5_WRONG")
    print(missing)
    exit("THERE ARE MISSING VARIANTS")
else: print("check complete")

