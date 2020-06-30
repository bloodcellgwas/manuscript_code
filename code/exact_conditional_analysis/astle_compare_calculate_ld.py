#!python
import os
import h5py
import argparse
import numpy as np
import pandas as pd
import general_functions

OUT_DIR = os.environ['OUT_DIR']
parser = argparse.ArgumentParser()
parser.add_argument('--chrom', type=int, help='chromosome', default=None)
args = parser.parse_args()
print(args)

print(args.chrom)
hd5_variants = pd.read_csv(OUT_DIR+"/astle_compare/genfiles/ids_"+str(args.chrom)+".tsv", header=None)[0]
chr_data = h5py.File(OUT_DIR+"/astle_compare/genfiles/chr"+str(args.chrom)+"_clean_id.h5", "r")

# the hd5 file should have an intercept
if int(sum(chr_data['add'][0,:])) == chr_data['add'].shape[1]:
    assert(hd5_variants.shape[0] == chr_data['add'].shape[0]-1)
else:
    exit("check to see if hd5 generation has changed because the intercept is missing")
    assert(hd5_variants.shape[0] == chr_data['add'].shape[0])
    
# remove intercept
chr_data_subset = np.delete(chr_data['add'].value, 0, axis=0)
chr_data.close()

chrom_condsig = general_functions.get_all_condsig_variants(args.chrom)
# load all astle variants on this chromosome
chrom_astle = general_functions.get_all_astle_variants(args.chrom)

variants_select = np.where(hd5_variants.isin(list(set(chrom_condsig + chrom_astle))))[0]
extracted_variants = hd5_variants[variants_select].astype(np.str).tolist()

# check we have all the condsig variants
assert(sum([True for s in chrom_condsig if s in extracted_variants]) == len(chrom_condsig))

# subset the dataset
chr_data_subset = chr_data_subset[variants_select,]

# crosscheck the mafs
if args.chrom != 23:
    snpstats_mafs = general_functions.get_variant_altfreq(extracted_variants, condsig_only=True)
    for vari in range(0, len(extracted_variants)):
        var = extracted_variants[vari]
        snpstats_maf = snpstats_mafs[vari]
        if np.isnan(snpstats_maf): continue
        calc_maf = np.mean(chr_data_subset[vari,]+1)/2
        diff = abs((snpstats_maf - calc_maf)/snpstats_maf) * 100
        if diff > 1: exit(str(vari) + "  " + extracted_variants[vari] + " difference too large: "+str(diff))
else: print("skipping MAF check")

all_data = chr_data_subset
varnames = extracted_variants

corrmat = np.corrcoef(all_data)
assert(sum(np.diagonal(corrmat)) == len(varnames))
print("calculated correlation matrix")

corrmat = pd.DataFrame(corrmat)
corrmat.columns = varnames
corrmat.index = varnames

# now melt back into table of correlations
corrdf = corrmat.stack().reset_index()
corrdf.columns = ['SNP_A', 'SNP_B', 'R2']
print("created correlation table")
a_split1 = corrdf['SNP_A'].str.split(":", n=1, expand=True)
a_split2 = a_split1[1].str.split("_", n=1, expand=True)
b_split1 = corrdf['SNP_B'].str.split(":", n=1, expand=True)
b_split2 = b_split1[1].str.split("_", n=1, expand=True)

corrdf['CHR_A'] = a_split1[0]
corrdf['BP_A'] = a_split2[0]
corrdf['CHR_B'] = b_split1[0]
corrdf['BP_B'] = b_split2[0]
print("splits complete")

corrdf = corrdf[["CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2"]]
corrdf.to_csv(OUT_DIR+"/astle_compare/ld/ld_chr"+str(args.chrom)+".ld", index=False)

