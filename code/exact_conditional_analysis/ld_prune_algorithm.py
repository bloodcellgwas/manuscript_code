#!python
import os
import numpy as np
import pandas as pd
import general_functions as gf

hg19_clumps = pd.read_csv(gf.out_dir+"/condout/ldclump_dosage/hg19_clump_table.tsv", sep="\t")
hg19_clumps['novel'] = None
condsig_vars = gf.get_all_condsig_variants()
condsig_df = gf.decompose_variantid(condsig_vars)

for chrom in range(1,23):
    print(chrom)
    condsig_df_sub = condsig_df.loc[condsig_df['chromosome_num'] == chrom,:]
    ld_table = pd.read_csv(gf.out_dir+"/astle_compare/ld/ld_chr"+str(chrom)+".ld")
    for condsig_var in condsig_df_sub['VARIANT'].tolist():
        ld_table_sub = ld_table.loc[ld_table['R2'] > 0.8,:]
        subset_a = (ld_table_sub['SNP_A'] == condsig_var) & (ld_table_sub['SNP_B'].isin(gf.get_all_astle_variants(chrom)))
        subset_b = (ld_table_sub['SNP_B'] == condsig_var) & (ld_table_sub['SNP_A'].isin(gf.get_all_astle_variants(chrom)))
        ld_table_sub = ld_table_sub.loc[subset_a | subset_b,:]

        hg19_clumps_sub = hg19_clumps.loc[hg19_clumps["hg19"] == condsig_var.split("_")[0], "clump_id"]
        clump_id = hg19_clumps_sub.unique()
        hg19_novel_sub = hg19_clumps.loc[hg19_clumps["hg19"] == condsig_var.split("_")[0], "novel"]
        clump_novel = hg19_novel_sub.unique()
        assert(len(clump_id) == 1)
        assert(len(clump_novel) == 1)

        if ld_table_sub.shape[0] > 0:
            hg19_clumps.loc[hg19_clumps['clump_id'] == clump_id[0], "novel"] = False
        elif clump_novel[0] != False and ld_table_sub.shape[0] == 0:
            hg19_clumps.loc[hg19_clumps['clump_id'] == clump_id[0], "novel"] = True

# remove the sex chromosomes
remove_a = hg19_clumps['hg19'].str.contains('X:') == False
remove_b = hg19_clumps['hg19'].str.contains('XY:') == False
hg19_clumps = hg19_clumps.loc[remove_a & remove_b,:]
assert(hg19_clumps.loc[hg19_clumps['novel'].isna(),:].shape[0] == 0)
hg19_clumps.to_csv(gf.out_dir+"/astle_compare/astle_var_ukbb500k.csv", index=False)
