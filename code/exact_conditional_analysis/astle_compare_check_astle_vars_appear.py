'''
A script to check all the astle condsig variants have a proxy or appear in the current ukbb500k GWAS
'''
import pdb
import numpy as np
import pandas as pd
import general_functions as gf

astle_condsig = pd.read_csv(gf.out_dir+"/astle_compare/mmc4.csv")
ukbb500k = pd.read_csv(gf.out_dir+"/BCX_final_output.tsv", sep="\t")
ukbb500k = ukbb500k.loc[~ukbb500k.iloc[:,0].isin(['MRV', 'MSCV']),:]

# check all ukbb500k traits exist in astle_condsig
ukbb500k_traits = ukbb500k.iloc[:,0].unique().tolist()
astle_traits = astle_condsig.iloc[:,0].unique().tolist()

dont_exist = np.setdiff1d(ukbb500k_traits, astle_traits)
assert(len(dont_exist) == 0)


astle_variant = []
astle_trait = []
ukbb500k_variant = []
r2 = []
for chrom in range(1,23):
    print(chrom)
    ld_table = pd.read_csv(gf.out_dir+"/astle_compare/ld/ld_chr"+str(chrom)+".ld")
    # loop through all of the ukbb500k traits
    for trait in ukbb500k.iloc[:,0].unique():
        # loop through all the astle condsigs associated with this trait and check that they have
        #  a proxy in the ukbb500k data
        astle_condsig_sub = astle_condsig.loc[(astle_condsig.iloc[:,5] == chrom) & (astle_condsig.iloc[:,0] == trait),:]
        subset_a = (ld_table.loc[:,"SNP_A"].isin(astle_condsig_sub.iloc[:,3])) & (ld_table.loc[:,"SNP_B"].isin(ukbb500k.iloc[:,3]))
        subset_b = (ld_table.loc[:,"SNP_B"].isin(astle_condsig_sub.iloc[:,3])) & (ld_table.loc[:,"SNP_A"].isin(ukbb500k.iloc[:,3]))
        ld_table_sub = ld_table.loc[subset_a | subset_b,:]
        for ind, row in astle_condsig_sub.iterrows():
            subset_a = (ld_table_sub.loc[:,"SNP_A"] == row.iloc[3]) & (ld_table_sub.loc[:, "SNP_B"].isin(ukbb500k.iloc[:,3]))
            subset_b = (ld_table_sub.loc[:,"SNP_B"] == row.iloc[3]) & (ld_table_sub.loc[:, "SNP_A"].isin(ukbb500k.iloc[:,3]))
            ld_table_row = ld_table_sub.loc[subset_a | subset_b,:]
            ld_table_row = ld_table_row.loc[ld_table_row['R2'].idxmax(),:]
            astle_variant.append(row.iloc[3])
            astle_trait.append(row.iloc[0])
            r2.append(ld_table_row["R2"])

            # if both snps are the same append thi
            if ld_table_row['SNP_A'] == ld_table_row['SNP_B']:
                ukbb500k_variant.append(ld_table_row['SNP_B'])
            elif ld_table_row['SNP_A'] == row.iloc[3]:
                ukbb500k_variant.append(ld_table_row['SNP_B'])
            elif ld_table_row['SNP_B'] == row.iloc[3]:
                ukbb500k_variant.append(ld_table_row['SNP_A']) 
            
data = {"trait": astle_trait, "astle_variant": astle_variant, "ukbb500k_variant": ukbb500k_variant, "r2": r2}
data = pd.DataFrame(data)

import pdb; pdb.set_trace()
# add other columns of results to the data
astle_condsig = astle_condsig.iloc[:,[0,3,12,18,19,22]]
data = data.merge(astle_condsig, left_on=['trait', 'astle_variant'], right_on=['Associated Blood Index', 'Unique Variant ID'])

data.to_csv(gf.out_dir+"/astle_compare/astle_var_ukbb500k.csv", index=False)
