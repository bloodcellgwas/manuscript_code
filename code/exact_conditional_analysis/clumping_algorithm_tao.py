#!python
import os
import numpy as np
import pandas as pd
import general_functions as gf
import pdb

OUT_DIR = os.environ['OUT_DIR']

condsig = None
for trait in gf.get_trait_names():
    ## load trait condsigs
    condsig_df = pd.read_csv(gf.out_dir+"/condout/results/condsig_"+trait+"_gwas_normalised.tsv", sep="\t")
    if condsig is None:
        condsig = gf.decompose_variantid(condsig_df['VARIANT'])
    else:
        condsig = condsig.append(gf.decompose_variantid(condsig_df['VARIANT']))

condsig = condsig.loc[~condsig['VARIANT'].duplicated(),:]

def get_partners(ld, var):
    partnersA = ld.loc[(ld['SNP_A'] == var) & (ld['R2'] > 0.1),"SNP_B"]
    partnersB = ld.loc[(ld['SNP_B'] == var) & (ld['R2'] > 0.1),"SNP_A"]
    partners = set(partnersB.tolist() + partnersA.tolist())
    return(list(partners))

def clump_function(partners, clumps, var):
    # if no clumps exist make the first one
    if clumps is None:
        clumps =  pd.DataFrame({'VARIANT': partners+[var]})
        clumps['clump_id'] = 1
        return(clumps)

    # check if any partners already exist in a clump
    matching_clumps = []
    if len(partners) != 0:
        clumps_s = clumps.loc[clumps['VARIANT'].isin(partners),:]
        matching_clumps = clumps_s['clump_id'].unique().tolist()

    # if there are no matching clumps add this var and partners to new clump
    clumps_add = pd.DataFrame({'VARIANT': partners+[var]})
    if len(matching_clumps) == 0:
        clumps_add['clump_id'] = clumps['clump_id'].max()+1
    # if there is one matching clump merge into that one
    elif len(matching_clumps) == 1:
        clumps_add['clump_id'] = matching_clumps[0]
    # if there are more than one matching clumps, merge them all into a new clump
    elif len(matching_clumps) > 1:
        clumps_add['clump_id'] = clumps['clump_id'].max()+1 
        clumps.loc[clumps['clump_id'].isin(matching_clumps),"clump_id"] = clumps['clump_id'].max()+1
        
    # append clumps_add into the dataframe
    clumps = clumps.append(clumps_add, ignore_index=True)
    return(clumps)

condsig.loc[condsig['chromosome'].isin(['X', 'XY']),"chromosome"] = 23
clumps = None
for chrom in condsig['chromosome'].unique():
    print(str(chrom)+" started")
    ld = pd.read_csv(gf.out_dir+"/condout/ldclump_dosage/dosage_"+str(chrom)+".ld")
    ld = ld.loc[ld['SNP_A'].isin(condsig['VARIANT']),:]
    ld = ld.loc[ld['SNP_B'].isin(condsig['VARIANT']),:]
    ld = ld.loc[ld['SNP_A'] != ld['SNP_B'],:]

    condsig_s = condsig.loc[condsig['chromosome'] == chrom,:]

    for index, row in condsig_s.iterrows():
        partners = get_partners(ld, row['VARIANT'])
        clumps = clump_function(partners, clumps, row['VARIANT'])

clumps = clumps.drop_duplicates()
clumps.to_csv(OUT_DIR+"/condout/ldclump_dosage/hg19_clump_table_tao.tsv", index=None)
