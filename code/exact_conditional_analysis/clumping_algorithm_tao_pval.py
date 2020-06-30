import pandas as pd
import numpy as np
import pdb
import os
OUT_DIR = os.environ['OUT_DIR']

clumps = pd.read_csv(OUT_DIR+"/condout/ldclump_dosage/hg19_clump_table_tao.tsv")
assert(sum(clumps['VARIANT'].duplicated()) == 0)
clumps['max_mlog10pval'] = np.nan
clumps['max_trait'] = np.nan
final_output = pd.read_csv(OUT_DIR+"/BCX_final_output_raw.csv")

for index, row in clumps.iterrows():
    final_output_sub = final_output.loc[final_output['VARIANT'] == row['VARIANT'],:]
    clumps.loc[index, "max_mlog10pval"] = final_output_sub.loc[final_output_sub['MLOG10P'].idxmax(),"MLOG10P"]
    clumps.loc[index, "max_trait"] = final_output_sub.loc[final_output_sub['MLOG10P'].idxmax(),"trait"]
    final_output = final_output.loc[final_output['VARIANT'] != row['VARIANT'],:]

clumps['sentinel'] = False
clumps.loc[clumps.groupby('clump_id')['max_mlog10pval'].idxmax(),"sentinel"] = True
assert(len(clumps.loc[clumps['sentinel'] == True,"clump_id"].unique()) == len(clumps['clump_id'].unique()) )
clumps.to_csv(OUT_DIR+"/condout/ldclump_dosage/hg19_clump_table_tao_sentinel.tsv", sep="\t", index=None)
