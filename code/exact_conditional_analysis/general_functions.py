import sys
gf = sys.modules[__name__]
out_dir = "/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_4"
def decompose_variantid(variants):
    import pandas as pd
    if type(variants) is pd.core.series.Series: variants = variants.tolist()
    if isinstance(variants, list) == False: variants = [variants]
    chromosome = []
    bp = []; ref = []; alt = []
    for variant in variants:
        chromosome.append(variant.split(":")[0])
        bp.append(variant.split(":")[1].split("_")[0])
        ref.append(variant.split(":")[1].split("_")[1])
        alt.append(variant.split(":")[1].split("_")[2])
    data = {'VARIANT':variants, 'chromosome':chromosome, 'bp':bp, 'ref':ref, 'alt':alt}
    result = pd.DataFrame(data=data)
    return(result)

def get_trait_names():
    import pandas as pd;
    traits = pd.read_csv(gf.out_dir+"/pheno_names.tsv", header=None)
    return(traits[0].tolist())

def get_condsig_variants(trait):
    import pandas as pd
    condsig_df = pd.read_csv(gf.out_dir+"/condout/results/condsig_"+trait+"_gwas_normalised.tsv", sep="\t")
    return(condsig_df["VARIANT"].tolist())

def get_all_condsig_variants(chrom=None):
    condsig = []
    for trait in gf.get_trait_names():
        condsig = list(set(gf.get_condsig_variants(trait) + list(condsig)))
    if chrom is not None:
        import re
        if chrom == 23:
            regexp = re.compile(r'^X')
        else:
            regexp = re.compile(r'^'+str(chrom)+':')
        return([s for s in condsig if regexp.search(s)])
    else: return(condsig)

def get_variant_altfreq(variants, condsig_only=False):
    import pandas as pd
    import numpy as np
    variant_df = gf.decompose_variantid(variants)
    variant_id = []
    altfreq = []
    for chrom in variant_df["chromosome"].unique():
        if condsig_only == True:
            file_loc = "/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_meta_analysis/snpstats/snpstats_"+str(chrom)+".csv"
            snpstats_file = pd.read_csv(file_loc, dtype=str)
        else:
            file_loc = "/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_2/snpstats/ukb_impv3_chr"+str(chrom)+"_snpstats.txt"
            snpstats_file = pd.read_csv(file_loc, dtype=str, sep=" ", usecols=['alternate_ids', 'alleleB_frequency'], skiprows=8)
        snpstats_file = snpstats_file.loc[snpstats_file["alternate_ids"].isin(variants),:]
        for i, row in variant_df.loc[variant_df["chromosome"] == chrom,:].iterrows():
            variant = row["VARIANT"]
            variant_id.append(variant)
            get_freq = snpstats_file.loc[snpstats_file["alternate_ids"] == variant, "alleleB_frequency"]
            if len(get_freq) == 0:
                altfreq.append(np.nan)
            elif len(get_freq) == 1:
                altfreq.append(snpstats_file.loc[snpstats_file["alternate_ids"] == variant, "alleleB_frequency"].values[0])
            else:
                exit(row+"there seems to be duplicates in the snpstats file")
    table = pd.DataFrame({'VARIANT': variant_id, 'ALT_FREQ': altfreq})
    # make sure order is same as input
    table.set_index('VARIANT', inplace=True)
    return(table["ALT_FREQ"].astype('float64').values)


def get_all_astle_variants(chrom=None):
    import pandas as pd
    astle_table = pd.read_csv(gf.out_dir+"/astle_compare/mmc4.csv")
    if chrom is not None:
        astle_table = astle_table.loc[astle_table['Chr (GRCh37)'] == chrom,:]
    return(astle_table['Unique Variant ID'].tolist())

def decompose_variantid(variants, sort=False):
    import pandas as pd
    if type(variants) is pd.core.series.Series: variants = variants.tolist()
    if type(variants) is pd.core.indexes.base.Index: variants = variants.tolist()
    if isinstance(variants, list) == False: variants = [variants]
    chromosome = []
    bp = []; ref = []; alt = []
    order = []; chromosome_num = []
    for variant in variants:
        chromosome.append(variant.split(":")[0])
        bp.append(variant.split(":")[1].split("_")[0])
        ref.append(variant.split(":")[1].split("_")[1])
        alt.append(variant.split(":")[1].split("_")[2])
        variant_chrom_num = variant.split(":")[0]
        if variant_chrom_num.isdigit() == False: variant_chrom_num = 23
        chromosome_num.append(variant_chrom_num)

    data = {'VARIANT':variants, 'chromosome':chromosome, 'chromosome_num': chromosome_num, 'bp':bp, 'ref':ref, 'alt':alt}
    result = pd.DataFrame(data=data)
    result[['bp']] = result[['bp']].apply(pd.to_numeric)
    if sort == True:
        result["sort"] = result["chromosome_num"].astype('int')*10**9 + result["bp"]
        result = result.sort_values(by=['sort'], ascending=False)
        result = result.drop(columns=['sort'])
    result['chromosome_num'] = result['chromosome_num'].astype('int')
    return(result)
