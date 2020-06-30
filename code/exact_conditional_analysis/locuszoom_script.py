#!/software/python-2.7.10/bin/python2.7
import os
import subprocess
import sys
print(sys.version)

OUT_DIR = os.environ['OUT_DIR']
PHENO_NAMES = os.environ['PHENO_NAMES']
SAMPLE_FILE = os.environ['SAMPLE_FILE']
IMPUTE_BGEN_FOLDER = os.environ['IMPUTE_BGEN_FOLDER']
SUPPORT_FILES = os.environ['SUPPORT_FILES']
LSB_JOBINDEX = int(os.environ['SLURM_ARRAY_TASK_ID'])
#LSB_JOBINDEX = 1
phen_names = open(PHENO_NAMES, 'r')
phen_arr = phen_names.read().split("\n")

def extract_region(chromosome, fromb, tob):
    # CANT RUN THIS FUNCTION IN PARALLEL YET (plink always outputs to the same file)
    fromkb = str(fromb)[:-3]
    tokb = str(tob)[:-3]
    cmd_str = """
    /nfs/team151/software/plink2_18_April_2015/plink \
        --bgen %(IMPUTE_BGEN_FOLDER)s/impute_%(chr)s_interval.bgen \
        --sample %(SAMPLE_FILE)s \
        --ld-window-r2 0 \
        --chr %(chr)s \
        --from-kb %(fromkb)s \
        --to-kb %(tokb)s \
        --r2 inter-chr""" % {'IMPUTE_BGEN_FOLDER': IMPUTE_BGEN_FOLDER, 'chr': chromosome, 'fromkb': fromkb, 'tokb': tokb,
                             'SAMPLE_FILE': SAMPLE_FILE}
    print(cmd_str)
    proc = subprocess.Popen(cmd_str, shell=True, executable='/usr/bin/bash')
    proc.wait()
    cmd_str = "./locuszoom_format_plink_ld.R"
    proc = subprocess.Popen(cmd_str, shell=True, executable='/usr/bin/bash')
    proc.wait()
    cmd_str = """mv ./plink.ld %(OUT_DIR)s/locuszoom/ld/%(chr)s_%(start)s_%(end)s.ld \
                 && mv ./plink.nosex %(OUT_DIR)s/locuszoom/ld/%(chr)s_%(start)s_%(end)s.nosex \
                 && mv ./plink.log %(OUT_DIR)s/locuszoom/ld%(chr)s_%(start)s_%(end)s.log""" % {'OUT_DIR': OUT_DIR,
                                                                                                          'start': fromb, 'end': tob,
                                                                                                          'chr': chromosome}
    proc = subprocess.Popen(cmd_str, shell=True, executable='/usr/bin/bash')
    proc.wait()

def download_region(chromosome, fromb, tob):
    chromosome = str(chromosome)
    fromb = str(fromb)
    tob = str(tob)
    cmd_str = "tabix -f ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr"+chromosome+".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz "+chromosome+":"+fromb+"-"+tob+" > "+OUT_DIR+"/locuszoom/vcf/"+chromosome+"_"+fromb+"_"+tob+".vcf && rm ./*.vcf.gz.tbi && bgzip -c "+OUT_DIR+"/locuszoom/vcf/"+chromosome+"_"+fromb+"_"+tob+".vcf > "+OUT_DIR+"/locuszoom/vcf/"+chromosome+"_"+fromb+"_"+tob+".vcf.gz && tabix -p vcf "+OUT_DIR+"/locuszoom/vcf/"+chromosome+"_"+fromb+"_"+tob+".vcf.gz"
    print(cmd_str)
    proc = subprocess.Popen(cmd_str, shell=True, executable='/usr/bin/bash')
    proc.wait()

def locuszoom_plot(phen_arr, chromosome, start, end, job_name, ld_flag=False, refsnp=None):
    for phen in phen_arr:
        
        # if pdf exists don't regenerate
        import os
        if os.path.exists(OUT_DIR+"/locuszoom/output/"+job_name+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)+"/"+phen+".pdf"):
            print("already done "+phen+".pdf")
            return(True)
        # create folder
        cmd_str = "mkdir -p "+OUT_DIR+"/locuszoom/output/"+job_name+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)
        proc = subprocess.Popen(cmd_str, shell=True, executable='/usr/bin/bash')
        proc.wait()

        refsnp_add = ""
        print(refsnp)
        if refsnp is not None:
            refsnp_add = " --refsnp \""+refsnp+"\""

        fileloc = OUT_DIR+"/locuszoom/input/"+phen+".tsv"
        cmd_str = SUPPORT_FILES+"/locuszoom/bin/locuszoom --metal "+fileloc+" --pop EUR --build hg19 --source 1000G_Nov2014 --chr "+str(chromosome)+" --start "+str(start)+" --end "+str(end)+" --prefix "+OUT_DIR+"/locuszoom/output/"+job_name+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)+"/"+phen+" --plotonly"+refsnp_add
        #cmd_str = cmd_str + " showRefsnpAnnot=FALSE colorCol=colour"

        #--denote-markers-file "+OUT_DIR+"/locuszoom/label_file.tsv"
        if ld_flag != False:
            cmd_str = cmd_str + str(ld_flag)
        print(cmd_str)
        proc = subprocess.Popen(cmd_str, shell=True, executable='/usr/bin/bash')
        proc.wait()

def locuszoom_plot_custom(phen_arr, chromosome, start, end, ld_flag=False):
    for phen in phen_arr:
        # create folder
        cmd_str = "mkdir -p "+OUT_DIR+"/locuszoom/output/"+str(chromosome)+"_"+str(start)+"_"+str(end)
        proc = subprocess.Popen(cmd_str, shell=True, executable='/usr/bin/bash')
        proc.wait()

        fileloc = OUT_DIR+"/locuszoom/input/"+phen+"_col.tsv"
        cmd_str = SUPPORT_FILES+"/locuszoom/bin/locuszoom --metal "+fileloc+" --pop EUR --build hg19 --source 1000G_Nov2014 --chr "+str(chromosome)+" --start "+str(start)+" --end "+str(end)+" --prefix "+OUT_DIR+"/locuszoom/output/"+str(chromosome)+"_"+str(start)+"_"+str(end)+"/"+phen+" --plotonly ymax=350 showGenes=FALSE"
        cmd_str = cmd_str + " showRefsnpAnnot=FALSE colorCol=colour"

        #--denote-markers-file "+OUT_DIR+"/locuszoom/label_file.tsv"
        if ld_flag != False:
            cmd_str = cmd_str + str(ld_flag)
        print(cmd_str)
        proc = subprocess.Popen(cmd_str, shell=True, executable='/usr/bin/bash')
        proc.wait()


monocyte_traits = ['mo_x_ch', 'mo_y_ch', 'mo_z_ch', 'mo_wx', 'mo_wy', 'mo_wz', 'mono', 'mono_p']
lymphocyte_traits = ['ly_x_ch', 'ly_y_ch', 'ly_z_ch', 'ly_wx', 'ly_wy', 'ly_wz', 'lymph', 'lymph_p']

jobs = []
job_dict = {}

job_dict["chromosome"] = 17
job_dict["start"] = int(int(29096406)-200e3)
job_dict["end"] = int(int(29151794)+200e3)
job_dict["phen_arr"] = phen_arr
job_dict["name"]="CRLF3"
job_dict["refsnp"] = "rs6505211"
jobs.append(job_dict.copy())

job_dict["chromosome"] = 6
job_dict["start"] = int(int(36461669)-200e3)
job_dict["end"] = int(int(36515247)+200e3)
job_dict["phen_arr"] = phen_arr
job_dict["name"]="STK38"
job_dict["refsnp"] = "rs141301223"
jobs.append(job_dict.copy())

job_dict["chromosome"] = 6
job_dict["start"] = int(int(74379655)-300e3)
job_dict["end"] = int(int(74406025)+300e3)
job_dict["phen_arr"] = phen_arr
job_dict["name"]="MOB1A"
job_dict["refsnp"] = "rs111623588"
jobs.append(job_dict.copy())



# execute
print(jobs[LSB_JOBINDEX-1]["name"])
locuszoom_plot(jobs[LSB_JOBINDEX-1]["phen_arr"], jobs[LSB_JOBINDEX-1]["chromosome"], jobs[LSB_JOBINDEX-1]["start"], jobs[LSB_JOBINDEX-1]["end"], jobs[LSB_JOBINDEX-1]["name"], refsnp=jobs[LSB_JOBINDEX-1]["refsnp"])

# bsub -q normal -n 1 -J "locuszoomplot[1]" -o ./output/niche/locuszoom_plot_%J_%I.out -e ./output/niche/locuszoom_plot_%J_%I.err -R"select[mem>10000] rusage[mem=10000] span[hosts=1]" -M10000 < /software/python-2.7.10/bin/python2.7 ./locuszoom_script_wrapper.sh

job_dict["chromosome"] = 12
job_dict["start"] = int(47841537)-100000
job_dict["end"] = int(47943048)+100000
job_dict["phen_arr"] = phen_arr
job_dict["name"]="VDR"
jobs.append(job_dict.copy())


#job_dict["chromosome"] = 17
#job_dict["start"] = int(29096406)-300000
#job_dict["end"] = int(29096406)+300000
#job_dict["phen_arr"] = ['pdw']
#job_dict["name"]="CRLF3"
#jobs.append(job_dict.copy())
##

#job_dict["chromosome"] = 14
#job_dict["start"] = int( 21511385)-300000
#job_dict["end"] = int(21511385)+300000
#job_dict["phen_arr"] = ['mo_x_ch', 'mo_y_ch', 'mo_z_ch', 'mo_wx', 'mo_wy', 'mo_wz', 'mono', 'mono_p']
#job_dict["name"]="RNASE7"
#jobs.append(job_dict.copy())


#ob_dict["chromosome"] = 14
#ob_dict["start"] = int(25043724)-300000
#ob_dict["end"] = int(25043724)+300000
#ob_dict["phen_arr"] = ['ne_ssc_ch', 'ne_fsc_ch', 'ne_sfl_ch', 'ne_wx', 'ne_wy', 'ne_wz', 'neut']
#ob_dict["name"]="CTSG"
#obs.append(job_dict.copy())

#job_dict["chromosome"] = 11
#job_dict["start"] = int(88046760)-300000
#job_dict["end"] = int(88046760)+300000
#job_dict["phen_arr"] = ['ne_ssc_ch', 'ne_fsc_ch', 'ne_sfl_ch', 'ne_wx', 'ne_wy', 'ne_wz', 'neut']
#job_dict["name"]="CTSC"
#jobs.append(job_dict.copy())


#job_dict["chromosome"] = 14
#job_dict["start"] = int(21250124)-300000
#job_dict["end"] = int(21250124)+300000
#job_dict["phen_arr"] = ['mo_y_ch', 'mono']
#job_dict["name"]="RNASE6"
#jobs.append(job_dict.copy())


## NBEAL1
#job_dict["chromosome"] = 2
#job_dict["start"] = 203879602-(450000/2)
#job_dict["end"] = 203879602+(450000/2)
#job_dict["phen_arr"] = ['ne_wx', 'ne_wy', 'ne_wz', 'neut', 'ne_ssc_ch', 'ne_sfl_ch', 'ne_fsc_ch']
#job_dict["name"]="NBEAL1"
#jobs.append(job_dict.copy())

# NBEAL2
#job_dict["chromosome"] = 3
#job_dict["start"] = int(47021173-(450000/2))
#job_dict["end"] = int(47021173+(450000/2))
##job_dict["phen_arr"] = ['pdw', 'mpv', 'plt', 'pct']
#job_dict["phen_arr"] = ['PLT_F_Xch', 'PLT_F_Ych', 'PLT_F_Zch', 'PF_PLTF_WX', 'PF_PLTF_WY', 'PF_PLTF_WZ', 'adj_PLT_F_Xch', 'adj_PLT_F_Ych', 'adj_PLT_F_Zch', #'adj_PF_PLTF_WX', 'adj_PF_PLTF_WY', 'adj_PF_PLTF_WZ']
#job_dict["name"]="NBEAL2"
#jobs.append(job_dict.copy())
#
# CRLF3
#job_dict["chromosome"] = 17
#job_dict["start"] = int(29124100-(450000/2))
#job_dict["end"] = int(29124100+(450000/2))
#job_dict["phen_arr"] = ['pdw', 'mpv', 'plt', 'pct']
#job_dict["phen_arr"] = ['PLT_F_Xch', 'PLT_F_Ych', 'PLT_F_Zch', 'PF_PLTF_WX', 'PF_PLTF_WY', 'PF_PLTF_WZ', 'adj_PLT_F_Xch', 'adj_PLT_F_Ych', 'adj_PLT_F_Zch', #'adj_PF_PLTF_WX', 'adj_PF_PLTF_WY', 'adj_PF_PLTF_WZ']
#job_dict["name"]="CRLF3"
#jobs.append(job_dict.copy())


job_dict["chromosome"] = 6
job_dict["start"] = int(147486533-(450000/2))
job_dict["end"] = int(147486533+(450000/2))
job_dict["phen_arr"] = ['PLT_F_Xch']
job_dict["name"]="STXBP5-AS1"
jobs.append(job_dict.copy())

job_dict["chromosome"] = 6
job_dict["start"] = int(31695368-(450000/2))
job_dict["end"] = int(31695368+(450000/2))
job_dict["phen_arr"] = ['h_ipf']
job_dict["name"]="DDAH2"
jobs.append(job_dict.copy())

job_dict["chromosome"] = 6
job_dict["start"] = int(47622573-(450000/2))
job_dict["end"] = int(47622573+(450000/2))
job_dict["phen_arr"] = ['h_ipf']
job_dict["name"]="GPR111"
jobs.append(job_dict.copy())

job_dict["chromosome"] = 9
job_dict["start"] = int(99234329-(450000/2))
job_dict["end"] = int(99234329+(450000/2))
job_dict["phen_arr"] = ['PF_PLTF_WX']
job_dict["name"]="HABP4"
jobs.append(job_dict.copy())

job_dict["chromosome"] = 12
job_dict["start"] = int(29393102-(450000/2))
job_dict["end"] = int(29393102+(450000/2))
job_dict["phen_arr"] = ['h_ipf']
job_dict["name"]="FAR2"
jobs.append(job_dict.copy())

job_dict["chromosome"] = 12
job_dict["start"] = int(29393102-(450000/2))
job_dict["end"] = int(29393102+(450000/2))
job_dict["phen_arr"] = ['h_ipf']
job_dict["name"]="FAR2"
jobs.append(job_dict.copy())

job_dict["chromosome"] = 3
job_dict["start"] = int(151054836-(450000/2))
job_dict["end"] = int(151054836+(450000/2))
job_dict["phen_arr"] = ['adj_PLT_F_Ych']
job_dict["name"]="MED12L"
jobs.append(job_dict.copy())

job_dict["chromosome"] = 2
job_dict["start"] = int(55212337-(450000/2))
job_dict["end"] = int(55212337+(450000/2))
job_dict["phen_arr"] = ['PLT_F_Xch']
job_dict["name"]="RTN4"
jobs.append(job_dict.copy())
