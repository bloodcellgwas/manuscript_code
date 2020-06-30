#!/software/python-2.7.10/bin/python2.7
# give it a range of ids
import os
import sys
import subprocess
import re
import time
import pdb
#import pandas as pd
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

def generate_submit_gen_extract(lsb_jobindex, memory, chromosome):
    if chromosome in ['X', 'XY']:
        generate_submit_gen_extract_sex(lsb_jobindex, memory, chromosome)
        return(True)
    if memory == 0: memory = 4
    gen_script = """#!/bin/bash
#SBATCH -A BRIDGE-CORE-SL2-CPU
#SBATCH --time 15:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task {memory}
#SBATCH --mail-type=FAIL
#SBATCH -p skylake
#SBATCH --output=/home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_extract_gen/{chromosome}/block_{LSB_JOBINDEX}_%j.out
#SBATCH --error=/home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_extract_gen/{chromosome}/block_{LSB_JOBINDEX}_%j.err

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
#! Modify the settings below to specify the application's environment, location 
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment
export MODULEPATH=/rds/project/who1000-1/rds-who1000-cbrc/user/cbrcmod/modules/out/modulefiles:$MODULEPATH 
module load qctool/beta

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS={memory}

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.

mkdir -p /home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_extract_gen/{chromosome}
source /home/pa354/2018_10_11_ukbb500k_2/envars
export BLOCK_CHR_ARR=($(awk 'NR>1 {{print $2}}' ${{OUT_DIR}}/blocks/block_chrs_order.tsv))
export BLOCK_ARR=($(awk 'NR>1 {{print $1}}' ${{OUT_DIR}}/blocks/block_chrs_order.tsv))

export CHR=${{BLOCK_CHR_ARR[{LSB_JOBINDEX}-1]}}
export BLOCK=${{BLOCK_ARR[{LSB_JOBINDEX}-1]}}
wc -l ${{OUT_DIR}}/blocks/pull_ids_block_${{BLOCK}}.tsv
qctool_v2.0-rc9 \\
    -g ${{OUT_DIR}}/genfiles/chr${{CHR}}_gwsig_clean_id.gen \\
    -s ${{SAMPLE_FILE}} \\
    -incl-variants ${{OUT_DIR}}/blocks/pull_ids_block_${{BLOCK}}.tsv  \\
    -threads {memory} \\
    -og ${{OUT_DIR}}/tmp_gen/block_${{BLOCK}}.gen
if [ $? -ne 0 ]; then {{ echo "QCTOOL Failed, aborting." ; exit 1; }} fi
sstat  -j   $SLURM_JOB_ID.batch   --format=JobID,MaxVMSize
echo qctool done

awk '{{print $1":"$4"_"$5"_"$6}}' ${{OUT_DIR}}/tmp_gen/block_${{BLOCK}}.gen > ${{OUT_DIR}}/tmp_gen/ids_${{BLOCK}}.tsv
if [ $? -ne 0 ]; then {{ echo "AWK Failed, aborting." ; exit 1; }} fi
echo awk done

export TYPE=gen_extract
Rscript /home/pa354/2018_10_11_ukbb500k_2/check_block_subset.R ${{BLOCK}}
if [ $? -ne 0 ]; then {{ echo "check_block_subset.R failed, aborting." ; exit 1; }} fi
sstat  -j   $SLURM_JOB_ID.batch   --format=JobID,MaxVMSize
echo check block done
    """.format(LSB_JOBINDEX=lsb_jobindex, memory=memory, chromosome=chromosome)
    f = open(OUT_DIR+"/block_scripts/"+str(lsb_jobindex)+"_gen_extract.sh", 'w')
    f.write(gen_script)
    f.close()
    directory = "/home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_extract_gen/"+chromosome
    if not os.path.exists(directory):
        os.makedirs(directory)
    #cmd_str = "bsub -q normal -n 1 -J gen_extract_"+str(lsb_jobindex)+" < "+OUT_DIR+"/block_scripts/"+str(lsb_jobindex)+"_gen_extract.sh"
    cmd_str = "sbatch -J gen_extract_"+str(lsb_jobindex) + " < "+OUT_DIR+"/block_scripts/"+str(lsb_jobindex)+"_gen_extract.sh"
    print(cmd_str)
    proc = subprocess.Popen(cmd_str, shell=True, executable='/bin/bash')
    proc.wait()

def generate_submit_gen_extract_sex(lsb_jobindex, memory, chromosome):
    if memory == 0: memory = 1
    gen_script = """#!/bin/bash
#SBATCH -A BRIDGE-CORE-SL2-CPU
#SBATCH --time 07:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task {memory}
#SBATCH --mail-type=FAIL
#SBATCH -p skylake
#SBATCH --output=/home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_extract_gen/{chromosome}/block_{LSB_JOBINDEX}_%j.out
#SBATCH --error=/home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_extract_gen/{chromosome}/block_{LSB_JOBINDEX}_%j.err

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
#! Modify the settings below to specify the application's environment, location 
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment
export MODULEPATH=/rds/project/who1000-1/rds-who1000-cbrc/user/cbrcmod/modules/out/modulefiles:$MODULEPATH 
module load qctool/beta

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS={memory}

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.

mkdir -p /home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_extract_gen/{chromosome}
source /home/pa354/2018_10_11_ukbb500k_2/envars
export BLOCK_CHR_ARR=($(awk 'NR>1 {{print $2}}' ${{OUT_DIR}}/blocks/block_chrs_order.tsv))
export BLOCK_ARR=($(awk 'NR>1 {{print $1}}' ${{OUT_DIR}}/blocks/block_chrs_order.tsv))

export CHR=${{BLOCK_CHR_ARR[{LSB_JOBINDEX}-1]}}
export BLOCK=${{BLOCK_ARR[{LSB_JOBINDEX}-1]}}
wc -l ${{OUT_DIR}}/blocks/pull_ids_block_${{BLOCK}}.tsv
qctool_v2.0-rc9 \\
    -g ${{OUT_DIR}}/genfiles/chrXY_merged.gen \\
    -s ${{OUT_DIR}}/genfiles/chrXY_merged.sample \\
    -incl-variants ${{OUT_DIR}}/blocks/pull_ids_block_${{BLOCK}}.tsv  \\
    -threads {memory} \\
    -og ${{OUT_DIR}}/tmp_gen/block_${{BLOCK}}.gen
if [ $? -ne 0 ]; then {{ echo "QCTOOL Failed, aborting." ; exit 1; }} fi
sstat  -j   $SLURM_JOB_ID.batch   --format=JobID,MaxVMSize
echo qctool done

awk '{{print $1":"$4"_"$5"_"$6}}' ${{OUT_DIR}}/tmp_gen/block_${{BLOCK}}.gen > ${{OUT_DIR}}/tmp_gen/ids_${{BLOCK}}.tsv
if [ $? -ne 0 ]; then {{ echo "AWK Failed, aborting." ; exit 1; }} fi
echo awk done

export TYPE=gen_extract
Rscript /home/pa354/2018_10_11_ukbb500k_2/check_block_subset.R ${{BLOCK}}
if [ $? -ne 0 ]; then {{ echo "check_block_subset.R failed, aborting." ; exit 1; }} fi
sstat  -j   $SLURM_JOB_ID.batch   --format=JobID,MaxVMSize
echo check block done
    """.format(LSB_JOBINDEX=lsb_jobindex, memory=memory, chromosome="XY")
    f = open(OUT_DIR+"/block_scripts/"+str(lsb_jobindex)+"_gen_extract.sh", 'w')
    f.write(gen_script)
    f.close()
    directory = "/home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_extract_gen/"+chromosome
    if not os.path.exists(directory):
        os.makedirs(directory)
    #cmd_str = "bsub -q normal -n 1 -J gen_extract_"+str(lsb_jobindex)+" < "+OUT_DIR+"/block_scripts/"+str(lsb_jobindex)+"_gen_extract.sh"
    cmd_str = "sbatch -J gen_extract_"+str(lsb_jobindex) + " < "+OUT_DIR+"/block_scripts/"+str(lsb_jobindex)+"_gen_extract.sh"
    print(cmd_str)
    proc = subprocess.Popen(cmd_str, shell=True, executable='/bin/bash')
    proc.wait()


def generate_submit_hd5_create(lsb_jobindex, gen_submitted, memory, chromosome):
    if memory == 0: memory = 3
    hd5_script = """#!/bin/bash
#SBATCH --time 01:00:00
#SBATCH -A BRIDGE-CORE-SL2-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task {memory}
#SBATCH --mail-type=FAIL
#SBATCH -p skylake
#SBATCH --error=/home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_create_h5/{chromosome}/block_{LSB_JOBINDEX}_%j.err
#SBATCH --output=/home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_create_h5/{chromosome}/block_{LSB_JOBINDEX}_%j.out

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
#! Modify the settings below to specify the application's environment, location 
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS={memory}

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=comopact.

mkdir -p /home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_create_h5/{chromosome}
source /home/pa354/2018_10_11_ukbb500k_2/envars
export BLOCK_CHR_ARR=($(awk 'NR>1 {{print $2}}' ${{OUT_DIR}}/blocks/block_chrs_order.tsv))
export BLOCK_ARR=($(awk 'NR>1 {{print $1}}' ${{OUT_DIR}}/blocks/block_chrs_order.tsv))
export CHR=${{BLOCK_CHR_ARR[{LSB_JOBINDEX}-1]}}
export BLOCK=${{BLOCK_ARR[{LSB_JOBINDEX}-1]}}

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_2/gen2hd5_lib

/rds/project/who1000-1/rds-who1000-cbrc/projects/gwas_bcx_ukbb500k/conditional_analysis_2/gen2hd5.git/src/gen2hd5 \\
   -g ${{OUT_DIR}}/tmp_gen/block_${{BLOCK}}.gen \\
   -o ${{OUT_DIR}}/tmp_hd5/block_${{BLOCK}}.h5
if [ $? -ne 0 ]; then {{ echo "gen2hd5 Failed, aborting." ; exit 1; }} fi
sstat  -j   $SLURM_JOB_ID.batch   --format=JobID,MaxVMSize
# confirmation that script ran to completion
echo {LSB_JOBINDEX} >> ${{OUT_DIR}}/tmp_hd5/block_${{BLOCK}}.txt
sstat  -j   $SLURM_JOB_ID.batch   --format=JobID,MaxVMSize
""".format(LSB_JOBINDEX=lsb_jobindex, memory=memory, chromosome=chromosome)
    f = open(OUT_DIR+"/block_scripts/"+str(lsb_jobindex)+"_hd5_create.sh", 'w')
    f.write(hd5_script)
    f.close()
    directory = "/home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_create_h5/"+chromosome
    if not os.path.exists(directory):
        os.makedirs(directory)
    if gen_submitted == True:
        cmd_str = "sbatch -J hd5_create_"+str(lsb_jobindex)+" --dependency=afterok:$(squeue --noheader --format %i --name gen_extract_"+str(lsb_jobindex)+") < "+OUT_DIR+"/block_scripts/"+str(lsb_jobindex)+"_hd5_create.sh"
    else :
        cmd_str = "sbatch -J hd5_create_"+str(lsb_jobindex)+" < "+OUT_DIR+"/block_scripts/"+str(lsb_jobindex)+"_hd5_create.sh"
    print(cmd_str)
    proc = subprocess.Popen(cmd_str, shell=True, executable='/bin/bash')
    proc.wait()

def generate_submit_conditional_analysis(lsb_jobindex, hd5_submitted, memory, chromosome, queue="normal"):
    #if memory == 0: memory = 99000 # default memory
    if memory == 0: memory = 15
    conditional_script = """#!/bin/bash
#SBATCH --time 07:00:00
#SBATCH -A BRIDGE-CORE-SL2-CPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task {memory}
#SBATCH --mail-type=FAIL
#SBATCH -p skylake
#SBATCH --error=/home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_run/{chromosome}/block_{LSB_JOBINDEX}_%j.err
#SBATCH --output=/home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_run/{chromosome}/block_{LSB_JOBINDEX}_%j.out
echo CORES USED: {memory}
#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
#! Modify the settings below to specify the application's environment, location 
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel7/default-peta4            # REQUIRED - loads the basic environment

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this
#! safe value to no more than 32:
export OMP_NUM_THREADS={memory}

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=comopact.

mkdir -p /home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_run/{chromosome}
source /home/pa354/2018_10_11_ukbb500k_2/envars
export BLOCK_CHR_ARR=($(awk 'NR>1 {{print $2}}' ${{OUT_DIR}}/blocks/block_chrs_order.tsv))
export BLOCK_ARR=($(awk 'NR>1 {{print $1}}' ${{OUT_DIR}}/blocks/block_chrs_order.tsv))
export CHR=${{BLOCK_CHR_ARR[{LSB_JOBINDEX}-1]}}
export BLOCK=${{BLOCK_ARR[{LSB_JOBINDEX}-1]}}
if [ $? -ne 0 ]; then {{ echo "cond analysis setup script failed" ; exit 1; }} fi
Rscript /home/pa354/2018_10_11_ukbb500k_2/univariate_conditional.R {memory}
if [ $? -ne 0 ]; then {{ echo "cond analysis script failed" ; sstat  -j   $SLURM_JOB_ID.batch   --format=JobID,MaxVMSize; exit 1; }} fi
#rm ${{OUT_DIR}}/tmp_gen/block_${{BLOCK}}.gen
#rm ${{OUT_DIR}}/tmp_hd5/block_${{BLOCK}}.h5
sstat  -j   $SLURM_JOB_ID.batch   --format=JobID,MaxVMSize
""".format(LSB_JOBINDEX=lsb_jobindex, memory=memory, queue=queue, chromosome=chromosome)
    f = open(OUT_DIR+"/block_scripts/"+str(lsb_jobindex)+"_cond_analysis.sh", 'w')
    f.write(conditional_script)
    f.close()
    directory = "/home/pa354/2018_10_11_ukbb500k_2/output/cond_analysis_run/"+chromosome
    if not os.path.exists(directory):
        os.makedirs(directory)
    # if the memory requirement is higher than 250G change the queue to hugemem
    if int(memory) > 196000:
        print("hugemem submission: "+memory)
#        exit()
        queue = "hugemem"
    if hd5_submitted == True:
        cmd_str = "sbatch -J cond_"+str(lsb_jobindex)+" --dependency=afterok:$(squeue --noheader --format %i --name hd5_create_"+str(lsb_jobindex)+") < "+OUT_DIR+"/block_scripts/"+str(lsb_jobindex)+"_cond_analysis.sh"
    else :
        cmd_str = "sbatch -J cond_"+str(lsb_jobindex)+" < "+OUT_DIR+"/block_scripts/"+str(lsb_jobindex)+"_cond_analysis.sh"
    print(cmd_str)
    proc = subprocess.Popen(cmd_str, shell=True, executable='/bin/bash')
    proc.wait()

def entire_pipeline_m(lsb_jobindex):
    if check_chromosome_exists(lsb_jobindex, blocks_array_chr) == False: return(True)
    block_id = str(blocks_array[lsb_jobindex])
    print("block: "+block_id+" lsb jobindex: "+str(lsb_jobindex)+" chromosome "+blocks_array_chr[lsb_jobindex])
    block_num_lines = sum(1 for line in open(OUT_DIR+"/blocks/pull_ids_block_"+block_id+".tsv")) - 1 # -1 because this file has header
    
    # if we have already done conditional analysis on this block then make sure we loaded in the right number of
    # variants and if we did skip this iteration
    if os.path.exists (OUT_DIR+"/condout/blocks/ind_var_block"+block_id+".tsv") == True:
        cond_num_lines = sum(1 for line in open(OUT_DIR+"/condout/blocks_done/block_"+block_id+".tsv"))
        if cond_num_lines != block_num_lines:
            sys.exit("lsb "+str(lsb_jobindex)+" or block:"+str(block_id)+" cond file analysed "+str(cond_num_lines)+" block has: "+str(block_num_lines))
        else:
            delete_log_files(lsb_jobindex, True, True, True)
            return(True)
    # if hd5 file exists than submit conditional analysis job and skip rest of loop because we dont need to extract gen or hd5
    if os.path.exists(OUT_DIR+"/tmp_hd5/block_"+block_id+".txt") == True:
        run_job_status_wrapper(lsb_jobindex, "cond", False)
        return(True)
    # if gen file doesn't already exist submit job
    if os.path.exists(OUT_DIR+"/tmp_gen/ids_"+block_id+".tsv") == False:
        # check if gen extraction failed previously and if so was it due to lack of memory?
        run_job_status_wrapper(lsb_jobindex, "gen", False)
        gen_submitted = True
    else:
        # if gen file already exists check the output
        gen_submitted = False
        # check gen file output
        gen_num_lines = sum(1 for line in open(OUT_DIR+"/tmp_gen/ids_"+block_id+".tsv"))
        if gen_num_lines != block_num_lines:
            print("lsb "+str(lsb_jobindex)+" or block:"+str(block_id)+" gen file has "+str(gen_num_lines)+" block has: "+str(block_num_lines))
            os.remove(OUT_DIR+"/tmp_gen/ids_"+block_id+".tsv")
            os.remove(OUT_DIR+"/tmp_gen/block_"+block_id+".gen")
            return(False)

    # we are going to run a hd5 generation script so make sure the .gen file for this block exists
    if os.path.exists(OUT_DIR+"/tmp_gen/block_"+block_id+".gen") == False and gen_submitted == False:
        sys.exit("gen id file exists, but not the .gen file and hd5 hasn't been run")
    # check if hd5 extraction failed and if so was it due to lack of memory?
    run_job_status_wrapper(lsb_jobindex, "hd5", gen_submitted)
    # if script has gotten to this point it means conditional analysis hasn't already done so submit
    run_job_status_wrapper(lsb_jobindex, "cond", True)

def run_job_status_wrapper(lsb_jobindex, job_type, wait_previous_job):
    job_status = get_job_status(lsb_jobindex, job_type)
    if job_status == "run":
        # job hasn't been run before run with default memory
        if job_type == "cond": generate_submit_conditional_analysis(lsb_jobindex, wait_previous_job, 0, get_chromosome(lsb_jobindex))
        elif job_type == "gen": generate_submit_gen_extract(lsb_jobindex, 0, get_chromosome(lsb_jobindex))
        elif job_type == "hd5": generate_submit_hd5_create(lsb_jobindex, wait_previous_job, 0, get_chromosome(lsb_jobindex))
    elif job_status == True:
        sys.exit("according to job status output file lsb jobindex: "+str(lsb_jobindex)+" "+job_type+" is apparently finished but file doesn't exist! (output file has now been deleted)")
    elif job_status == "running":
        # job already running so dont do anything
        print("job already running")
        return(True)
    elif isinstance(job_status, tuple):
        if job_type == "cond": generate_submit_conditional_analysis(lsb_jobindex, wait_previous_job, memory = job_status[1], chromosome = get_chromosome(lsb_jobindex), queue="long")
        else: sys.exit("lsb jobindex: "+str(lsb_jobindex) + " job type: "+job_type+" ran out of time")
    elif job_status != False:
        # job hasn't been run before run with default memory
        if job_type == "cond": generate_submit_conditional_analysis(lsb_jobindex, wait_previous_job, job_status, get_chromosome(lsb_jobindex))
        elif job_type == "gen": generate_submit_gen_extract(lsb_jobindex, job_status, get_chromosome(lsb_jobindex))
        elif job_type == "hd5": generate_submit_hd5_create(lsb_jobindex, wait_previous_job, job_status, get_chromosome(lsb_jobindex))

def get_job_status(lsb_jobindex, job_type):
    # check if job is running, if so return running
    if grep_job_bsub(lsb_jobindex, job_type) == True:
        return("running")

    if job_type == "cond":
        log_path = LOG_DIR+"/cond_analysis_run/"+get_chromosome(lsb_jobindex)
    elif job_type == "hd5":
        log_path = LOG_DIR+"/cond_analysis_create_h5/"+get_chromosome(lsb_jobindex)
    elif job_type == "gen":
        log_path = LOG_DIR+"/cond_analysis_extract_gen/"+get_chromosome(lsb_jobindex)
    elif job_type == "test":
        log_path = LOG_DIR+"/test"
   
    if os.path.exists(log_path) == False: return("run")

    # list all files in directory
    output_files = os.listdir(log_path)
    matching = [s for s in output_files if "block_"+str(lsb_jobindex)+"_" in s]
    matching = [s for s in matching if ".out" in s]
    # job doesn't exit, if so return "run"
    if len(matching) == 0: return("run")
    # extract the latest one
    old_versions = []
    latest = 0
    for output in matching:
        m = re.match("block_"+str(lsb_jobindex)+"_(\d+).out", output)
        if latest == 0:
            latest = int(m.group(1))
        elif latest < int(m.group(1)):
            old_versions.append(str(latest))
            latest = int(m.group(1))
        else:
            old_versions.append(str(m.group(1)))
            
    # delete old versions
    for version in old_versions:
        os.remove(log_path+"/block_"+str(lsb_jobindex)+"_"+version+".out")
        os.remove(log_path+"/block_"+str(lsb_jobindex)+"_"+version+".err")

    # read output file TODO remove this when fixed
#    if os.path.exists(log_path+"/block_"+str(lsb_jobindex)+"_"+str(latest)+".out"):
#        return(False)
    output_file = open(log_path+"/block_"+str(lsb_jobindex)+"_"+str(latest)+".out", 'r')
    output_file_text = output_file.read()
    output_file.close()

    error_file = open(log_path+"/block_"+str(lsb_jobindex)+"_"+str(latest)+".err", 'r')
    error_file_text = error_file.read()
    error_file.close()

    if "Successfully completed." in output_file_text:
        if job_type == "gen": delete_log_files(lsb_jobindex, gen_extract = True)
        elif job_type == "hd5": delete_log_files(lsb_jobindex, hd5 = True)
        elif job_type == "cond": delete_log_files(lsb_jobindex, conditional = True)
        return(True) # do we want to delete log files?
    elif "TERM_MEMLIMIT: job killed after reaching LSF memory usage limit." in output_file_text:
        new_memory_limit = str(round(int(re.findall("Max Memory :\s+([\d\.]+) MB", output_file_text)[0])*2.5)).replace(".0", "")
        print(str(lsb_jobindex)+" new mem limit: "+new_memory_limit)
        return(new_memory_limit)
    elif ("Cannot allocate memory" in error_file_text or "std::bad_alloc" in error_file_text) and 1 == 0:
        new_memory_limit = str(round(int(re.findall("Total Requested Memory :\s+([\d\.]+) MB", output_file_text)[0].replace(".00", ""))*2.5)).replace(".0", "")
        print(str(lsb_jobindex)+" new mem limit: "+new_memory_limit)
        return(new_memory_limit)
    elif "TERM_RUNLIMIT" in output_file_text:
        new_memory_limit = re.findall("Total Requested Memory :\s+([\d\.]+) MB", output_file_text)[0].replace(".00", "")
        print(new_memory_limit)
        return(("time", new_memory_limit))
    elif "slurmstepd: error: Exceeded step memory limit at some point." in error_file_text or "Cannot allocate memory" in error_file_text or "add more memory" in error_file_text:
        # calculate how many cpus the job had before
        new_memory_limit = slurm_job_cpus(lsb_jobindex, job_type) + 8
        print(str(new_memory_limit) + "----------")
        if new_memory_limit > 32: return(32)
        return(new_memory_limit)
    else:
        error_file = open(log_path+"/block_"+str(lsb_jobindex)+"_"+str(latest)+".err", 'r')
        print(str(lsb_jobindex) + " " + job_type + " failed\n"+error_file.read())
        return(False)

def slurm_job_cpus(lsb_jobindex, job_type):
    if job_type == "gen": postfix = "gen_extract"
    elif job_type == "hd5": postfix = "hd5_create"
    elif job_type == "cond": postfix = "cond_analysis"

    script_file = open(OUT_DIR+"/block_scripts/"+str(lsb_jobindex)+"_"+postfix+".sh", 'r')
    script_file_text = script_file.read()
    script_file.close()
    cpus = int(re.findall("#SBATCH --cpus-per-task (\d+)", script_file_text)[0])
    return(cpus)

def grep_job_bsub(lsb_jobindex, job_type):
    if job_type == "cond": job_name = "cond_"+str(lsb_jobindex)
    if job_type == "hd5": job_name = "hd5_create_"+str(lsb_jobindex)
    if job_type == "gen": job_name = "gen_extract_"+str(lsb_jobindex)
    cmd_str = ["squeue", "-u", "pa354", "--format=\"%.20j\""]
    proc = subprocess.Popen(cmd_str, stdout=subprocess.PIPE)
    proc.wait()
    output = str(proc.stdout.read())
    cmd_str = ["squeue", "-u", "wja24", "--format=\"%.20j\""]
    proc = subprocess.Popen(cmd_str, stdout=subprocess.PIPE)
    proc.wait()
    output += "\n"+str(proc.stdout.read())
    cmd_str = ["squeue", "-u", "tj241", "--format=\"%.20j\""]
    proc = subprocess.Popen(cmd_str, stdout=subprocess.PIPE)
    proc.wait()
#    output += "\n"+str(proc.stdout.read())
#    cmd_str = ["squeue", "-u", "asb38", "--format=\"%.18i %.10P %.20j %.8u %.8T %.10M %.9l %.6D %R\""]
#    proc = subprocess.Popen(cmd_str, stdout=subprocess.PIPE)
#    proc.wait()
#    output += "\n"+str(proc.stdout.read())
#    cmd_str = ["squeue", "-u", "dg333", "--format=\"%.18i %.10P %.20j %.8u %.8T %.10M %.9l %.6D %R\""]
#    proc = subprocess.Popen(cmd_str, stdout=subprocess.PIPE)
#    proc.wait()
#    output += "\n"+str(proc.stdout.read())
#    cmd_str = ["squeue", "-u", "et341", "--format=\"%.18i %.10P %.20j %.8u %.8T %.10M %.9l %.6D %R\""]
#    proc = subprocess.Popen(cmd_str, stdout=subprocess.PIPE)
#    proc.wait()
#    output += "\n"+str(proc.stdout.read())

    if job_name in output:
        return(True)
    else:
        return(False)

def delete_log_files(lsb_jobindex, conditional = False, hd5 = False, gen_extract = False):
    if conditional == True:
        cmd_str = "rm -f "+LOG_DIR+"/cond_analysis_run/"+get_chromosome(lsb_jobindex)+"/block_"+str(lsb_jobindex)+"*"
        proc = subprocess.Popen(cmd_str, shell=True, executable='/bin/bash')
        proc.wait()
    if hd5 == True:
        cmd_str = "rm -f "+LOG_DIR+"/cond_analysis_create_h5/"+get_chromosome(lsb_jobindex)+"/block_"+str(lsb_jobindex)+"*"
        proc = subprocess.Popen(cmd_str, shell=True, executable='/bin/bash')
        proc.wait()
    if gen_extract == True:
        cmd_str = "rm -f "+LOG_DIR+"/cond_analysis_extract_gen/"+get_chromosome(lsb_jobindex)+"/block_"+str(lsb_jobindex)+"*"
        proc = subprocess.Popen(cmd_str, shell=True, executable='/bin/bash')
        proc.wait()

def rerun_gen_script(lsb_jobindex, block_id):
    # read the existing log file
    log_path = LOG_DIR+"/cond_analysis_extract_gen"
    output_files = os.listdir(log_path)
    matching = [s for s in output_files if "block_"+str(lsb_jobindex) in s]
    matching = [s for s in matching if ".out" in s]

    # job doesn't exit, if so return "run"
    if len(matching) == 0: return("run")

    # extract the latest one
    old_versions = []
    latest = 0
    for output in matching:
        m = re.match("block_"+str(lsb_jobindex)+"_(\d+).out", output)
        if latest == 0:
            latest = int(m.group(1))
        elif latest < int(m.group(1)):
            old_versions.append(str(latest))
            latest = int(m.group(1))
        else:
            old_versions.append(str(m.group(1)))
            
    # delete old versions
    #for version in old_versions:
    #    os.remove(log_path+"/block_"+str(lsb_jobindex)+"_"+version+".out")
    #    os.remove(log_path+"/block_"+str(lsb_jobindex)+"_"+version+".err")

    # read output file
    output_file = open(log_path+"/block_"+str(lsb_jobindex)+"_"+str(latest)+".err", 'r')
    output_file_text = output_file.read()
    output_file.close()
    print(output_file_text)
    os.remove(log_path+"/block_"+str(lsb_jobindex)+"_"+str(latest)+".out")

    # delete output files
    delete_log_files(lsb_jobindex, conditional = False, hd5 = False, gen_extract = True)

    # delete the ids file
    os.remove(OUT_DIR+"/tmp_gen/ids_"+str(block_id)+".tsv")

    #time.sleep(5)

chromosomes = []
for chromosome in range(1,25):
    chromosomes.append([])
    
for lsb_jobindex in range(1, len(blocks_array_chr_numeric)):
    chromosomes[int(get_chromosome_numeric(lsb_jobindex))-1].append(lsb_jobindex)

def del_if_exists(fpath):
    if os.path.isfile(fpath):
        os.remove(fpath)

# get data frame of block sizes
#block_sizes_db = pd.read_csv(OUT_DIR+"/blocks/block_sizes_db.csv")
#block_sizes_db = block_sizes_db.loc[block_sizes_db["block_size"] < 200,:]
#ind = 1
#for index, row in block_sizes_db.iterrows():
#    ind += 1
#    if ind < 700: print("next"); next
#    if row["chr_num"] == 23: next
#    entire_pipeline_m(row["lsb_jobindex"])


#chromosome = 23
#f = open(OUT_DIR+"/empty_blocks.csv", 'r')
#missing_blocks = []
#for line in f:
#    row = line.split(" ")
#    missing_blocks.append(int(row[0].strip()))

#import pdb
#pdb.set_trace()
#for lsb_jobindex in [item for sublist in chromosomes for item in sublist]:
for chromosome in range(23,24):    
    for lsb_jobindex in chromosomes[chromosome-1]:
        entire_pipeline_m(lsb_jobindex)
#import pdb; pdb.set_trace()
#entire_pipeline_m(4702)

def check_chromosome_done(chromosome, chromosomes):
    remaining_blocks = 0
    if check_chromosome_exists(chromosomes[chromosome-1][0], blocks_array_chr) == False:
        print("Chromosome gen file "+str(chromosome)+" doesn't exist")
    for lsb_jobindex in chromosomes[chromosome-1]:
        block_id = str(blocks_array[lsb_jobindex])
        block_num_lines = sum(1 for line in open(OUT_DIR+"/blocks/pull_ids_block_"+block_id+".tsv")) - 1
        # -1 because this file has header
        if os.path.exists(OUT_DIR+"/condout/blocks/ind_var_block"+block_id+".tsv") == True:
            if os.path.exists(OUT_DIR+"/condout/blocks_done/block_"+block_id+".tsv") == False:
                os.remove(OUT_DIR+"/condout/blocks/ind_var_block"+block_id+".tsv")
                continue
            cond_num_lines = sum(1 for line in open(OUT_DIR+"/condout/blocks_done/block_"+block_id+".tsv"))
            if cond_num_lines != block_num_lines:
                sys.exit("(lsb "+str(lsb_jobindex)+" or block: "+ str(block_id)+") cond file analysed has "+
                         str(cond_num_lines)+" variants, but block has: "+str(block_num_lines)) + " variants"
        else:
#            print("\t "+str(lsb_jobindex))
            remaining_blocks = remaining_blocks + 1

    if remaining_blocks == 0:
        print("Chromosome "+str(chromosome)+" is COMPLETE!")
        return(True)
    else:
         print("Chromosome "+str(chromosome)+" isn't complete, " + str(remaining_blocks)+ " blocks haven't been analysed!")
         return(False)

for chromosome in range(1,24):
    check_chromosome_done(chromosome, chromosomes)
