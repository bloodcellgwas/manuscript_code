#!/bin/bash
#SBATCH -A BRIDGE-CORE-SL2-CPU
#SBATCH --time 36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task 1
#SBATCH --mail-type=FAIL
#SBATCH -p skylake
#SBATCH --output=./output/locuszoom_plot_%A_%a.out
#SBATCH --error=./output/locuszoom_plot_%A_%a.err

module load plink/1.90beta
source /home/pa354/software/miniconda3/bin/activate cond_analysis
export PATH=/software/R-3.3.0/bin:$PATH 
python ${SCRIPT_DIR}/locuszoom_script.py
echo aadsfd > data/locuszoom/output/done.txt
# bsub -q normal -n 1 -J "locuszoomplot[1-5]" -o ./output/niche/locuszoom_plot_%J_%I.out -e ./output/niche/locuszoom_plot_%J_%I.err -R"select[mem>10000] rusage[mem=10000] span[hosts=1]" -M10000 < /software/python-2.7.10/bin/python2.7 ./locuszoom_script_wrapper.sh
#  --array=1-15

# sbatch -J locuszoom --array=2 < locuszoom_script_wrapper.sh

