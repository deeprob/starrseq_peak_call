#!/bin/bash
#SBATCH --account=girirajan
#SBATCH --partition=girirajan
#SBATCH --job-name=starr_peaks
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=400:0:0
#SBATCH --mem-per-cpu=10G
#SBATCH --chdir /data5/deepro/starrseq/main_library/3_peak_call/src
#SBATCH -o /data5/deepro/starrseq/main_library/3_peak_call/slurm/logs/out_%a.log
#SBATCH -e /data5/deepro/starrseq/main_library/3_peak_call/slurm/logs/err_%a.log
#SBATCH --array 9

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data5/deepro/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data5/deepro/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/data5/deepro/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/data5/deepro/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate starrseq

echo `date` starting job on $HOSTNAME
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /data5/deepro/starrseq/main_library/3_peak_call/slurm/data/smap.txt)

echo $LINE
python /data5/deepro/starrseq/main_library/3_peak_call/src/0_call_peaks.py $LINE

echo `date` ending job
