#!/bin/bash
# set -ue # don't set ue because it interrrupts with conda activate

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

# activate conda environment
conda activate starrseq-da


infile=$1
outfile=$2
echo $PWD

Rscript /data5/deepro/starrseq/main_library/3_peak_call/src/shell_scripts/deseq2.R $infile $outfile
