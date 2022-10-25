#!/bin/bash
set -ue

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
conda activate starrpeaker

## Get file paths to required arguments ##
# Prefix to store peaks
prefix=$1

# starrseq input file
input=$2

#starrseq output file
output=$3

# chromosome sizes file
chromsize=$4 

# blacklist regions file
blacklist=$5

# covariate files
cov1=$6
cov2=$7
cov3=$8

# starrpeaker peak calling command
starrpeaker --prefix $prefix --chromsize $chromsize --blacklist $blacklist --cov $cov1 $cov2 $cov3 -i $input -o $output
