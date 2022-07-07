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

conda activate macs2

## Get file paths to required arguments ##
# Prefix to store peaks
peaks_dir=$1

# starrseq input file
input=$2

#starrseq output file
output=$3


macs2 callpeak -t ${output} \
	           -c ${input} \
 	           -f BAMPE \
	           --outdir ${peaks_dir}

