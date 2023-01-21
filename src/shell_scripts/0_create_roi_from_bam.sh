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
conda activate starrseq

IN_BAM=$1
OUT_BED=$2

# convert bam to bed
bedtools bamtobed -i ${IN_BAM} > ${OUT_BED}
# keep only first three columns
cut -f 1,3 -d \t ${OUT_BED} > ${OUT_BED}_tmp
# sort the bed file
bedtools sort -i ${OUT_BED}_tmp > ${OUT_BED}

mv ${OUT_BED}_tmp ${OUT_BED}
# sort the bed file
