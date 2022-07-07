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
conda activate cradle

# get all arguments
while getopts i:r: flag
do
    case "${flag}" in
        i) INPUT_PREFIX=${OPTARG};;
        r) BIOL_REPS=(${OPTARG});;
    esac
done

# example command: 
# bash bamtobw.sh -i /data5/deepro/starrseq/main_library/1_dedup_align_filter/data/filtered/IN/Input_SeqReady -r A1_R1_S1 B1_R2_S2 C1_R3_S3

for rep in "${BIOL_REPS[@]}" # ${AR[0]}
do
    bamfile=${INPUT_PREFIX}_${rep}.bam
    outfile=${INPUT_PREFIX}_${rep}.bw

    samtools index $bamfile
    bamCoverage -b $bamfile -o $outfile -of bigwig -p 64
done
