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

# starrpeaker preprocessing files dir
starrdir=$4

# chromosome sizes file
chromsize="${starrdir}GRCh38.chrom.sizes.simple.sorted" 

# blacklist regions file
blacklist="${starrdir}ENCODE_blacklist_GRCh38_ENCFF419RSJ_merged.bed"

# covariate files
cov1="${starrdir}covariates/STARRPeaker_cov_GRCh38_gem-mappability-100mer.bw"
cov2="${starrdir}covariates/STARRPeaker_cov_GRCh38_ucsc-gc-5bp.bw"
cov3="${starrdir}covariates/STARRPeaker_cov_GRCh38_linearfold-folding-energy-100bp.bw"

# starrpeaker peak calling command
starrpeaker --prefix $prefix --chromsize $chromsize --blacklist $blacklist --cov $cov1 $cov2 $cov3 -i $input -o $output
