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
while getopts i:o:p:r:g:c:s: flag
do
    case "${flag}" in
        i) INPUT_BWS=(${OPTARG});;
        o) OUTPUT_BWS=(${OPTARG});;
        p) PEAKDIR=${OPTARG};;
        r) ROIFILE=${OPTARG};;
        g) GENOMEFILE=${OPTARG};;
        c) CRADLE_DATADIR=${OPTARG};;
        s) OUT_SHORTFORM=${OPTARG};;
    esac
done

# blacklist regions file
blacklist="${CRADLE_DATADIR}ENCODE_blacklist_GRCh38_ENCFF419RSJ_merged.bed"

# mappability bias file
mapfile="${CRADLE_DATADIR}mappability/hg38_mappability_100mer.bw" 
# gquad bias file
gquadfile="${CRADLE_DATADIR}gquadruplex/GSE63874_Na_K_PDS_plus_hits_intersect_hg38_uniq_K.bw"

# covariate directory stored file
cov_dir="${CRADLE_DATADIR}hg38_fragLen500_kmer100"

# cradle correctBias command
cradle correctBias -ctrlbw ${INPUT_BWS[@]} -expbw ${OUTPUT_BWS[@]} -l 500 -r ${ROIFILE} -biasType shear pcr map gquad -genome ${GENOMEFILE} -bl ${blacklist} -p 64 -o ${PEAKDIR} -kmer 100 -mapFile $mapfile -gquadFile $gquadfile

# cradle correctBias_stored command
# cradle correctBias_stored -ctrlbw ${input1}.bw ${input2}.bw ${input3}.bw -expbw ${output1}.bw ${output2}.bw ${output3}.bw -r ${ROIFILE} -biasType shear pcr map gquad -genome ${GENOMEFILE} -bl $blacklist -p 64 -o $outprefix_corrected -covariDir $cov_dir

input_corrected_bws=($(find $PEAKDIR -name "Input*corrected.bw" | sort)) # assumes input file prefix starts with "Input"
output_corrected_bws=($(find $PEAKDIR -name "${OUT_SHORTFORM}*corrected.bw" | sort))

echo ${input_corrected_bws[@]}
# # cradle peakcall command
cradle callPeak -ctrlbw ${input_corrected_bws[@]} -expbw ${output_corrected_bws[@]} -r $ROIFILE -fdr 0.05 -o $PEAKDIR -p 64


## command 
# bash 3b_1_call_peaks_cradle.sh -i "/data5/deepro/starrseq/main_library/1_dedup_align_filter/data/filtered/IN/Input_SeqReady_A1_R1_S1.bw /data5/deepro/starrseq/main_library/1_dedup_align_filter/data/filtered/IN/Input_SeqReady_B1_R2_S2.bw /data5/deepro/starrseq/main_library/1_dedup_align_filter/data/filtered/IN/Input_SeqReady_C1_R3_S3.bw" -o "/data5/deepro/starrseq/main_library/1_dedup_align_filter/data/filtered/CC/CC_R1_STARR_Seq_S4.bw /data5/deepro/starrseq/main_library/1_dedup_align_filter/data/filtered/CC/CC_R2_STARR_Seq_S5.bw /data5/deepro/starrseq/main_library/1_dedup_align_filter/data/filtered/CC/CC_R3_STARR_Seq_S6.bw" -p /data5/deepro/starrseq/main_library/3_peak_call/data/CC/cradle -r /data5/deepro/starrseq/main_library/0_region_selection/data/master/master.bed -g /data5/deepro/genomes/hg38/hg38.2bit -c /data5/deepro/starrseq/main_library/3_peak_call/data/cradle_data/ -s CC     
