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
while getopts i:o:p:r:g:x:y:l:m:q:s: flag
do
    case "${flag}" in
        i) INPUT_BWS=(${OPTARG});;
        o) OUTPUT_BWS=(${OPTARG});;
        p) PEAKDIR=${OPTARG};;
        r) ROIFILE=${OPTARG};;
        g) GENOMEFILE=${OPTARG};;
        x) IN_PREFIX=${OPTARG};;
        y) OUT_PREFIX=${OPTARG};;
        l) BLACKLIST=${OPTARG};;
        m) MAPFILE=${OPTARG};;
        q) GQUADFILE=${OPTARG};;
        s) STORED_COV=${OPTARG};;
    esac
done


# cradle correctBias command
# cradle correctBias -ctrlbw ${INPUT_BWS[@]} -expbw ${OUTPUT_BWS[@]} -l 500 -r ${ROIFILE} -biasType shear pcr map gquad -genome ${GENOMEFILE} -bl ${BLACKLIST} -p 64 -o ${PEAKDIR} -kmer 100 -mapFile $MAPFILE -gquadFile $GQUADFILE

# cradle correctBias_stored command
cradle correctBias_stored -ctrlbw ${INPUT_BWS[@]} -expbw ${OUTPUT_BWS[@]} -r ${ROIFILE} -biasType shear pcr map gquad -genome ${GENOMEFILE} -bl ${BLACKLIST} -p 64 -o ${PEAKDIR} -covariDir ${STORED_COV} -rngSeed 775948695

input_corrected_bws=($(find $PEAKDIR -name "${IN_PREFIX}*corrected.bw" | sort)) 
output_corrected_bws=($(find $PEAKDIR -name "${OUT_PREFIX}*corrected.bw" | sort))

echo ${input_corrected_bws[@]}
# # cradle peakcall command
cradle callPeak -ctrlbw ${input_corrected_bws[@]} -expbw ${output_corrected_bws[@]} -r $ROIFILE -fdr 0.05 -o $PEAKDIR -p 64


## command 
# bash 3b_1_call_peaks_cradle.sh -i "/data5/deepro/starrseq/main_library/1_dedup_align_filter/data/filtered/IN/Input_SeqReady_A1_R1_S1.bw /data5/deepro/starrseq/main_library/1_dedup_align_filter/data/filtered/IN/Input_SeqReady_B1_R2_S2.bw /data5/deepro/starrseq/main_library/1_dedup_align_filter/data/filtered/IN/Input_SeqReady_C1_R3_S3.bw" -o "/data5/deepro/starrseq/main_library/1_dedup_align_filter/data/filtered/CC/CC_R1_STARR_Seq_S4.bw /data5/deepro/starrseq/main_library/1_dedup_align_filter/data/filtered/CC/CC_R2_STARR_Seq_S5.bw /data5/deepro/starrseq/main_library/1_dedup_align_filter/data/filtered/CC/CC_R3_STARR_Seq_S6.bw" -p /data5/deepro/starrseq/main_library/3_peak_call/data/CC/cradle -r /data5/deepro/starrseq/main_library/0_region_selection/data/master/master.bed -g /data5/deepro/genomes/hg38/hg38.2bit -c /data5/deepro/starrseq/main_library/3_peak_call/data/cradle_data/ -s CC     
