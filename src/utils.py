import os
import subprocess
import json
from argparse import Namespace

#### GLOBALS ####
CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))

#####################
# argument creation #
#####################

def create_args(meta_file, lib_name):
    with open(meta_file, "r") as f: 
        meta_dict = json.load(f)

    args = Namespace(
        # from metadata file
        input_library_prefix = meta_dict["input"]["prefix"],
        input_library_reps = meta_dict["input"]["replicates"],
        input_library_short = meta_dict["input"]["shortform"],
        output_library_prefix = meta_dict[lib_name]["prefix"],
        output_library_reps = meta_dict[lib_name]["replicates"],
        output_library_short = meta_dict[lib_name]["shortform"],
        reference_genome = meta_dict["genome"]["ref_fasta"],
        reference_genome_twobit = meta_dict["genome"]["ref_twobit"],
        roi_file = meta_dict["roi"]["roi_sorted"]
    )

    args.input_library_filtered_prefix = args.input_library_prefix.replace("raw_data", "filtered")
    args.output_library_filtered_prefix = args.output_library_prefix.replace("raw_data", "filtered")
    return args

###################
# filename parser #
###################

def get_peak_dir_path(peak_call_dir, out_lib_short, rep, peak_call_method):
    """
    Returns the directory where output peak files are stored
    """
    # create the output peaks path where it will be stored
    output_peaks_prefix = os.path.join(peak_call_dir, out_lib_short, rep, peak_call_method)
    # create the output peaks path where it will be stored
    os.makedirs(output_peaks_prefix, exist_ok=True)    
    return output_peaks_prefix

#######################
# starrpeaker helpers #
#######################

# starrpeaker peak calling functions
def call_starrpeaker_peaks_helper(
    peaks_prefix, 
    input_filtered_bam, output_filtered_bam, 
    starrpeaker_data_dir
    ):
    cmd = ["bash", f"{CURRENT_DIR_PATH}/shell_scripts/3a_call_peaks_starrpeaker.sh", 
            f"{peaks_prefix}", 
            f"{input_filtered_bam}", f"{output_filtered_bam}", 
            f"{starrpeaker_data_dir}"]
    subprocess.run(cmd)
    return

def call_starrpeaker_peaks(
    input_library_filtered_prefix, output_library_filtered_prefix,
    starrpeaker_data_dir,
    peak_call_dir,
    output_library_short,
    ):
    """Call peaks for output library using starrpeaker"""
    # create the output peaks path where it will be stored
    output_peaks_prefix = get_peak_dir_path(peak_call_dir, output_library_short, "", "starrpeaker")
    # prepare the input files, starrpeaker requires merged bam
    input_library_filtered_bam =  input_library_filtered_prefix + ".bam"
    output_library_filtered_bam =  output_library_filtered_prefix + ".bam"
    # add peaks to prefix to follow starrpeaker dir creation convention
    output_peaks_prefix = os.path.join(output_peaks_prefix, "peaks")
    call_starrpeaker_peaks_helper(output_peaks_prefix, input_library_filtered_bam, output_library_filtered_bam, starrpeaker_data_dir)
    return

def call_starrpeaker_peaks_for_each_replicate(
    input_library_filtered_prefix, input_library_replicates, 
    output_library_filtered_prefix, output_library_replicates,
    starrpeaker_data_dir, 
    peak_call_dir,
    output_library_short,
    ):
    """Call peaks for each output library replicates using starrpeaker"""
    ireps = input_library_replicates.split()
    oreps = output_library_replicates.split()
    for irep,orep in zip(ireps,oreps):
        # create the output peaks path where it will be stored
        output_peaks_prefix = get_peak_dir_path(peak_call_dir, output_library_short, orep, "starrpeaker")
        # prepare the input files per replicate for starrpeaker
        input_library_filtered_bam =  "_".join([input_library_filtered_prefix, f"{irep}.bam"])
        output_library_filtered_bam =  "_".join([output_library_filtered_prefix, f"{orep}.bam"])
        output_peaks_prefix = os.path.join(output_peaks_prefix, "peaks")
        call_starrpeaker_peaks_helper(output_peaks_prefix, input_library_filtered_bam, output_library_filtered_bam, starrpeaker_data_dir)
    return

##################
# cradle helpers #
##################

# cradle peak calling functions
def create_bigwig_helper(library_prefix, library_replicates):
    cmd = ["bash", f"{CURRENT_DIR_PATH}/shell_scripts/3b_0_bamtobw.sh", "-i", f"{library_prefix}", "-r", f"{library_replicates}"]
    subprocess.run(cmd)
    return

def create_bigwig(
    input_library_prefix, input_library_replicates,
    output_library_prefix, output_library_replicates
    ):
    """
    Create bigwig files from the filtered and filtered bamfiles for visualization and peak calling using cradle
    """
    create_bigwig_helper(input_library_prefix, input_library_replicates)
    create_bigwig_helper(output_library_prefix, output_library_replicates)
    return

def call_cradle_peaks_helper(
    input_bigwigs, output_bigwigs,
    peak_call_dir, 
    roi_file, reference_genome_twobit,
    cradle_data_dir, output_library_short
    ):
    
    cmd = ["bash", f"{CURRENT_DIR_PATH}//shell_scripts/3b_1_call_peaks_cradle.sh"]
    cmd += ["-i", f"{input_bigwigs}", "-o", f"{output_bigwigs}", "-p", peak_call_dir]
    cmd += ["-r", roi_file, "-g", reference_genome_twobit, "-c", cradle_data_dir, "-s", output_library_short]
    subprocess.run(cmd)
    return

def call_cradle_peaks(
    input_library_prefix, input_library_replicates,
    output_library_prefix, output_library_replicates,
    reference_genome_twobit, roi_file, 
    cradle_data_dir,
    peak_call_dir, output_library_short
    ):
    """
    Call peaks for output library using cradle
    """
    output_peaks_prefix = get_peak_dir_path(peak_call_dir, output_library_short, "", "cradle")
    # get the paths where the bigwig files are stores from the library prefixes and replicates
    input_bigwigs = " ".join(["_".join([input_library_prefix, f"{rep}.bw"]) for rep in input_library_replicates.split()])
    output_bigwigs = " ".join(["_".join([output_library_prefix, f"{rep}.bw"]) for rep in output_library_replicates.split()])

    call_cradle_peaks_helper(
        input_bigwigs, output_bigwigs,
        output_peaks_prefix, roi_file, reference_genome_twobit,
        cradle_data_dir, output_library_short
        )
    return

#################
# macs2 helpers #
#################

# macs2 peak calling functions
def call_macs2_peaks_helper(peaks_prefix, input_filtered_bam, output_filtered_bam):
    cmd = ["bash", f"{CURRENT_DIR_PATH}/shell_scripts/3c_call_peaks_macs2.sh", 
            f"{peaks_prefix}", 
            f"{input_filtered_bam}", f"{output_filtered_bam}"]
    subprocess.run(cmd)
    return

def call_macs2_peaks(
    input_library_filtered_prefix, 
    output_library_filtered_prefix,
    peak_call_dir,
    output_library_short,
    ):
    """
    Call peaks for output library using macs2
    """
    # create the output peaks path where it will be stored
    output_peaks_prefix = get_peak_dir_path(peak_call_dir, output_library_short, "", "macs2")
    # prepare the input files, macs2 requires merged bam
    input_library_filtered_bam =  input_library_filtered_prefix + ".bam"
    output_library_filtered_bam =  output_library_filtered_prefix + ".bam"
    # call macs2 peak call function
    call_macs2_peaks_helper(output_peaks_prefix, input_library_filtered_bam, output_library_filtered_bam)
    return
