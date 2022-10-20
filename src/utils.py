import os
import subprocess
import json
from argparse import Namespace
import pybedtools


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
        roi_file = meta_dict["roi"]["sorted"]
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


#################
# deseq2 helpers #
#################

def make_windows(in_bed, out_bed, window_size=500, window_stride=50):
    """
    Break the ROIs into fragments of an user defined window size and stride
    """
    window = pybedtools.BedTool().window_maker(b=in_bed ,w=window_size, s=window_stride)
    window_df = window.to_dataframe()
    # get rid of windows which have the same end point
    last_end = None
    rows_to_omit = []
    for i, row in enumerate(window_df.itertuples()):
        if row.end == last_end:
            rows_to_omit.append(i)
        last_end = row.end
    window_df = window_df.loc[~window_df.index.isin(rows_to_omit)]
    os.makedirs(os.path.dirname(out_bed), exist_ok=True)
    window_df.to_csv(out_bed, sep="\t", header=None, index=None)
    return

def get_roi_coverage(filtered_bam, roi_sorted_bed, bed_out):
    """
    Calculates the coverage of each region in the roi file by looking at the bam file
    """
    bam = pybedtools.BedTool(filtered_bam)
    roi = pybedtools.BedTool(roi_sorted_bed)
    c = roi.coverage(bam)
    os.makedirs(os.path.dirname(bed_out), exist_ok=True)
    c.moveto(bed_out)
    return

def get_deseq_compatible_files_from_bed(in_cov_files, out_cov_files, deseq_file):
    # read the bed files
    in_df = pd.concat([pybedtools.BedTool(icf).to_dataframe(disable_auto_names=True, header=None).iloc[:, [0,1,2,-4]].set_index([0,1,2]) for icf in in_cov_files], axis=1)
    out_df = pd.concat([pybedtools.BedTool(icf).to_dataframe(disable_auto_names=True, header=None).iloc[:, [0,1,2,-4]].set_index([0,1,2]) for ocf in out_cov_files], axis=1)
    in_df.columns = [f"Input_Rep_{i}" for i in range(1, len(in_bed_files) + 1)]
    out_df.columns = [f"Output_Rep_{i}" for i in range(1, len(in_bed_files) + 1)]
    df = pd.concat([in_df, out_df], axis=1)
    df.index = ["_".join(list(map(str, i))) for i in df.index]
    df.index.rename("unique_id", inplace=True)
    df.to_csv(deseq_file, index=True, header=True)
    return

def convert_deseq_file_to_bed(deseq_outfile, bed_outdir):
    df = pd.read_csv(deseq_out)
    df["name"] = df.index
    df = df.reset_index(drop=True)
    df = df.merge(df.name.str.split("_", expand=True).rename(columns={0: "chr", 1: "start", 2: "end"}), left_index=True, right_index=True)
    df["strand"] = "."
    df = df.loc[:, ["chr", "start", "end", "name", "strand", "stat", "log2FoldChange", "baseMean", "lfcSE", "pvalue", "padj"]]
    df_active_peaks = df.loc[(df["log2FoldChange"]>0)&(df["padj"]<0.05)]
    df_active_peaks_file = os.path.join(bed_outdir, "peaks.bed")
    df_active_peaks.to_csv(df_active_peaks_file, index=False, header=False, sep="\t")
    df_active_peaks_file_merged = os.path.join(bed_outdir, "peaks.merged.bed")
    df_active_peaks_bed = pybedtools.BedTool.from_dataframe(df_active_peaks)
    df_active_peaks_merged_bed = df_active_peaks_bed.merge(c=[7, 10, 11], o=["max", "min", "min"])
    df_active_peaks_merged_bed.moveto(df_active_peaks_file_merged)
    return


# deseq2 peak calling functions
def call_deseq2_peaks_helper(infile, outfile):
    cmd = ["bash", f"{CURRENT_DIR_PATH}/shell_scripts/1d_call_peaks_deseq2.sh", 
            f"{infile}", 
            f"{outfile}"]
    subprocess.run(cmd)
    return

def call_deseq2_peaks(
    input_library_prefix, input_library_replicates,
    output_library_prefix, output_library_replicates,
    roi_file, bam_dir, peak_call_dir, 
    input_library_short, output_library_short
    ):
    """
    Call peaks for output library using deseq2
    """
    # create the output peaks path where it will be stored
    output_peaks_prefix = get_peak_dir_path(peak_call_dir, output_library_short, "", "deseq2")
    # break the roi file into windows
    roi_window_file = os.path.join(output_peaks_prefix, "roi_windows.bed")
    make_windows(roi_file, roi_window_file)
    # calculate coverage of the roi regions 
    input_library_filtered_bams =  [os.path.join(bam_dir, input_library_short, f"{input_library_prefix}_{rep}.bam") for rep in input_library_replicates.split()]
    output_library_filtered_bams =  [os.path.join(bam_dir, output_library_short, f"{output_library_prefix}_{rep}.bam") for rep in output_library_replicates.split()]
    input_roi_cov_files = [os.path.join(output_peaks_prefix, f"{input_library_prefix}_{rep}.bed") for rep in input_library_replicates.split()]
    output_roi_cov_files = [os.path.join(output_peaks_prefix, f"{output_library_prefix}_{rep}.bed") for rep in output_library_replicates.split()]
    cov_iter = [(ib,ic) for ib,ic in zip(input_library_filtered_bams, roi_window_file, input_roi_cov_files)] + [(ob,oc) for ob,oc in zip(output_library_filtered_bams, roi_window_file, output_roi_cov_files)]
    run_multiargs_pool_job(get_roi_coverage, cov_iter)
    # create deseq compatible file
    deseq_infile = os.path.join(output_peaks_prefix, "deseq_in.csv")
    get_deseq_compatible_files_from_bed(input_roi_cov_files, output_roi_cov_files, deseq_in)
    # call deseq2 peak call function
    deseq_outfile = os.path.join(output_peaks_prefix, "deseq_out.csv")
    call_deseq2_peaks_helper(deseq_infile, deseq_outfile)
    # convert deseq output to compatible bed file
    return


################
# multiprocess #
################

def run_singleargs_pool_job(pool_function, pool_iter):
    pool = mp.Pool(len(pool_iter))
    results = pool.map(pool_function, pool_iter)
    pool.close()
    pool.join()
    return results

def run_multiargs_pool_job(pool_function, pool_iter):
    pool = mp.Pool(len(pool_iter))
    results = pool.starmap(pool_function, pool_iter)
    pool.close()
    pool.join()
    return results
