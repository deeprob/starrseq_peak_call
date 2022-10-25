import os
import argparse
import utils as ut


def call_peaks(
    input_library_prefix, output_library_prefix, 
    input_library_replicates, output_library_replicates,
    bam_dir, peak_call_dir, 
    input_library_short, output_library_short,
    reference_genome_twobit, roi_file,
    starrpeaker_meta_data, cradle_meta_data,
    starrpeaker_peak_flag, bigwig_peak_flag, cradle_peak_flag, macs2_peak_flag, deseq2_peak_flag, replicate_peak_flag,
    ):

    # get unique fragments from the merged bam file
    if not roi_file:
        input_bam = os.path.join(bam_dir, input_library_short, f"{input_library_prefix}.bam")
        roi_file = os.path.join(peak_call_dir, output_library_short, "roi_file.bed")
        ut.get_unique_fragments(input_bam, uniq_frag_bed)

    # starrpeaker peak call
    if starrpeaker_peak_flag:
        ut.call_starrpeaker_peaks(
            input_library_prefix, 
            output_library_prefix,
            bam_dir, 
            starrpeaker_meta_data,
            peak_call_dir, 
            input_library_short, output_library_short,
            )

    # starrpeaker replicate peak call
    if replicate_peak_flag:
        ut.call_starrpeaker_peaks_for_each_replicate(
            input_library_prefix, input_library_replicates, 
            output_library_prefix, output_library_replicates,
            bam_dir,
            starrpeaker_meta_data,
            peak_call_dir,
            input_library_short, output_library_short
            )
            
    # cradle peak call
    if cradle_peak_flag:
        # create bigwig files from filtered bam files, required for cradle
        if bigwig_peak_flag:
            ut.create_bigwig(
                input_library_prefix, input_library_replicates,
                output_library_prefix, output_library_replicates,
                bam_dir,
                input_library_short, output_library_short
                )

        ut.call_cradle_peaks(
            input_library_prefix, input_library_replicates,
            output_library_prefix, output_library_replicates,
            bam_dir,
            reference_genome_twobit, roi_file, 
            cradle_meta_data,
            peak_call_dir, 
            input_library_short, output_library_short,
            )

    # macs2 peak call
    if macs2_peak_flag:
        ut.call_macs2_peaks(
            input_library_prefix, output_library_prefix,
            bam_dir, 
            peak_call_dir,
            input_library_short, output_library_short,
            )

    # deseq2 peak call
    if deseq2_peak_flag:
        ut.call_deseq2_peaks(
            input_library_prefix, input_library_replicates,
            output_library_prefix, output_library_replicates,
            roi_file, 
            bam_dir, 
            peak_call_dir, 
            input_library_short, output_library_short,          
        )

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq analysis")
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored")
    parser.add_argument("input_lib", type=str, help="The input library name as given in the meta file")
    parser.add_argument("output_lib", type=str, help="The output library name as given in the meta file")
    parser.add_argument("bam_dir", type=str, help="Directory where the bam files are stored")
    parser.add_argument("peak_call_dir", type=str, help="The root dir where called peaks will be stored")
    parser.add_argument("starrpeaker_meta_data", type=str, help="Starrpeaker required files path as a json file")
    parser.add_argument("cradle_meta_data", type=str, help="Cradle required files path as a json file")
    parser.add_argument("-s", "--starrpeaker", action="store_true", help="Run starrpeaker pipeline")
    parser.add_argument("-b", "--bigwig", action="store_true", help="Create bigwig files")
    parser.add_argument("-c", "--cradle", action="store_true", help="Run cradle pipeline")
    parser.add_argument("-m", "--macs", action="store_true", help="Run macs2 pipeline")
    parser.add_argument("-d", "--deseq", action="store_true", help="Run deseq2 pipeline")
    parser.add_argument("-r", "--rep_peaks", action="store_true", help="Run starrpeaker pipeline and call replicate peaks")
    parser.add_argument("-t", "--threads", type=int, help="number of threads to use", default=64) # TODO: include threads in main argument

    cli_args = parser.parse_args()
    in_lib_args = ut.create_args(cli_args.meta_file, cli_args.input_lib)
    out_lib_args = ut.create_args(cli_args.meta_file, cli_args.output_lib)
    call_peaks(
        in_lib_args.library_prefix,
        out_lib_args.library_prefix,
        in_lib_args.library_reps,
        out_lib_args.library_reps, 
        cli_args.bam_dir,
        cli_args.peak_call_dir,
        in_lib_args.library_short,
        out_lib_args.library_short,
        in_lib_args.reference_genome_twobit,
        in_lib_args.roi_file,
        cli_args.starrpeaker_meta_data,
        cli_args.cradle_meta_data,
        cli_args.starrpeaker,
        cli_args.bigwig,
        cli_args.cradle,
        cli_args.macs,
        cli_args.deseq,
        cli_args.rep_peaks,  
        )
