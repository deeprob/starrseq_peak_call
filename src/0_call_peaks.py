import argparse
import utils as ut


def call_peaks(
    input_library_filtered_prefix, output_library_filtered_prefix, 
    input_library_replicates, output_library_replicates,
    peak_call_dir, output_library_short,
    reference_genome_twobit, roi_file,
    starrpeaker_data_dir, cradle_data_dir,
    starrpeaker_peak_flag, bigwig_peak_flag, cradle_peak_flag, macs2_peak_flag, replicate_peak_flag
    ):
    # starrpeaker peak call
    if starrpeaker_peak_flag:
        ut.call_starrpeaker_peaks(
            input_library_filtered_prefix, 
            output_library_filtered_prefix,
            starrpeaker_data_dir,
            peak_call_dir, output_library_short,
            )
    # cradle peak call
    if cradle_peak_flag:
        # create bigwig files from filtered bam files, required for cradle
        if bigwig_peak_flag:
            ut.create_bigwig(
                input_library_filtered_prefix, input_library_replicates,
                output_library_filtered_prefix, output_library_replicates
                )
        ut.call_cradle_peaks(
            input_library_filtered_prefix, input_library_replicates,
            output_library_filtered_prefix, output_library_replicates,
            reference_genome_twobit, roi_file, cradle_data_dir,
            peak_call_dir, output_library_short
            )
    # macs2 peak call
    if macs2_peak_flag:
        ut.call_macs2_peaks(
            input_library_filtered_prefix, 
            output_library_filtered_prefix,
            peak_call_dir, 
            output_library_short,
            )
    # starrpeaker replicate peak call
    if replicate_peak_flag:
        ut.call_starrpeaker_peaks_for_each_replicate(
            input_library_filtered_prefix, input_library_replicates, 
            output_library_filtered_prefix, output_library_replicates,
            starrpeaker_data_dir,
            peak_call_dir,
            output_library_short
            )
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq analysis")
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored")
    parser.add_argument("lib", type=str, help="The library name as given in the meta file")
    parser.add_argument("peak_call_dir", type=str, help="The root dir where called peaks will be stored")
    parser.add_argument("starrpeaker_data_dir", type=str, help="Starrpeaker required files data dir")
    parser.add_argument("cradle_data_dir", type=str, help="Cradle required files data dir")
    parser.add_argument("-s", "--starrpeaker", action="store_true", help="Run starrpeaker pipeline")
    parser.add_argument("-b", "--bigwig", action="store_true", help="Create bigwig files")
    parser.add_argument("-c", "--cradle", action="store_true", help="Run cradle pipeline")
    parser.add_argument("-m", "--macs", action="store_true", help="Run macs2 pipeline")
    parser.add_argument("-r", "--rep_peaks", action="store_true", help="Run starrpeaker pipeline and call replicate peaks")
    parser.add_argument("-t", "--threads", type=int, help="number of threads to use", default=64) # TODO: include threads in main argument

    cli_args = parser.parse_args()
    lib_args = ut.create_args(cli_args.meta_file, cli_args.lib)
    call_peaks(
        lib_args.input_library_filtered_prefix,
        lib_args.output_library_filtered_prefix,
        lib_args.input_library_reps,
        lib_args.output_library_reps, 
        cli_args.peak_call_dir,
        lib_args.output_library_short,
        lib_args.reference_genome_twobit,
        lib_args.roi_file,
        cli_args.starrpeaker_data_dir,
        cli_args.cradle_data_dir,
        cli_args.starrpeaker,
        cli_args.bigwig,
        cli_args.cradle,
        cli_args.macs,
        cli_args.rep_peaks,      
        )
