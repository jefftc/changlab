#!/usr/bin/env python

def make_peakseq_preproc_command(bam_file, out_path):
    from genomicode import config
    from genomicode import filelib
    from genomicode import shell

    # samtools view bam11.bam | PeakSeq -preprocess SAM stdin bam12
    samtools = filelib.which_assert(config.samtools)
    peakseq = filelib.which_assert(config.peakseq)
    sq = shell.quote
    cmd = [
        sq(samtools),
        "view", sq(bam_file),
        "|",
        sq(peakseq),
        "-preprocess",
        "SAM",
        "stdin",
        sq(out_path),
        ]
    return " ".join(cmd)

def make_peakseq_run_command(config_file):
    import os
    from genomicode import config
    from genomicode import filelib
    from genomicode import shell

    assert os.path.exists(config_file)
    config_file = os.path.realpath(config_file)

    # PeakSeq -peak_select <config_file>
    peakseq = filelib.which_assert(config.peakseq)
    sq = shell.quote
    cmd = [
        sq(peakseq),
        "-peak_select", config_file,
        ]
    return " ".join(cmd)

def make_config_file(
    filename, treatment_path, control_path, mapability_file,
    experiment_name=None, fragment_length=None):
    import os

    experiment_name = experiment_name or "experiment"
    fragment_length = fragment_length or 146

    assert os.path.exists(treatment_path)
    assert os.path.isdir(treatment_path)
    assert os.path.exists(control_path)
    assert os.path.isdir(control_path)
    assert os.path.exists(mapability_file)
    assert fragment_length > 0 and fragment_length < 10000

    # Make full paths in case running directory changes.
    treatment_path = os.path.realpath(treatment_path)
    control_path = os.path.realpath(control_path)
    mapability_file = os.path.realpath(mapability_file)
    
    
    handle = open(filename, 'w')
    w = handle.write

    w("Experiment_id %s\n" % experiment_name)
    
    # Enrichment fragment length For tag extension, this is the value
    # of average fragment length.
    w("Enrichment_mapped_fragment_length %d\n" % fragment_length)

    # Target FDR in the simulations.
    w("target_FDR 0.05\n")

    # Number of simulations performed while estimating the putative
    # peaks.
    w("N_Simulations 10\n")

    # Minimum distance between consecutive peaks
    w("Minimum_interpeak_distance 200\n")

    # Mappability file that includes the uniquely mappable number of
    # nucleotides per window for each chromosome.
    w("Mappability_map_file %s\n" % mapability_file)

    # The directory that contains the preprocessed ChIP-Seq reads, can
    # specify multiple directories to pool reads from multiple source
    # (e.g. replicates)
    w("ChIP_Seq_reads_data_dirs %s\n" % treatment_path)

    # The directory that contains the preprocessed Input (control)
    # experiment reads. (Multiple directories allowed)
    w("Input_reads_data_dirs %s\n" % control_path)

    # Seed for pseudo-random number generator. This is necessary for
    # simulated background option (specified below).
    #Simulation_seed 1234567
    
    w("max_Qvalue 0.05\n")

    # Background_model Poisson
    w("Background_model Simulated\n")


def main():
    import os
    import shutil
    import argparse
    from genomicode import filelib
    from genomicode import shell

    p = filelib.tswrite
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("mapability_file", help="PeakSeq mapability file.")
    parser.add_argument("treatment_bam", help="BAM file of treated sample.")
    parser.add_argument("control_bam", help="BAM file of background sample.")
    parser.add_argument("outpath", help="Directory to store the results.")
    
    parser.add_argument("--experiment_name", help="Name of experiment.")
    parser.add_argument("--fragment_length", type=int, help="")
    #group.add_argument(
    #    "--noclobber", action="store_true",
    #    help="Don't overwrite files if they already exist.")

    args = parser.parse_args()
    filelib.assert_exists_nz(args.mapability_file)
    filelib.assert_exists_nz(args.treatment_bam)
    filelib.assert_exists_nz(args.control_bam)

    if args.fragment_length:
        assert args.fragment_length > 0 and args.fragment_length < 10000

    # Set up directories to run it on.
    p("Setting up directories.\n")
    if not os.path.exists(args.outpath):
        os.mkdir(args.outpath)

    # Copy the mapability file to the outpath.
    shutil.copy2(args.mapability_file, args.outpath)
    
    # Do preprocessing for PeakSeq.
    p("Preprocessing.\n")
    treatment_preproc_path = os.path.join(args.outpath, "preprocess.treatment")
    control_preproc_path = os.path.join(args.outpath, "preprocess.control")
    if not os.path.exists(treatment_preproc_path):
        os.mkdir(treatment_preproc_path)
    if not os.path.exists(control_preproc_path):
        os.mkdir(control_preproc_path)
    x1 = make_peakseq_preproc_command(
        args.treatment_bam, treatment_preproc_path)
    x2 = make_peakseq_preproc_command(args.control_bam, control_preproc_path)
    x = shell.parallel([x1, x2])
    print x
    # Make sure expected files exist.
    x1 = os.path.join(treatment_preproc_path, "chr_ids.txt")
    x2 = os.path.join(control_preproc_path, "chr_ids.txt")
    filelib.assert_exists_nz(x1)
    filelib.assert_exists_nz(x2)


    # Make configuration file.
    p("Making configuration file.\n")
    config_file = os.path.join(args.outpath, "config.dat")
    make_config_file(
        config_file, treatment_preproc_path, control_preproc_path,
        args.mapability_file, experiment_name=args.experiment_name,
        fragment_length=args.fragment_length)

    # Run PeakSeq.
    p("Running PeakSeq in %s.\n" % args.outpath)
    cmd = make_peakseq_run_command(config_file)
    x = shell.single(cmd, path=args.outpath)
    print x

    p("Done.\n")
    

if __name__ == '__main__':
    main()
