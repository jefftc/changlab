#!/usr/bin/env python

def find_sppscript():
    import os
    from genomicode import config
    
    file_ = "sppscript.R"
    filename = os.path.join(config.changlab_Rlib, file_)
    assert os.path.exists(filename), "I could not find %s." % file_
    filename = os.path.realpath(filename)
    return filename


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
    import argparse
    from genomicode import filelib
    from genomicode import shell

    p = filelib.tswrite
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("treatment_bam", help="BAM file of treated sample.")
    parser.add_argument("control_bam", help="BAM file of background sample.")
    parser.add_argument("outpath", help="Directory to store the results.")

    parser.add_argument(
        "-j", dest="num_procs", type=int, default=1,
        help="Number of jobs to run in parallel.")
    parser.add_argument(
        "--fdr_cutoff", default=0.05, type=float, help="")

    args = parser.parse_args()
    filelib.assert_exists_nz(args.treatment_bam)
    filelib.assert_exists_nz(args.control_bam)
    args.treatment_bam = os.path.realpath(args.treatment_bam)
    args.control_bam = os.path.realpath(args.control_bam)

    assert args.num_procs >= 1 and args.num_procs < 100, \
           "Please specify between 1 and 100 processes."
    assert args.fdr_cutoff > 0.0 and args.fdr_cutoff < 1.0

    # Set up directories to run it on.
    p("Setting up directories.\n")
    if not os.path.exists(args.outpath):
        os.mkdir(args.outpath)

    # Run SPP.
    p("Running spp in %s.\n" % args.outpath)
    sq = shell.quote
    sppscript = find_sppscript()
    x = sq(args.treatment_bam), sq(args.control_bam), args.fdr_cutoff, \
        args.num_procs
    x = " ".join(map(str, x))
    cmd = "cat %s | R --vanilla %s" % (sppscript, x)
    x = shell.single(cmd, path=args.outpath)
    print x

    p("Done.\n")
    

if __name__ == '__main__':
    main()
