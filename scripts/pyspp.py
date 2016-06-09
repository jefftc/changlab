#!/usr/bin/env python

def find_sppscript():
    import os
    from genomicode import config
    
    file_ = "sppscript.R"
    filename = os.path.join(config.changlab_Rlib, file_)
    assert os.path.exists(filename), "I could not find %s." % file_
    filename = os.path.realpath(filename)
    return filename


def main():
    import os
    import argparse
    from genomicode import filelib
    from genomicode import parallel

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
    sq = parallel.quote
    sppscript = find_sppscript()
    x = sq(args.treatment_bam), sq(args.control_bam), args.fdr_cutoff, \
        args.num_procs
    x = " ".join(map(str, x))
    cmd = "cat %s | R --vanilla %s" % (sppscript, x)
    x = parallel.sshell(cmd, path=args.outpath)
    print x

    p("Done.\n")
    

if __name__ == '__main__':
    main()
