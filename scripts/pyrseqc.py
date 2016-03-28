#!/usr/bin/env python

ARG_BAM_FILE = "PYRSEQC_ARG___BAM_FILE"
ARG_FASTA_FILE = "PYRSEQC_ARG___FASTA_FILE"
ARG_OUT_PREFIX = "PYRSEQC_ARG___OUTPUT_PREFIX"
ARG_READ_LENGTH = "PYRSEQC_ARG___READ_LENGTH"
ARG_REF_GENE_MODEL = "PYRSEQC_ARG___REF_GENE_MODEL"
ARG_HOUSEKEEPING_GENE_MODEL = "PYRSEQC_ARG___HOUSEKEEPING_GENE_MODEL"


class Module:
    def __init__(self, script_name, params, outfile=None):
        self.script_name = script_name
        self.params = params
        self.outfile = outfile


def pairedend2PESE(is_paired):
    layout = "SE"
    if is_paired:
        layout = "PE"
    return layout


def script2cmd(
    rseqc_path, script_name, outpath, out_prefix, outfile, args, params):
    import os
    from genomicode import parallel

    assert os.path.exists(outpath)

    sq = parallel.quote

    # Convert script_name to full path.
    assert os.path.exists(rseqc_path)
    filename = os.path.join(rseqc_path, script_name)
    assert os.path.exists(filename)
    filename = sq(filename)

    # Clean up each of the params.
    clean = []
    for x in params:
        if x == ARG_BAM_FILE:
            x = args.bam_file
        elif x == ARG_FASTA_FILE:
            x = args.fasta_file
        elif x == ARG_OUT_PREFIX:
            x = os.path.join(outpath, out_prefix)
        elif x == ARG_READ_LENGTH:
            x = args.read_length
        elif x == ARG_REF_GENE_MODEL:
            x = args.ref_gene_model
        elif x == ARG_HOUSEKEEPING_GENE_MODEL:
            x = args.housekeeping_gene_model
        if type(x) is type(""):
            assert not x.startswith("PYRSEQC_ARG"), "Unhandled: %s" % x
        x = sq(str(x))
        clean.append(x)

    # RSeQC scripts use #!/usr/bin/python, which may not be the right
    # one.  Use the python on the path.
    cmd = ["python", filename] + clean
    cmd = " ".join(cmd)
    if outfile is not None:
        outfilename = os.path.join(outpath, outfile)
        cmd = "%s >& %s" % (cmd, outfilename)
        #cmd = "%s > %s" % (cmd, outfilename)
    return cmd

    

def main():
    import os
    import argparse
    import shutil
    from genomicode import config
    from genomicode import filelib
    from genomicode import parallel

    p = filelib.tswrite
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("bam_file", help="Should be sorted and indexed.")
    parser.add_argument("fasta_file", help="FASTA version of the bam_file.")
    parser.add_argument("read_length", help="Length of reads.")
    parser.add_argument("ref_gene_model", help="Gene model, in BED format.")
    parser.add_argument(
        "housekeeping_gene_model",
        help="Gene model for housekeeping genes, in BED format.")
    parser.add_argument("outpath", help="Directory to store the results.")

    parser.add_argument(
        "-j", dest="num_procs", type=int, default=1,
        help="Number of jobs to run in parallel.")
    parser.add_argument(
        "--paired_end", action="store_true",
        help="BAM file contains paired end reads.")
    parser.add_argument(
        "--dry_run", action="store_true",
        help="Show the commands to run and exit.")


    args = parser.parse_args()

    filelib.assert_exists_nz(args.bam_file)
    filelib.assert_exists_nz(args.fasta_file)
    assert args.read_length is not None, "--read_length must be provided."
    args.read_length = int(args.read_length)
    assert args.read_length >= 10 and args.read_length <= 1000, \
           "Invalid --read_length: %d" % args.read_length
    filelib.assert_exists_nz(args.ref_gene_model)
    filelib.assert_exists_nz(args.housekeeping_gene_model)

    assert args.num_procs >= 1 and args.num_procs < 100, \
           "Please specify between 1 and 100 processes."


    # Set up directories to run it on.
    p("Setting up directories.\n")
    if not os.path.exists(args.outpath):
        os.mkdir(args.outpath)

    # Run the libraries in RSeQC.
    rseqc_path = filelib.which_assert(config.rseqc_path)

    TIN_MIN_COVERAGE = 10

    # Set up each of the modules.
    M = Module
    MODULES = [
        M("clipping_profile.py", (
            "-i", ARG_BAM_FILE, "-o", ARG_OUT_PREFIX,
            "-s", pairedend2PESE(args.paired_end))),
        M("deletion_profile.py", (
            "-i", ARG_BAM_FILE, "-o", ARG_OUT_PREFIX,
            "-l", ARG_READ_LENGTH)),
        M("insertion_profile.py", (
            "-i", ARG_BAM_FILE, "-o", ARG_OUT_PREFIX,
            "-s", pairedend2PESE(args.paired_end))),
        # Error: object 'A2C' not found
        # This means BAM file doesn't have MD tags.
        M("mismatch_profile.py", (
            "-i", ARG_BAM_FILE, "-o", ARG_OUT_PREFIX,
            "-l", ARG_READ_LENGTH)),
        #M("geneBody_coverage.py", (
        #    "-i", ARG_BAM_FILE, "-o", ARG_OUT_PREFIX,
        #    "-r", ARG_HOUSEKEEPING_GENE_MODEL)),
        M("geneBody_coverage2.py", (
            "-i", ARG_BAM_FILE, "-o", ARG_OUT_PREFIX,
            "-r", ARG_HOUSEKEEPING_GENE_MODEL)),
        M("infer_experiment.py", (
            "-i", ARG_BAM_FILE, "-r", ARG_REF_GENE_MODEL),
          "infer_experiment.txt"),
        M("inner_distance.py", (
            "-i", ARG_BAM_FILE, "-o", ARG_OUT_PREFIX,
            "-r", ARG_REF_GENE_MODEL)),
        M("read_distribution.py", (
            "-i", ARG_BAM_FILE, "-r", ARG_REF_GENE_MODEL),
          "read_distribution.txt"),
        M("read_duplication.py", (
            "-i", ARG_BAM_FILE, "-o", ARG_OUT_PREFIX)),
        M("read_GC.py", (
            "-i", ARG_BAM_FILE, "-o", ARG_OUT_PREFIX)),
        M("read_hexamer.py", ("-i", ARG_FASTA_FILE), "hexamer.txt"),
        M("read_NVC.py", (
            "-i", ARG_BAM_FILE, "-o", ARG_OUT_PREFIX)),
        # Doesn't work.  Uses R process that takes 100's of Gb RAM.
        #M("read_quality.py", (
        #    "-i", ARG_BAM_FILE, "-o", ARG_OUT_PREFIX)),
        M("RNA_fragment_size.py", (
            "-i", ARG_BAM_FILE, "-r", ARG_REF_GENE_MODEL), "fragment.txt"),
        M("tin.py", (
            "-i", ARG_BAM_FILE, "-r", ARG_REF_GENE_MODEL,
            "-c", TIN_MIN_COVERAGE), "tin.txt"),
        
        # Not implemented
        # RPKM_saturation.  Needs pairedness info
        
        # Splicing stuff.
        # Not working.
        #M("junction_annotation.py", (
        #    "-i", ARG_BAM_FILE, "-o", ARG_OUT_PREFIX,
        #    "-r", ARG_REF_GENE_MODEL)),
        # Not working.
        #M("junction_saturation.py", (
        #    "-i", ARG_BAM_FILE, "-o", ARG_OUT_PREFIX,
        #    "-r", ARG_REF_GENE_MODEL)),
        ]

    commands = []
    for x in MODULES:
        out_prefix = "rseqc"
        cmd = script2cmd(
            rseqc_path, x.script_name, args.outpath, out_prefix, x.outfile,
            args, x.params)
        commands.append(cmd)

    if args.dry_run:
        for x in commands:
            print x
        return

    x1 = ""
    if len(commands) > 1:
        x1 = "s"
    x2 = "."
    if args.num_procs > 1:
        x2 = " on %d processors." % args.num_procs
    x3 = "Running %d command%s%s" % (len(commands), x1, x2)
    p("%s\n" % x3)
    x = parallel.pshell(commands, max_procs=args.num_procs)

    # tin.py will save some files in the current directory.  Look for
    # them and move them into the outpath.
    x, f = os.path.split(args.bam_file)
    f, x = os.path.splitext(f)
    TIN_FILES = [
        "%s.summary.txt" % f,
        "%s.tin.xls" % f,
        ]
    for f in TIN_FILES:
        if os.path.exists(f):
            shutil.move(f, args.outpath)

    p("Done.\n")
    

if __name__ == '__main__':
    main()
