#!/usr/bin/env python

import sys, os

# read_signatures
# process_gp_imod_all_vars
#
# make_file_layout
# init_paths
# 
# make_pybinreg_cmd
# run_one_pybinreg
# run_many_pybinreg
# check_pybinreg
#
# extract_reports
# summarize_probabilities
# summarize_heatmap
# summarize_signatures
# summarize_report
#
# _hash_name

# The structure of "signatures" is based on the signatures.txt file.
# process_gp_imod_all_vars can introduce a boolean variable "Changed"
# that indicates whether the default variables of the signature were
# changed.
def read_signatures(
    sigdb_path, desired_normalization, desired_ids, desired_tags):
    # Read the signatures and return the ones that match the
    # specifications.  Always obey desired_normalization.  If
    # desired_ids is specified, then do the ones with that ID and
    # ignore desired_tags.  Otherwise, obey desired_tags.
    
    from genomicode import filelib
    
    opr = os.path.realpath
    opj = os.path.join

    filename = opj(sigdb_path, "signatures.txt")
    assert os.path.exists(filename), "Missing signatures.txt file."

    x = [x.upper() for x in desired_normalization]
    desired_normalization = {}.fromkeys(x)

    x = [x.upper() for x in desired_ids]
    desired_ids = {}.fromkeys(x)

    x = [x.upper() for x in desired_tags]
    desired_tags = {}.fromkeys(x)

    ds = []
    for d in filelib.read_row(filename, header=1):
        # xls2txt converts all values to floats.  Convert them back to
        # integers.
        d.xID = int(float(d.xID))
        d.Genes = int(float(d.Genes))
        d.Metagenes = int(float(d.Metagenes))

        # Skip if not the right normalization.
        if d.Normalization.upper() not in desired_normalization:
            continue

        # If the IDs are specified, then make sure the id matches.  If
        # no IDs are specified, then make sure the tags match.
        if desired_ids:
            if str(d.xID) not in desired_ids:
                continue
        else:
            tags = [x.upper() for x in d.Tags.split()]
            for tag in tags:
                if tag in desired_tags:
                    break
            else:
                # None of these tags matched any of the desired ones.
                continue
        
        # Skip if not all parameters supplied.
        if not d.Normalization:
            continue
        if not d.Genes or not d.Metagenes:
            continue
        if not d.Metagenes or not d.Quantile:
            continue
        if not d.Train0 or not d.Train1:
            continue

        # Find the training files.  If not found, then skip.
        train0 = opr(opj(sigdb_path, d.Train0))
        train1 = opr(opj(sigdb_path, d.Train1))
        if not os.path.exists(train0):
            train0 = train0 + ".gz"
        if not os.path.exists(train1):
            train1 = train1 + ".gz"
        if not os.path.exists(train0) or not os.path.exists(train1):
            continue
        d.Train0 = train0
        d.Train1 = train1

        ds.append(d)
    return ds

def process_gp_imod_all_vars(gp_imod_all_vars, signatures, why_dropped):
    import copy
    import urllib

    why_dropped = copy.copy(why_dropped)
    
    # Make a list of the options.
    options = []  # list of (key, value)
    assert gp_imod_all_vars.find("&") >= 0
    keyvalue = gp_imod_all_vars.split("&")
    for x in keyvalue:
        k, v = x.split("=", 1)
        k, v = k.strip(), v.strip()
        # Unquote after we split on "&" and "=".
        k = urllib.unquote_plus(k)
        v = urllib.unquote_plus(v)
        options.append((k, v))

    # See if the user wanted to do all signatures (without modification).
    do_all_signatures = False
    for key, value in options:
        if key == "which_signatures" and value == "all of them":
            do_all_signatures = True
            break
    # If the user wanted to do all the signatures, then just return
    # all the original signatures.
    if do_all_signatures:
        return signatures, why_dropped
    
    # Filter out some parameters to ignore.
    options_clean = []
    for key, value in options:
        if key == "which_signatures":
            continue
        if key in [
            "rma_expression_file_cb", "rma_expression_file_url",
            "mas5_expression_file_cb", "mas5_expression_file_url"]:
            continue
        options_clean.append((key, value))
    options = options_clean

    # Options are in the format:
    # sig_<Name>                                   "yes (default parameters)"
    #                                              "yes (custom parameters)"
    #                                              "no"
    # sig_<Name>_num_genes
    # sig_<Name>_num_metagenes
    # sig_<Name>_apply_quantile_normalization      "yes" "no"
    # sig_<Name>_apply_shiftscale_normalization    "yes" "no"
    #
    # If "default parameters" or "no" are requested, then the other
    # parameters will not be given.

    # Reformat the parameters for easy processing.
    # Name -> "yesno", "num_genes", "num_metagenes", ... -> value
    # Any of these values can be missing.
    # Set some constants for convenience.
    YESNO = "yesno"
    NUM_GENES = "num_genes"
    NUM_METAGENES = "num_metagenes"
    QUANTILE_NORM = "apply_quantile_normalization"
    SHIFTSCALE_NORM = "apply_shiftscale_normalization"
    PARAMETER_NAMES = [
        YESNO, NUM_GENES, NUM_METAGENES, QUANTILE_NORM, SHIFTSCALE_NORM]
    name2options = {}
    for key, value in options:
        assert key.startswith("sig_"), "Unknown key: %s" % key
        x = key[4:]
        x = x.split("_", 1)   # e.g. <name>, "num_genes"
        name = x[0]   # signature name

        # Make sure the parameters are valid.
        parameter = YESNO
        if len(x) == 2:
            parameter = x[1]
        assert parameter in PARAMETER_NAMES
        
        if name not in name2options:
            name2options[name] = {}
        assert parameter not in name2options[name], \
               "%s has duplicate values for %s." % (name, parameter)
        name2options[name][parameter] = value

    # Process the parameters for each of the signatures.
    signatures_clean = []
    for sig in signatures:
        sig = copy.copy(sig)
        options = name2options.get(sig.Name, {})
        # By default, let users choose custom values for each signature.
        yesno = options.get(YESNO, "custom")
        if yesno.find("default") >= 0:
            # If the user wants to use the default parameters, then
            # ignore the values of the parameters provided.
            signatures_clean.append(sig)
            continue
        if yesno == "no":
            # If the user wanted to skip this signature, then ignore
            # all other parameters provided.
            #assert len(name2options[sig.Name]) == 1  # should be just "yesno"
            print "Signature %s disabled by user.  Skipping" % sig.Name
            x = "User chose to skip %s signature." % sig.Name
            why_dropped[sig.xID] = x
            continue

        # I do no checking on the values of the variables.  I will
        # let pybinreg take care of that.  If there is a problem,
        # pybinreg will raise an Exception and cause the
        # processing for this analysis to be aborted.
        changed = False
        if NUM_GENES in options:
            x = int(options[NUM_GENES])
            if sig.Genes != x:
                print "Changing %s '%s' from %d to %d." % (
                    sig.Name, NUM_GENES, sig.Genes, x)
                sig.Genes = x
                changed = True
        if NUM_METAGENES in options:
            x = int(options[NUM_METAGENES])
            if sig.Metagenes != x:
                print "Changing %s '%s' from %d to %d." % (
                    sig.Name, NUM_METAGENES, sig.Metagenes, x)
                sig.Metagenes = x
                changed = True
        if QUANTILE_NORM in options:
            assert options[QUANTILE_NORM].lower() in ["yes", "no"]
            # sig.Quantile should be "Yes" or "No"
            x = options[QUANTILE_NORM].capitalize()
            if sig.Quantile != x:
                print "Changing %s '%s' from '%s' to '%s'." % (
                    sig.Name, QUANTILE_NORM, sig.Quantile, x)
                sig.Quantile = x
                changed = True
        if SHIFTSCALE_NORM in options:
            assert options[SHIFTSCALE_NORM].lower() in ["yes", "no"]
            # sig.Shift_Scale should be "Yes" or "No"
            x = options[SHIFTSCALE_NORM].capitalize()
            if sig.Shift_Scale != x:
                print "Changing %s '%s' from '%s' to '%s'." % (
                    sig.Name, SHIFTSCALE_NORM, sig.Shift_Scale, x)
                sig.Shift_Scale = x
                changed = True
        if changed:
            # Flag this as being changed.
            setattr(sig, "Changed", True)
            #print "Signature %s run with altered parameters." % sig.Name
            
        signatures_clean.append(sig)
    signatures = signatures_clean

    return signatures, why_dropped

def make_file_layout(outpath):
    from genomicode import filelayout

    outpath = outpath or "."
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    outpath = os.path.realpath(outpath)

    Path, File = filelayout.Path, filelayout.File
    
    # Have only set set of these files for the whole analysis.
    FILES = [
        Path.ATTIC("attic"),
        File.PROBABILITIES_PCL("probabilities.pcl"),
        File.PROBABILITIES_PNG("probabilities.png"),
        File.PROBABILITIES_CDT("probabilities.cdt"),  # made by clustering
        File.PROBABILITIES_GTR("probabilities.gtr"),  # made by clustering
        File.PROBABILITIES_ATR("probabilities.atr"),  # made by clustering
        File.DATASET_RMA("dataset.rma.gct"),
        File.DATASET_MAS5("dataset.mas5.gct"),
        File.PARAMETERS("parameters.txt"),
        File.REPORT("REPORT.html"),
        File.FILES("FILES"),
        ]
    file_layout = Path.OUTPATH(outpath, *FILES)
    return file_layout

def init_paths(file_layout):
    from genomicode import filelayout

    for x in filelayout.walk(file_layout):
        dirpath, dirnames, filenames = x
        if os.path.exists(dirpath):
            continue
        os.mkdir(dirpath)

def make_pybinreg_cmd(
    pybinreg, python, binreg_path, matlab, arrayplot, povray, cluster, libpath,
    outpath, archive, num_genes, num_metagenes,
    quantile_normalize, shift_scale_normalize,
    train0_file, train1_file, test_file):
    pybinreg = pybinreg or "pybinreg.py"

    cmd = []
    # If the python path was explicitly specified, then make sure to
    # use this one to launch pybinreg.  Somehow GenePattern is not
    # getting the path correctly anymore (100228).
    if python:
        cmd.append(python)
    cmd.append(pybinreg)
    if python:
        cmd.append("--python=%s" % python)
    if binreg_path:
        cmd.append("--binreg=%s" % binreg_path)
    if matlab:
        cmd.append("--matlab=%s" % matlab)
    if arrayplot:
        cmd.append("--arrayplot=%s" % arrayplot)
    if povray:
        cmd.append("--povray=%s" % povray)
    if cluster:
        cmd.append("--cluster=%s" % cluster)
    for path in libpath:
        cmd.append("--libpath=%s" % path)
    cmd.append("-o %s" % outpath)
    if archive:
        cmd.append("-z")
    cmd.append("-g %d" % num_genes)
    cmd.append("-m %d" % num_metagenes)
    
    if quantile_normalize:
        cmd.append("-q")
    if shift_scale_normalize:
        cmd.append("-s")
    cmd.append("'%s'" % train0_file)
    cmd.append("'%s'" % train1_file)
    cmd.append("'%s'" % test_file)
    cmd = " ".join(cmd)

    return cmd

def run_one_pybinreg(cmd, outpath, outfile):
    import subprocess

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    print cmd
    p = subprocess.Popen(
        cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    w, r = p.stdin, p.stdout
    w.close()

    binreg_out = r.read()
    open(outfile, 'w').write(binreg_out)
    check_pybinreg(outpath, outfile)

def run_many_pybinreg(jobs, num_procs):
    # jobs should be a list of cmd, outpath, outfile.
    from genomicode import parallel

    commands = []
    for x in jobs:
        cmd, outpath, outfile = x
        
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        cmd = "%s >& %s" % (cmd, outfile)
        commands.append(cmd)

    parallel.run_commands(commands, num_procs)

    for x in jobs:
        cmd, outpath, outfile = x
        check_pybinreg(outpath, outfile)

def check_pybinreg(outpath, outfile):
    if os.path.exists(os.path.join(outpath, "probabilities.txt")):
        # BinReg completed successfully.
        return
    # Problem with BinReg.
    print >>sys.stderr, "I could not find a probabilities.txt file in %s." % (
        outpath)

    if not os.path.exists(outfile):
        print >>sys.stderr, "I could not find any output in %s." % outpath
        return
    
    # Try to find a Traceback.
    found_tb = False
    traceback = []
    for line in open(outfile):
        if line.startswith("Traceback (most recent call last):"):
            found_tb = True
        if found_tb:
            traceback.append(line)
    if traceback:
        print >>sys.stderr, "I found a Traceback in the pybinreg output."
        for line in traceback:
            print >>sys.stderr, "pybinreg>  %s" % line.rstrip()

    error = "It looks like pybinreg failed."
    if traceback:
        x = traceback[-1]
        if x.find(":") >= 0:
            x = x[x.index(":")+1:].strip()
        error = "pybinreg error: %s" % x
            
    raise AssertionError, "%s\nProcessing aborted." % error

def extract_reports(names, paths, file_layout):
    import shutil
    import re
    
    # From each signature, pull out the following files:
    # REPORT.html
    # signature.png
    # predictions.png
    # probabilities.txt
    #
    # Copy to:
    # SIG19_AKT_REPORT.html
    # SIG19_AKT_signature.png
    # ... etc
    #
    # Fix the links in SIG19_AKT_REPORT.html.
    # Can be either HREFs or IMG.
    files = [
        "REPORT.html",
        "signature.png",
        "predictions.png",
        "probabilities.txt",
        ]
    outpath = file_layout.OUTPATH
    report_files = []
    for name, path in zip(names, paths):
        infiles = files
        outfiles = ["%s_%s" % (name, x) for x in files]
        infiles_full = [os.path.join(path, x) for x in infiles]
        outfiles_full = [os.path.join(outpath, x) for x in outfiles]
        for x in infiles_full:
            assert os.path.exists(x)
        for src, dst in zip(infiles_full, outfiles_full):
            shutil.copy2(src, dst)

        # Fix the title and links in REPORT.html.
        x = [x for x in outfiles_full if x.endswith("_REPORT.html")]
        assert len(x) == 1, "I could not find the report file."
        report_file = x[0]
        report_files.append(report_file)
        data = open(report_file).read()

        # Fix the title.
        x = re.sub(r"^SIG\d+_", "", name)
        src = r"<H1><EM>CreateSignatures</EM> Report</H1>"
        dst = r"<H1><EM>CreateSignatures</EM> Report for %s</H1>" % x
        assert data.find(src) >= 1, "I could not find the title of the report"
        data = data.replace(src, dst)

        # Fix the links.
        for infile, outfile in zip(infiles, outfiles):
            infile = infile.replace(".", "\\.")
            data = re.sub(
                infile, outfile, data, re.IGNORECASE|re.MULTILINE)
        open(report_file, 'w').write(data)
        
        
    assert len(report_files) == len(names)
    return report_files

def summarize_probabilities(signatures, names, paths, file_layout):
    from genomicode import filelib
    from genomicode import binreg

    _hash = binreg._hash_sampleid
    
    sample_names = []   # list of sample names
    probabilities = []  # matrix of probabilities
    for i, sig in enumerate(signatures):
        name, outpath = names[i], paths[i]
        filename = os.path.join(outpath, "probabilities.txt")
        assert os.path.exists(filename), \
               "Could not find probability file for %s." % name
        ds = [d for d in filelib.read_row(filename, header=1)]
        ds = [d for d in ds if d.Type == "test"]

        # Assign and check the sample names.
        if not sample_names:
            sample_names = [d.Sample for d in ds]
        assert len(sample_names) == len(ds)
        for i, d in enumerate(ds):
            assert _hash(d.Sample) == _hash(sample_names[i])

        # Bug: What is there are nan's here?
        probs = [float(d.Probability) for d in ds]
        probabilities.append(probs)

    # Write out the probability file.
    handle = open(file_layout.PROBABILITIES_PCL, 'w')
    x = ["SigID", "NAME"] + sample_names
    print >>handle, "\t".join(x)
    for i, sig in enumerate(signatures):
        name = sig.Name
        if getattr(sig, "Changed", False):
            name = "%s*" % name
        # If the signature was modified, make a notation.
        x = [sig.xID, name] + probabilities[i]
        print >>handle, "\t".join(map(str, x))
    handle.close()

def summarize_heatmap(python, arrayplot, cluster, libpath, file_layout):
    from genomicode import graphlib

    # Bug: what if there are nan's in the probabilities?
    xpix, ypix = 30, 30
    x = graphlib.plot_heatmap(
        file_layout.PROBABILITIES_PCL, file_layout.PROBABILITIES_PNG,
        xpix, ypix, color="bild", show_colorbar=True, show_grid=True,
        gene_label=True, cluster_genes=True,
        array_label=True, cluster_arrays=True, scale=-0.5, gain=2.0,
        no_autoscale=True,
        python=python, arrayplot=arrayplot, cluster=cluster,
        libpath=libpath)
    print x
    
    # Clean up some extra files.
    if os.path.exists(file_layout.PROBABILITIES_CDT):
        src = file_layout.PROBABILITIES_CDT
        x = os.path.split(file_layout.PROBABILITIES_CDT)[1]
        dst = os.path.join(file_layout.ATTIC, x)
        os.rename(src, dst)
    if os.path.exists(file_layout.PROBABILITIES_GTR):
        src = file_layout.PROBABILITIES_GTR
        x = os.path.split(file_layout.PROBABILITIES_GTR)[1]
        dst = os.path.join(file_layout.ATTIC, x)
        os.rename(src, dst)
    if os.path.exists(file_layout.PROBABILITIES_ATR):
        src = file_layout.PROBABILITIES_ATR
        x = os.path.split(file_layout.PROBABILITIES_ATR)[1]
        dst = os.path.join(file_layout.ATTIC, x)
        os.rename(src, dst)
    

def summarize_signatures(signatures, file_layout):
    handle = open(file_layout.PARAMETERS, 'w')
    x = ("xID", "Name", "Customized", "Genes", "Metagenes", "Preprocessing",
         "Quantile Norm", "Shift-Scale Norm", "Train0", "Train1")
    print >>handle, "\t".join(x)
    for sig in signatures:
        modified = "No"
        if getattr(sig, "Changed", False):
            modified = "Yes"
        train0 = os.path.split(sig.Train0)[1]
        train1 = os.path.split(sig.Train1)[1]
        x = (sig.xID, sig.Name, modified, sig.Genes, sig.Metagenes,
             sig.Normalization, sig.Quantile, sig.Shift_Scale, train0, train1)
        print >>handle, "\t".join(map(str, x))
    handle.close()

def summarize_report(
    signatures, orig_signatures, report_files, start_time, why_dropped,
    file_layout):
    import time
    import subprocess
    from genomicode import parselib
    from genomicode import htmllib

    def highlight(s):
        return htmllib.SPAN(s, style="background-color:yellow")
    def smaller(s):
        return htmllib.FONT(s, size=-1)

    id2orig = {}
    for sig in orig_signatures:
        id2orig[sig.xID] = sig

    id2new = {}
    for sig in signatures:
        id2new[sig.xID] = sig

    assert len(signatures) == len(report_files)
    id2reportfile = {}
    for sig, file in zip(signatures, report_files):
        # The report_file in the HTML should be a relative path.
        x, file = os.path.split(file)
        id2reportfile[sig.xID] = file

    # Figure out which of the signatures were dropped.
    missing_ids = []
    for sig in orig_signatures:
        if sig.xID in id2new:
            continue
        missing_ids.append(sig.xID)

    # Make a list of all the signatures.
    all_ids = {}.fromkeys(id2orig.keys() + id2new.keys())
    schwartz = [(id2orig[x].Name, x) for x in all_ids]
    schwartz.sort()
    all_ids = [x[-1] for x in schwartz]
    

    lines = []
    w = lines.append
    w("<HTML>")
    w(htmllib.HEAD(htmllib.TITLE("ScoreSignatures Report")))
    w("<BODY>")
    w(htmllib.CENTER(htmllib.H1(htmllib.EM("ScoreSignatures") + " Report")))

    w(htmllib.H3("I.  Signatures"))

    # Make a table with each of the signatures.
    rows = []

    x = htmllib.TR(
        htmllib.TH("ID", align="LEFT") +
        htmllib.TH("Signature", align="LEFT") +
        htmllib.TH("Preprocessing", align="LEFT") +
        htmllib.TH("Genes", align="LEFT") +
        htmllib.TH("Metagenes", align="LEFT") +
        htmllib.TH("Normalization", align="LEFT")
        )
    rows.append(x)
    
    which_changed = {}  # ID -> 1
    for id in all_ids:
        orig = id2orig[id]
        sig = id2new.get(id)

        cols = []

        # ID
        cols.append(htmllib.TD(orig.xID))

        # Name
        name = orig.Name
        report_file = None
        if sig:
            report_file = id2reportfile.get(sig.xID)
        if report_file:
            name = htmllib.A(name, href=report_file)
        cols.append(htmllib.TD(name))

        # If this signature was not run, then skip the rest of the columns.
        if not sig:
            x = why_dropped.get(orig.xID, "Skipped for unknown reason.")
            x = htmllib.TD(highlight(x), colspan=4)
            cols.append(x)
            rows.append(
                htmllib.TR("\n".join(cols)))
            continue
        
        # Preprocessing
        x = sig.Normalization
        if sig.Normalization != orig.Normalization:
            which_changed[sig.xID] = 1
            x = "%s<BR>%s" % (
                highlight(sig.Normalization),
                smaller(htmllib.EM("default: %s" % orig.Normalization)))
        cols.append(htmllib.TD(x))

        # Genes
        x = sig.Genes
        if sig.Genes != orig.Genes:
            which_changed[sig.xID] = 1
            x = "%s<BR>%s" % (
                highlight(sig.Genes),
                smaller(htmllib.EM("default: %s" % orig.Genes)))
        cols.append(htmllib.TD(x))
        
        # Metagenes
        x = sig.Metagenes
        if sig.Metagenes != orig.Metagenes:
            which_changed[sig.xID] = 1
            x = "%s<BR>%s" % (
                highlight(sig.Metagenes),
                smaller(htmllib.EM("default: %s" % orig.Metagenes)))
        cols.append(htmllib.TD(x))

        # Normalization
        norm = []
        if sig.Quantile.upper() == "YES":
            norm.append("Quantile")
        if sig.Shift_Scale.upper() == "YES":
            norm.append("Shift-Scale")
        norm_str = "None"
        if norm:
            norm_str = " and ".join(norm)
        if sig.Quantile.upper() != orig.Quantile.upper() or \
           sig.Shift_Scale.upper() != orig.Shift_Scale.upper():
            which_changed[sig.xID] = 1
            norm = []
            if orig.Quantile.upper() == "YES":
                norm.append("Quantile")
            if orig.Shift_Scale.upper() == "YES":
                norm.append("Shift-Scale")
            x = "None"
            if norm:
                x = " and ".join(norm)
            norm_str = "%s<BR>%s" % (
                highlight(norm_str), smaller(htmllib.EM("default: %s" % x)))
        cols.append(htmllib.TD(norm_str))

        #assert sig_changed == getattr(sig, "Changed", False), "%s %s %s" % (
        #   sig.Name, sig_changed, getattr(sig, "Changed", "missing"))
        x = htmllib.TR("\n".join(cols))
        rows.append(x)
        
    w(htmllib.TABLE("\n".join(rows), border=1, cellpadding=3, cellspacing=0))

    w(htmllib.P())
    w(htmllib.B("Table 1: Signatures Analyzed."))
    if not which_changed:
        w("All signatures were run with the default parameters, "
          "as shown above.")
    else:
        w("The customized parameters are highlighted in yellow.")

    w(htmllib.P())

    w(htmllib.H3("II.  Results"))

    prob_file = os.path.split(file_layout.PROBABILITIES_PNG)[1]
    w(htmllib.A(htmllib.IMG(height=768, src=prob_file), href=prob_file))

    w(htmllib.P())

    w(htmllib.B("Figure 1: Predictions."))
    w("In this heatmap, each row contains a signature and each column "
      "contains a sample from your data set.")
    if which_changed:
        #names = sorted([id2orig[x].Name for x in which_changed])
        w("The asterisks denote the signatures that were run with "
          "customized parameters.")
    w("The color corresponds to the probability that a pathway is activated "
      "in a sample.")
    w("Warm colors represent high probabilities, and cool colors low.\n")

    w(htmllib.P())
    prob_file = os.path.split(file_layout.PROBABILITIES_PCL)[1]
    w("The raw values from this plot are available as a "
      'PCL-formatted file: %s' % htmllib.A(prob_file, href=prob_file))

    # Format the current time.
    end_time = time.time()
    time_str = parselib.pretty_date(start_time)
    x = int(end_time-start_time)
    num_min = x / 60
    num_secs = x % 60
    if num_min == 0:
        run_time = "%ss" % parselib.pretty_int(num_secs)
    else:
        run_time = "%sm %ss" % (parselib.pretty_int(num_min), num_secs)

    # Get the hostname.
    cmd = "hostname"
    p = subprocess.Popen(
        cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    wh, r = p.stdin, p.stdout
    wh.close()
    hostname = r.read().strip()
    assert hostname, "I could not get the hostname."

    w(htmllib.P())
    w(htmllib.HR())
    w(htmllib.EM(
        "This analysis was run on %s on %s.  It took %s to complete.\n" %
        (time_str, hostname, run_time)))
    
    w("</BODY>")
    w("</HTML>")

    x = "\n".join(lines) + "\n"
    outfile = file_layout.REPORT
    open(outfile, 'w').write(x)

def _hash_name(name):
    import re
    # Fix the header to be a python variable.
    x = name
    # Replace all non-word character with _.
    x = re.sub(r"\W", "_", x)
    # Replace initial numbers with Xnumber.
    x = re.sub(r"^(\d)", r"X\1", x)
    return x

def main():
    from optparse import OptionParser, OptionGroup
    
    usage = "usage: %prog [options] sigdb_path"
    parser = OptionParser(usage=usage, version="%prog 01")

    parser.add_option(
        "-r", "--rma", dest="rma_dataset", type="string", default=None,
        help="Specify the RMA-normalized data to analyze.")
    parser.add_option(
        "-m", "--mas5", dest="mas5_dataset", type="string", default=None,
        help="Specify the MAS5-normalized data to analyze.")
    parser.add_option(
        "", "--sigid", dest="signature_ids", default=[], action="append",
        help="Specify a specific signature to use.")
    parser.add_option(
        "", "--max_signatures", dest="max_signatures", type="int",
        default=None,
        help="Maximum number of signatures to run (for DEBUGGING).")
    parser.add_option(
        "-j", "", dest="num_procs", type="int", default=1,
        help="Number of jobs to run in parallel.")
    parser.add_option(
        "-z", "", dest="archive", action="store_true", default=False,
        help="Archive the individual signatures.  Helpful for GenePattern.")
    parser.add_option(
        "", "--libpath", dest="libpath", action="append", default=[],
        help="Add to the Python library search path.")
    parser.add_option(
        "-o", "--outpath", dest="outpath", type="string", default=None,
        help="Save files in this path.")
    parser.add_option(
        "", "--gp_imod_all_vars", dest="gp_imod_all_vars", type="string",
        default=None,
        help="Special internal variable for use with GenePattern "
        "interactive modules.")
    parser.add_option(
        "", "--debug_gp_imod_all_vars", action="store_true", default=False, 
        dest="debug_gp_imod_all_vars",
        )
    
    #group = OptionGroup(parser, "Normalization")
    #group.add_option(
    #    "", "--normalization", dest="normalization", default="MAS5",
    #    help="How was the data set normalized (default MAS5).")
    #group.add_option(
    #    "-l", "--log_data", dest="log_data", action="store_true",
    #    default=False,
    #    help="Log the MAS5 data before analyzing.")
    #parser.add_option_group(group)

    group = OptionGroup(parser, "Pybinreg")
    group.add_option(
        "", "--pybinreg", dest="pybinreg", default=None,
        help="Specify the command to run pybinreg.py.")
    group.add_option(
        "", "--python", dest="python", default=None,
        help="Specify the command to run python.")
    group.add_option(
        "", "--binreg", dest="binreg_path", default=None,
        help="Specify the path to the BinReg2.0 code.")
    group.add_option(
        "", "--matlab", dest="matlab", default=None,
        help="Specify the command to run matlab.")
    group.add_option(
        "", "--arrayplot", dest="arrayplot", default=None,
        help="Specify the command to run arrayplot.")
    group.add_option(
        "", "--povray", dest="povray", default=None,
        help="Specify the command to run povray.")
    group.add_option(
        "", "--cluster", dest="cluster", default=None,
        help="Specify the command to run cluster.")
    parser.add_option_group(group)

    options, args = parser.parse_args()
    if len(args) < 1:
        #print sys.argv
        #print len(args), args
        parser.error("Please specify sigdb_path.")
    elif len(args) > 1:
        parser.error("Too many arguments.")
    sigdb_path, = args
    assert os.path.exists(sigdb_path), \
           "I could not find the signatures database: %s." % sigdb_path
    sigdb_path = os.path.realpath(sigdb_path)

    # DEBUG the gp_imod_all_vars variable.
    if options.debug_gp_imod_all_vars:
        assert not options.gp_imod_all_vars
        options.gp_imod_all_vars = (
            "mas5_expression_file_cb=file&mas5_expression_file_url=&"
            "rma_expression_file_cb=file&rma_expression_file_url=&"
            # Skip AKT signature.
            "sig_AKT=no&"
            # Change BCAT normalization.
            "sig_BCAT=yes (custom parameters)&"
            "sig_BCAT_apply_quantile_normalization=no&"
            "sig_BCAT_apply_shiftscale_normalization=no&"
            "sig_BCAT_num_genes=85&sig_BCAT_num_metagenes=2&"
            # No changes in E2F1.
            "sig_E2F1=yes (custom parameters)&"
            "sig_E2F1_apply_quantile_normalization=yes&"
            "sig_E2F1_apply_shiftscale_normalization=yes&"
            "sig_E2F1_num_genes=150&sig_E2F1_num_metagenes=2&"
            # Change genes in EGFR.
            "sig_EGFR=yes (custom parameters)&"
            "sig_EGFR_apply_quantile_normalization=no&"
            "sig_EGFR_apply_shiftscale_normalization=yes&"
            #"sig_EGFR_num_genes=50000&sig_EGFR_num_metagenes=2&"
            "sig_EGFR_num_genes=501&sig_EGFR_num_metagenes=2&"
            # Change quantile, genes, metagenes in ER.
            "sig_ER=yes (custom parameters)&"
            "sig_ER_apply_quantile_normalization=no&"
            "sig_ER_apply_shiftscale_normalization=yes&"
            "sig_ER_num_genes=150&sig_ER_num_metagenes=3&"
            "sig_HER2=yes (default parameters)&"
            "sig_IFNalpha=yes (default parameters)&"
            "sig_IFNgamma=yes (default parameters)&"
            "sig_MYC=yes (default parameters)&"
            "sig_P53=yes (default parameters)&"
            "sig_P63=yes (default parameters)&"
            "sig_PI3K=yes (default parameters)&"
            "sig_PR=yes (default parameters)&"
            "sig_RAS=yes (default parameters)&"
            "sig_SRC=yes (default parameters)&"
            "sig_STAT3=yes (default parameters)&"
            "sig_TGFB=yes (default parameters)&"
            "sig_TNFa=yes (default parameters)&"
            "which_signatures=I choose myself"
            )
        
    datafile_rma = datafile_mas5 = None
    if options.rma_dataset is not None:
        assert os.path.exists(options.rma_dataset), "RMA file not found."
        datafile_rma = os.path.realpath(options.rma_dataset)
    if options.mas5_dataset is not None:
        assert os.path.exists(options.mas5_dataset), "MAS5 file not found."
        datafile_mas5 = os.path.realpath(options.mas5_dataset)
    assert datafile_rma or datafile_mas5, \
           "Please specify an RMA and/or a MAS5 normalized data set."

    if options.libpath:
        sys.path = options.libpath + sys.path
    # Import after the library path is set.
    import time
    import arrayio
    from genomicode import jmath
    from genomicode import parallel
    from genomicode import filelib
    from genomicode import archive
    from genomicode import binreg
    from genomicode import genepattern

    start_time = time.time()
    
    genepattern.fix_environ_path()
    
    file_layout = make_file_layout(options.outpath)
    init_paths(file_layout)

    # Read the signatures and select the ones to score.
    # BUG: Should allow this to be specified on the command line.
    desired_tags = ["Production"]
    all_normalization = ["RMA", "MAS5"]
    desired_normalization = []
    if datafile_rma is not None:   # RMA datafile is specified.
        desired_normalization.append("RMA")
    if datafile_mas5 is not None:  # MAS5 datafile is specified.
        desired_normalization.append("MAS5")
        
    # If any signature IDs are specified, then use only those IDs and
    # ignore the desired tags.
    print "Reading signature database."
    desired_ids = []
    if options.signature_ids:
        desired_ids = options.signature_ids[:]
    x = read_signatures(
        sigdb_path, all_normalization, desired_ids, desired_tags)
    signatures = x
    orig_signatures = signatures[:]
    assert signatures, "No signatures available."

    # Filter for just the normalization that we have data files for.
    # Keep track of why we filtered out certain signatures.
    why_dropped = {}  # ID -> explanation as string
    good = []
    for sig in signatures:
        if sig.Normalization.upper() in desired_normalization:
            good.append(sig)
            continue
        x = "Signature requires a %s file, but one was not provided." % (
            sig.Normalization.upper())
        why_dropped[sig.xID] = x
    signatures = good

    # Process additional parameters from GenePattern.
    # o Do this before max_signatures, so that the maximum signatures
    #   is selected only out of the ones that the user specified.
    # o Do this before names and paths, so the variables will be
    #   aligned.
    # gp_imod_all_vars can be None or "".
    if options.gp_imod_all_vars:
        x = process_gp_imod_all_vars(
            options.gp_imod_all_vars, signatures, why_dropped)
        signatures, why_dropped = x

    sys.stdout.flush()
    DATA_rma = DATA_mas5 = None
    if datafile_rma is not None:
        print "Reading RMA file."
        DATA_rma = arrayio.read(datafile_rma)
        DATA_rma = arrayio.convert(DATA_rma, to_format=arrayio.gct_format)
    if datafile_mas5 is not None:
        print "Reading MAS5 file."
        DATA_mas5 = arrayio.read(datafile_mas5)
        DATA_mas5 = arrayio.convert(DATA_mas5, to_format=arrayio.gct_format)
    # Don't handle the log.  Let pybinreg do it.
    # Make sure the data sets contain the same samples.  Align them if
    # necessary.
    if DATA_rma and DATA_mas5:
        assert DATA_rma.ncol() == DATA_mas5.ncol(), \
               "RMA/MAS5 data sets have different numbers of samples."
        if not binreg.are_cols_aligned(DATA_rma, DATA_mas5):
            x = binreg.align_cols(DATA_rma, DATA_mas5)
            DATA_rma, DATA_mas5 = x
        assert binreg.are_cols_aligned(DATA_rma, DATA_mas5)
    if DATA_rma:
        arrayio.gct_format.write(
            DATA_rma, open(file_layout.DATASET_RMA, 'w'))
    if DATA_mas5:
        arrayio.gct_format.write(
            DATA_mas5, open(file_layout.DATASET_MAS5, 'w'))

    # Figure out the names and paths for each signature.
    names = [None] * len(signatures)   # SIG19_AKT[_modified]
    paths = [None] * len(signatures)   # <path>/SIG19_AKT[_modified]
    for i, sig in enumerate(signatures):
        name = "SIG%02d_%s" % (sig.xID, _hash_name(sig.Name))
        # If the user has modified the signature from the default
        # parameters, then make a note of it.
        if getattr(sig, "Changed", False):
            name = "%s_modified" % name
        outpath = os.path.join(file_layout.OUTPATH, name)
        names[i] = name
        paths[i] = outpath

    if options.max_signatures is not None:
        signatures = signatures[:options.max_signatures]

    # Make a list of the jobs.
    jobs = []  # list of cmd, outpath, outfile
    for i, sig in enumerate(signatures):
        name, outpath = names[i], paths[i]
        #print "Generating signature %s [%d:%d]" % (
        #    name, i+1, len(signatures))
        #sys.stdout.flush()
        
        quantile_normalize = False
        assert sig.Quantile.upper() in ["YES", "NO"]
        if sig.Quantile.upper() == "YES":
            quantile_normalize = True
        shift_scale_normalize = False
        assert sig.Shift_Scale.upper() in ["YES", "NO"]
        if sig.Shift_Scale.upper() == "YES":
            shift_scale_normalize = True
        
        #outfile = os.path.join(files.outpath, "%s.out.txt" % name)
        outfile = os.path.join(outpath, "out.txt")

        if sig.Normalization.upper() == "RMA":
            datafile = file_layout.DATASET_RMA
            assert DATA_rma
        elif sig.Normalization.upper() == "MAS5":
            datafile = file_layout.DATASET_MAS5
            assert DATA_mas5
        else:
            raise AssertionError, "Unknown normalization."

        # If the entire analysis should be archived, then go ahead and
        # archive each of the pybinreg runs too.  This will prevent
        # large analyses from taking up too much disk space.  The
        # drawback is that the files that are archived are no longer
        # available for use here.  Hopefully this won't be a problem.
        cmd = make_pybinreg_cmd(
            options.pybinreg, options.python, options.binreg_path,
            options.matlab, options.arrayplot, options.povray,
            options.cluster, options.libpath,
            outpath, options.archive, sig.Genes, sig.Metagenes,
            quantile_normalize, shift_scale_normalize,
            sig.Train0, sig.Train1, datafile)
        x = cmd, outpath, outfile
        jobs.append(x)

    # Run each of the jobs.
    if options.num_procs < 1 or options.num_procs > 100:
        parser.error("Please specify between 1 and 100 processes.")
    if options.num_procs > 1:
        if parallel.find():
            num_sigs = min(options.num_procs, len(jobs))
            print "Predicting %d signatures at a time." % num_sigs
        else:
            print("I could not find GNU parallel.  "
                  "Predicting 1 signature at a time.")
            options.num_procs = 1
        sys.stdout.flush()

    DEBUG = False   # Can disable pybinreg temporarily for debugging.
    if not DEBUG:  
        if options.num_procs <= 1:
            for x in jobs:
                cmd, outpath, outfile = x
                run_one_pybinreg(cmd, outpath, outfile)
        else:
            run_many_pybinreg(jobs, options.num_procs)

    if signatures:
        print "Extracting the reports from each signature."
        report_files = extract_reports(names, paths, file_layout)
        
        print "Combining probabilities from each of the signatures."
        summarize_probabilities(signatures, names, paths, file_layout)

        print "Making heatmap of the results."
        sys.stdout.flush()
        summarize_heatmap(
            options.python, options.arrayplot, options.cluster,
            options.libpath, file_layout)

        print "Summarizing signatures."
        summarize_signatures(signatures, file_layout)

        print "Making a report."
        summarize_report(
            signatures, orig_signatures, report_files, start_time, why_dropped,
            file_layout)

    if options.archive:
        print "Compressing results."
        sys.stdout.flush()
        archive.zip_path(file_layout.ATTIC)
        for i, sig in enumerate(signatures):
            name, outpath = names[i], paths[i]
            archive.zip_path(outpath)
    
    print "Done."
    
if __name__ == '__main__':
    #print "'%s'" % os.environ["PATH"]
    #print os.getcwd()
    #sys.exit(0)
    main()
