#!/usr/bin/env python

import os, sys

def find_annotation_file(chipname):
    from genomicode import config

    chipname2annotfile = {
        "HG-U133A" : "HG-U133A_annot.csv.gz",
        "HG-U133A_2" : "HG-U133A_2.na30.annot.csv.gz",
        "HG-U133B" : "HG-U133B_annot.csv.gz",
        "HG-U133_Plus_2" : "HG-U133_Plus_2_annot.csv.gz",
        "Mouse430A_2" : "Mouse430A_2_annot.csv.gz",
        "RAE230A" : "RAE230A.na26.annot.csv.gz",
        "HT_HG-U133A" : "HT_HG-U133A.na31.annot.csv.gz",
        "HG_U95Av2" : "HG_U95Av2_annot.csv.gz",
        }
    assert chipname in chipname2annotfile, "I have no annotations for %s" % \
           chipname
    file = chipname2annotfile[chipname]
    filename = os.path.join(config.normalize_AFFYMETRIX, file)
    assert os.path.exists(filename), "I could not find annotation file %s." % \
           filename
    filename = os.path.realpath(filename)
    return filename

def find_normscript():
    from genomicode import config
    
    file = "normscript.R"
    filename = os.path.join(config.normalize_RLIB, file)
    assert os.path.exists(filename), "I could not find normscript.R."
    filename = os.path.realpath(filename)
    return filename

def main():
    import tempfile
    from optparse import OptionParser, OptionGroup
    from genomicode import affyio

    usage = (
        "usage: %prog [options] algorithm path_to_cel_files\n\n"
        'algorithm should be "RMA" or "MAS5".')
    parser = OptionParser(usage=usage, version="%prog 01")

    parser.add_option(
        "--platform", dest="platform", type="string", default=None,
        help="Normalize only arrays from this platform.")
    parser.add_option(
        "-s", "--filestem", dest="filestem", type="string", default=None,
        help="Use this as the filestem for the normalized matrix.")
    parser.add_option(
        "-n", "--noclobber", action="store_true", dest="noclobber",
        default=False)

    options, args = parser.parse_args()
    if len(args) != 2:
        parser.error("I expected 2 arguments.")
    algorithm, path = args
    algorithm = algorithm.upper()

    assert os.path.exists(path), "Path %s does not exist." % path
    assert algorithm.upper() in ["RMA", "MAS5"], \
           "Algorithm must be RMA or MAS5."

    if options.filestem:
        filestem = options.filestem
    else:
        # Guess the filestem based on the name of the path.
        x = os.path.split(path)[1]
        x = os.path.splitext(x)[0]
        x = x.replace("_cel", "")   # For ArrayExpress paths
        x = x.replace(".CEL", "")   # For ArrayExpress paths
        filestem = x
    
    outfile = "%s.%s" % (filestem, algorithm.lower())
    if options.noclobber and os.path.exists(outfile):
        print "Outfile %s exists.  Will not overwrite." % outfile
        return

    if algorithm == "MAS5":
        log_signal, filter_25, filter_50 = 0, 0, 0
    elif algorithm == "RMA":
        log_signal, filter_25, filter_50 = 0, 0, 0
    else:
        raise AssertionError, "Unknown algorithm: %s" % algorithm

    # Figure out the chip for each CEL file.
    file2chipname = {}
    for i, file in enumerate(os.listdir(path)):
        filename = os.path.join(path, file)
        cn = affyio.extract_chip_name(filename)
        assert cn, "Unknown chip type: %s" % filename
        file2chipname[file] = cn

    # If the user specified a specific platform, then ignore all files
    # not of that platform.
    for (file, cn) in file2chipname.items():
        if options.platform and cn != options.platform:
            del file2chipname[file]

    # Make sure files exist.
    assert file2chipname, "no files in path: %s" % path

    # Make sure all the files are the same type.
    # Should normalize all the chip files separately.
    chipnames = sorted({}.fromkeys(file2chipname.values()).keys())
    if len(chipnames) != 1:
        # Count the number of times each chip appears.
        counts = {}
        for x in file2chipname.values():
            counts[x] = counts.get(x, 0) + 1
        x = ["%s (%d)" % x for x in counts.iteritems()]
        x = ", ".join(x)
        assert len(chipnames) == 1, "multiple platforms: %s" % x
    chipname = chipnames[0]

    # Choose the annotation file for the chip name.
    annotfile = find_annotation_file(chipname)
    # Should allow the user to generate a file that's not annotated.
    assert annotfile, "I don't know the annotation file for %s" % chipname
    assert os.path.exists(annotfile), "Missing %s [%s]" % (annotfile, chipname)

    temppath = None
    try:
        # Make a directory with only the CEL files with this type of annotation
        temppath = tempfile.mkdtemp(dir=".")
        for file, cn in file2chipname.iteritems():
            if cn != chipname:
                continue
            file1 = os.path.join(os.path.realpath(path), file)
            file2 = os.path.join(temppath, file)
            os.symlink(file1, file2)
            #cmd = 'ln -s "%s" "%s"' % (file1, file2)
            #os.system(cmd)

        # Format the arguments and call the normalize script.
        normscript = find_normscript()
        x = (temppath, annotfile, filestem, algorithm,
             log_signal, filter_25, filter_50)
        x = " ".join(map(str, x))
        cmd = "cat %s | R --vanilla %s" % (normscript, x)
        print "NORM: %s\n" % cmd
        os.system(cmd)
    finally:
        # Clear out stuff from the temppath.
        if temppath and os.path.exists(temppath):
            for file in os.listdir(temppath):
                filename = os.path.join(temppath, file)
                os.unlink(filename)
            os.rmdir(temppath)

    print "Done"

if __name__ == '__main__':
    main()
