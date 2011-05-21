#!/usr/bin/env python

import os

def main():
    from optparse import OptionParser, OptionGroup
    
    from genomicode import affyio

    usage = "usage: %prog [options] <cel_file>"
    parser = OptionParser(usage=usage, version="%prog 01")
    parser.add_option(
        "--no-clobber",  dest="clobber", action="store_false", default=True)
    
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.error("Please specify an infile.")
    filename, = args
    assert os.path.exists(filename), "I could not find file: %s" % filename

    version = affyio.guess_cel_version(filename)
    assert version == "cc1", "File does not appear to be cc1: %s" % filename

    outfile = "%s.v3" % filename
    assert options.clobber or not os.path.exists(outfile), "noclobber"
    affyio.convert_cel_cc1_to_3(filename, open(outfile, 'w'))

if __name__ == '__main__':
    main()
