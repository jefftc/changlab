#!/usr/bin/env python

# TODO:
# - handle "tmp..." files

def main():
    import os
    import sys
    import time
    import argparse
    from genomicode import parselib
    from Betsy import config
    from Betsy import rule_engine_bie3
    from Betsy import module_utils as mlib
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", "--verbose", action="store_true")
    args = parser.parse_args()


    output_path = config.OUTPUTPATH
    if not os.path.exists(output_path):
        return
    output_path = os.path.realpath(output_path)

    print "BETSY cache path: %s" % output_path
    print

    x = os.listdir(output_path)
    x = [x for x in x if not x.startswith("tmp")]
    x = [os.path.join(output_path, x) for x in x]
    x = [x for x in x if os.path.isdir(x)]
    analysis_paths = x

    # list of (full_path, start_time, module_name, version, hash, parameters)
    analyses = []
    for path in analysis_paths:
        p, f = os.path.split(path)
        # index_bam_folder__B006__617b92ee4d313bcd0148b1ab6a91b12f
        x = f.split("__")
        assert len(x) == 3, path
        module_name, version, hash_ = x

        # Read the parameter file.
        x = os.path.join(path, rule_engine_bie3.BETSY_PARAMETER_FILE)
        if not os.path.exists(x):
            # Broken analysis.
            continue
        params = rule_engine_bie3._read_parameter_file(x)
        assert params.get("module_name", module_name) == module_name

        start_time = params.get("start_time")
        assert start_time, "Missing: start_time"
        start_time = time.strptime(start_time, rule_engine_bie3.TIME_FMT)
        x = path, start_time, module_name, version, hash_, params
        analyses.append(x)
        
    # Sort by start time
    schwartz = [(x[1], x) for x in analyses]
    schwartz.sort()
    analyses = [x[-1] for x in schwartz]

    sizes = [None]*len(analyses)  # parallel to analyses
    for i, x in enumerate(analyses):
        path, start_time, module_name, version, hash_, params = x

        metadata = params.get("metadata", {})

        # Format the time.
        time_str = time.strftime("%a %m/%d %I:%M %p", start_time)
        run_time = params.get("elapsed_pretty")
        assert run_time, "Missing: elapsed_pretty"
        if run_time == "instant":
            x = "ran instantly"
        else:
            x = "took %s" % run_time

        # Format the directory size.
        size = mlib.get_dirsize(path)
        size_str = parselib.pretty_filesize(size)
        sizes[i] = size
            
        print "[%s]  %s (%s; %s) %s" % (
            time_str, module_name, size_str, run_time, hash_)

        tool = metadata.get("tool")
        if tool:
            print "  %s" % tool
        commands = metadata.get("commands")
        if commands and args.verbose:
            for x in commands:
                print "  %s" % x
        sys.stdout.flush()

    # BUG: Does not account for size in tmp directories.
    total_size = sum(sizes)
    x = parselib.pretty_filesize(total_size)
    print "Total size: %s" % x

    x = os.statvfs(output_path)
    free_size = x.f_bavail*x.f_frsize
    x = parselib.pretty_filesize(free_size)
    print "Free space: %s" % x
    


if __name__ == '__main__':
    main()
