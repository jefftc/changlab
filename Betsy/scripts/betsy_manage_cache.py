#!/usr/bin/env python

def main():
    import os
    import time
    from Betsy import config
    from Betsy import rule_engine_bie3
    
    output_path = config.OUTPUTPATH
    if not os.path.exists(output_path):
        return
    output_path = os.path.realpath(output_path)

    print "BETSY cache path: %s" % output_path

    x = os.listdir(output_path)
    x = [os.path.join(output_path, x) for x in x]
    x = [x for x in x if os.path.isdir(x)]
    analysis_paths = x

    # list of (start_time, module_name, version, hash, parameters)
    analyses = []
    for path in analysis_paths:
        p, f = os.path.split(path)
        # index_bam_folder__B006__617b92ee4d313bcd0148b1ab6a91b12f
        x = f.split("__")
        assert len(x) == 3
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
        x = start_time, module_name, version, hash_, params
        analyses.append(x)
    # Sort by start time
    analyses.sort()

    for x in analyses:
        start_time, module_name, version, hash_, params = x
        
        metadata = params.get("metadata", {})

        # Print out the module.
        time_str = time.strftime("%a %I:%M %p", start_time)
        run_time = params.get("elapsed_pretty")
        assert run_time, "Missing: elapsed_pretty"
        if run_time == "instant":
            x = "ran instantly"
        else:
            x = "took %s" % run_time
        print "[%s]  %s (%s) %s" % (time_str, module_name, run_time, hash_)

        tool = metadata.get("tool")
        if tool:
            print "  %s" % tool
        commands = metadata.get("commands")
        if 0 and commands:
            for x in commands:
                print "  %s" % x

        


if __name__ == '__main__':
    main()
