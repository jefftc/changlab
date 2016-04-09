from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        import os
        from genomicode import filelib
        from genomicode import alignlib

        count_path = in_data.identifier
        assert os.path.exists(count_path)
        assert os.path.isdir(count_path)
        result_files = filelib.list_files_in_path(
            count_path, endswith=".count")
        assert result_files, "No .count files found."

        # Parse the count files.
        name2results = {}
        for filename in result_files:
            x = os.path.split(filename)[1]
            x = os.path.splitext(x)[0]
            name = x
            assert name not in name2results
            x = alignlib.parse_htseq_count_output(filename)
            name2results[name] = x
        assert name2results, "No samples"

        # Assemble into a summary matrix.
        # Rows:
        # no_feature
        # ambiguous
        # too_low_aQual
        # not_aligned
        # alignment_not_unique
        # total_mapped
        # total_fragments
        # percent_mapped

        ROWS = [
            "no_feature",
            "ambiguous",
            "too_low_aQual",
            "not_aligned",
            "alignment_not_unique",
            ]
        all_names = sorted(name2results)
        
        matrix = []
        header = ["Feature"] + all_names
        matrix.append(header)
        for rn in ROWS:
            x = [rn] + [getattr(name2results[n], rn) for n in all_names]
            assert len(x) == len(header)
            matrix.append(x)
            
        # Count the total mapped and total_fragments.
        total_mapped = []
        total_fragments = []
        perc_mapped = []
        perc_no_feature = []
        perc_ambiguous = []
        for n in all_names:
            # Sum up the counts
            results = name2results[n]
            tm, tf, pm = "", "", ""
            pnf, pamb = "", ""
            if not results.errors:
                x1 = sum(results.counts.values())
                x2 = 0
                for rn in ROWS:
                    x2 += getattr(results, rn)
                tm = x1
                tf = x1+x2
                pm = tm/float(tf)
                pnf = results.no_feature / float(tf)
                pamb = results.ambiguous / float(tf)
            total_mapped.append(tm)
            total_fragments.append(tf)
            perc_mapped.append(pm)
            perc_no_feature.append(pnf)
            perc_ambiguous.append(pamb)

        x1 = ["total_mapped"] + total_mapped
        x2 = ["total_fragments"] + total_fragments
        x3 = ["perc_mapped"] + perc_mapped
        x4 = ["perc_no_feature"] + perc_no_feature
        x5 = ["perc_ambiguous"] + perc_ambiguous
        assert len(x1) == len(header)
        assert len(x2) == len(header)
        assert len(x3) == len(header)
        assert len(x4) == len(header)
        assert len(x5) == len(header)
        matrix.append(map(str, x1))
        matrix.append(map(str, x2))
        matrix.append(map(str, x3))
        matrix.append(map(str, x4))
        matrix.append(map(str, x5))

        # Write the data file.
        handle = open(outfile, 'w')
        for x in matrix:
            print >>handle, "\t".join(map(str, x))


    def name_outfile(self, antecedents, user_options):
        return "summary.txt"
    
