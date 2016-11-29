from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        outfile):
        #import subprocess
        import shutil
        import arrayio
        #from genomicode import config
        from genomicode import arrayplatformlib
        from genomicode import parallel
        #from genomicode import filelib
        from Betsy import module_utils as mlib

        DATA = arrayio.read(in_data.identifier)

        #chipname = arrayplatformlib.identify_platform_of_matrix(DATA)
        scores = arrayplatformlib.score_matrix(DATA)
        assert scores, "Unable to identify platform: %s" % in_data.identifier
        chipname = scores[0]

        platform = "HG_U133A"
        assert arrayplatformlib.get_bm_attribute(platform), \
               "Unrecognized platform: %s" % platform

        if chipname == platform:
            shutil.copyfile(in_data.identifier, outfile)
        else:
            Annot_BIN = mlib.get_config(
                "annotate_matrix", which_assert_file=True)
            sq = parallel.quote
            cmd = [
                "python",
                sq(Annot_BIN),
                sq(in_data.identifier),
                "--platform", sq(platform),
                '--min_match_score', 0.80,
                ]
            cmd = " ".join(map(str, cmd))
            cmd = "%s >& %s" % (cmd, sq(outfile))
            parallel.sshell(cmd)

            #f = file(outfile, 'w')
            #try:
            #    process = subprocess.Popen(
            #        command, shell=False, stdout=f, stderr=subprocess.PIPE)
            #finally:
            #    f.close()
            #error_message = process.communicate()[1]
            #if error_message:
            #    raise ValueError(error_message)
        #change the HG_U133A to the first column

        f = file(outfile, 'r')
        txt = f.readlines()
        f.close()
        header = txt[0].split('\t')
        index = header.index('HG_U133A')
        f = file(outfile, 'w')
        for line in txt:
            line = line.split('\t')
            newline = [line[index]] + line[0:index] + line[index + 1:]
            f.write('\t'.join(newline))
        f.close()


    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #filename = 'signal_' + original_file + '.tdf'
        #return filename
        return "signal.tdf"

