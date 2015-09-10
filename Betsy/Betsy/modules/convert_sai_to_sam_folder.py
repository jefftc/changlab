from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import config
        from genomicode import filelib
        from Betsy import module_utils

        fastq_node, sai_node, index_node, group_node = antecedents
        module_utils.safe_mkdir(out_path)

        # Technically, doesn't need the SampleGroupFile, since that's
        # already reflected in the sai data.  But better, because the
        # sai data might not always be generated by BETSY.

        fastq_path = fastq_node.identifier
        sai_path = sai_node.identifier
        index_path = index_node.identifier

        assert os.path.exists(fastq_path)
        assert os.path.exists(sai_path)
        assert os.path.exists(index_path)
        assert os.path.isdir(fastq_path)
        assert os.path.isdir(sai_path)
        assert os.path.isdir(index_path)

        # Find the merged fastq files.
        x = module_utils.find_merged_fastq_files(
            group_node.identifier, fastq_path)
        fastq_files = x

        # Find the sai files.
        x = filelib.list_files_in_path(sai_path)
        x = [x for x in x if x.lower().endswith(".sai")]
        sai_filenames = x
        assert sai_filenames, "No .sai files."

        # Find the indexed reference genome at:
        # <index_path>/<assembly>.fa
        x = os.listdir(index_path)
        x = [x for x in x if x.lower().endswith(".fa")]
        assert len(x) == 1, "Cannot find bwa index."
        x = x[0]
        reference_fa = os.path.join(index_path, x)
        
        bwa = module_utils.which_assert(config.bwa)
        # bwa samse <reference.fa> <input.sai> <input.fq> > <output.sam>
        # bwa sampe <reference.fa> <input_1.sai> <input_2.sai>
        #   <input_1.fq> <input_2.fq> > <output.sam>

        # list of (pair1.fq, pair1.sai, pair2.fq, pair2.sai, output.sam)
        # all full paths
        jobs = []
        for x in fastq_files:
            sample, pair1_fq, pair2_fq = x

            # The sai file should be in the format:
            # <sai_path>/<sample>.sai    Single end read
            # <sai_path>/<sample>_1.sai  Paired end read
            # <sai_path>/<sample>_2.sai  Paired end read
            # Look for pair1_sai and pair2_sai.
            pair1_sai = pair2_sai = None
            for sai_filename in sai_filenames:
                p, f = os.path.split(sai_filename)
                s, e = os.path.splitext(f)
                assert e == ".sai"
                if s == sample:
                    assert not pair1_sai
                    pair1_sai = sai_filename
                elif s == "%s_1" % (sample):
                    assert not pair1_sai
                    pair1_sai = sai_filename
                elif s == "%s_2" % (sample):
                    assert not pair2_sai
                    pair2_sai = sai_filename
            assert pair1_sai, "Missing .sai file: %s" % sample
            if pair2_fq:
                assert pair2_sai, "Missing .sai file 2: %s" % sample
            if pair2_sai:
                assert not pair2_fq, "Extra .sai file 2: %s" % sample
                
            sam_filename = os.path.join(out_path, "%s.sam" % sample)

            x = sample, pair1_fq, pair1_sai, pair2_fq, pair2_sai, sam_filename
            jobs.append(x)

        # Make a list of bwa commands.
        sq = module_utils.shellquote
        commands = []
        for x in jobs:
            sample, pair1_fq, pair1_sai, pair2_fq, pair2_sai, sam_filename = x

            samse = "samse"
            if pair2_fq:
                samse = "sampe"

            x = [
                bwa,
                samse,
                sq(reference_fa),
                sq(pair1_sai),
                sq(pair1_fq),
                ]
            if pair2_fq:
                x += [
                    sq(pair2_sai),
                    sq(pair2_fq),
                    ]
            x += [
                ">",
                sam_filename,
                ]
            x = " ".join(x)
            commands.append(x)
            
        module_utils.run_parallel(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        for x in jobs:
            sample, pair1_fq, pair1_sai, pair2_fq, pair2_sai, sam_filename = x
            assert module_utils.exists_nz(sam_filename), \
                   "Missing: %s" % sam_filename
    

    def name_outfile(self, antecedents, user_options):
        return "sam.bwa"
        #from Betsy import module_utils
        #original_file = module_utils.get_inputid(antecedents.identifier)
        #return 'bamFiles_' + original_file
