from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        from genomicode import parallel
        from genomicode import filelib
        from genomicode import alignlib
        from genomicode import hashlib
        from Betsy import module_utils as mlib

        fastq_node, sample_node, orient_node, reference_node = antecedents
        fastq_files = mlib.find_merged_fastq_files(
            sample_node.identifier, fastq_node.identifier)
        ref = alignlib.create_reference_genome(reference_node.identifier)
        assert os.path.exists(ref.fasta_file_full)
        orient = mlib.read_orientation(orient_node.identifier)
        filelib.safe_mkdir(out_path)

        metadata = {}
        metadata["tool"] = "bowtie2 %s" % alignlib.get_bowtie2_version()

        # Bowtie2 doesn't handle files with spaces in them.  Make
        # temporary files without spaces.
        

        # Make a list of the jobs to run.
        jobs = []
        for i, x in enumerate(fastq_files):
            sample, pair1, pair2 = x
            bam_filename = os.path.join(out_path, "%s.bam" % sample)
            log_filename = os.path.join(out_path, "%s.log" % sample)
            sample_h = hashlib.hash_var(sample)
            temp_pair1 = "%d_%s_1.fa" % (i, sample_h)
            temp_pair2 = None
            if pair2:
                temp_pair2 = "%d_%s_2.fa" % (i, sample_h)
            j = filelib.GenericObject(
                sample=sample,
                pair1=pair1,
                pair2=pair2,
                temp_pair1=temp_pair1,
                temp_pair2=temp_pair2,
                bam_filename=bam_filename,
                log_filename=log_filename)
            jobs.append(j)

        for j in jobs:
            os.symlink(j.pair1, j.temp_pair1)
            if pair2:
                os.symlink(j.pair2, j.temp_pair2)
        
        # Generate bowtie2 commands for each of the files.
        attr2orient = {
            "single" : None,
            "paired_fr" : "fr",
            "paired_rf" : "rf",
            "paired_ff" : "ff",
            }
        orientation = attr2orient[orient.orientation]
        #x = sample_node.data.attributes["orientation"]
        #orientation = attr2orient[x]

        # Takes ~4 Gb per job.
        samtools = mlib.findbin("samtools")
        sq = parallel.quote
        commands = []
        for j in jobs:
            #sample, pair1, pair2, bam_filename, log_filename = x
            nc = max(1, num_cores/len(jobs))

            # bowtie2 -p 8 -x <genome> -1 <.fq> -2 <.fq> --fr
            #  2> test.log | samtools view -bS -o test.bam -
            x1 = alignlib.make_bowtie2_command(
                ref.fasta_file_full, j.temp_pair1, fastq_file2=j.temp_pair2,
                orientation=orientation, num_threads=nc)
            x2 = [
                sq(samtools),
                "view",
                "-bS",
                "-o", sq(j.bam_filename),
                "-",
                ]
            x2 = " ".join(x2)
            x = "%s 2> %s | %s" % (x1, sq(j.log_filename), x2)
            #x = "%s >& %s" % (x, sq(log_filename))
            commands.append(x)
        metadata["commands"] = commands
        parallel.pshell(commands, max_procs=num_cores)

        # Make sure the analysis completed successfully.
        x1 = [x.bam_filename for x in jobs]
        x2 = [x.log_filename for x in jobs]
        filelib.assert_exists_nz_many(x1+x2)
            
        return metadata


    def name_outfile(self, antecedents, user_options):
        return "bowtie2.bam"
