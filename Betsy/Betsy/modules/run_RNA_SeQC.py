from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_path):
        import os
        import subprocess
        from Betsy import module_utils
        from genomicode import config

        module_utils.safe_mkdir(out_path)
        in_path = in_data.identifier
        bam_filenames = module_utils.find_bam_files(in_path)

        # java -jar /usr/local/bin/RNA-SeQC_v1.1.8.jar \
        #   -o <sample> -r <reference_file> -s "<sample>|<in_filename>|NA"
        #   -t <gtf_file> >& <log_filename>"
        # <in_filename>     BAM file
        # <reference_file>  /data/biocore/genomes/UCSC/mm10.fa
        # <gtf_file>   /data/biocore/rsem/mouse_refseq_mm10/UCSC_knownGenes.gtf
        #
        # <reference_file> must be indexed and has a dict file.

        rna_seqc_jar = module_utils.which_assert(config.rna_seqc_jar)
        
        
        REF = user_options['RNA_ref']
        GTF = user_options['RNA_gtf']
        assert os.path.exists(REF), ('cannot find the %s' % REF)
        assert os.path.exists(GTF), ('cannot find the %s' % GTF)
        
        for filename in filenames:
            infile = os.path.join(directory, filename)
            outname = os.path.splitext(filename)[0]

            command = ['java', '-jar', RNA_SeQC_path, '-o', outfile, '-r', REF,
                       '-s', outname + '|' + infile + '|NA', '-t', GTF]
            process = subprocess.Popen(command,
                                       shell=False,
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
            process.wait()
            error_message = process.communicate()
            # guess if there is error message in the output
            if 'error' in error_message[1]:
                raise ValueError(error_message)
    
    
    def name_outfile(self, antecedents, user_options):
        from Betsy import module_utils
        original_file = module_utils.get_inputid(antecedents.identifier)
        filename = 'bamFolder_' + original_file
        return filename
