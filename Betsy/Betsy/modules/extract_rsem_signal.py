from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, in_data, out_attributes, user_options, num_cores,
        out_filename):
        import os
        from genomicode import jmath
        from genomicode import AnnotationMatrix
        from genomicode import alignlib
        from Betsy import module_utils as mlib

        rsem_path = in_data.identifier
        assert os.path.exists(rsem_path)
        assert os.path.isdir(rsem_path)

        result_files = alignlib.find_rsem_result_files(rsem_path)
        assert result_files, "No .results files found."

        preprocess = out_attributes.get("preprocess")
        assert preprocess in ["tpm", "fpkm"]

        x = mlib.get_user_option(
            user_options, "genes_or_isoforms", not_empty=True,
            allowed_values=["genes", "isoforms"])
        get_genes = x == "genes"

        # For each of the gene files, get the expression data.
        sample2matrix = {}  # sample -> AnnotationMatrix
        for x in result_files:
            sample, gene_filename, isoform_filename = x
            # Get the gene results.
            # TODO: Implement isoforms.
            filename = gene_filename
            if not get_genes:
                filename = isoform_filename
            if filename is None:
                continue
            assert os.path.exists(filename)
            matrix = AnnotationMatrix.read(filename)
            # Do some checking on the matrix.
            assert "gene_id" in matrix.headers
            assert "transcript_id(s)" in matrix.headers
            assert "TPM" in matrix.headers
            assert "FPKM" in matrix.headers
            sample2matrix[sample] = matrix
        assert sample2matrix, "No samples"

        gene_id = transcript_id = None
        # Pull out the gene and transcript IDs.
        for matrix in sample2matrix.itervalues():
            x1 = matrix["gene_id"]
            x2 = matrix["transcript_id(s)"]
            if gene_id is None:
                gene_id = x1
            if transcript_id is None:
                transcript_id = x2
            assert x1 == gene_id
            assert x2 == transcript_id
        assert gene_id
        assert transcript_id
        assert len(gene_id) == len(transcript_id)

        # Assemble into a gene expression matrix.
        header = "TPM"
        if preprocess == "fpkm":
            header = "FPKM"
        t_data = []  # matrix, where each row is a sample.
        t_data.append(gene_id)
        t_data.append(transcript_id)
        samples = []
        for sample in sorted(sample2matrix):
            matrix = sample2matrix[sample]
            exp = matrix[header]
            assert len(exp) == len(gene_id)
            t_data.append(exp)
            samples.append(sample)

        data = jmath.transpose(t_data)
        header = ["gene_id", "transcript_id(s)"] + samples
        data = [header] + data

        # Write out the data file.
        handle = open(out_filename, 'w')
        for x in data:
            print >>handle, "\t".join(map(str, x))


    def name_outfile(self, antecedents, user_options):
        #from Betsy import module_utils
        #data_node, group_node = antecedents
        #original_file = module_utils.get_inputid(data_node.identifier)
        #filename = 'Samfolder_' + original_file
        #return filename
        return "alignments.bowtie"


## def concatenate_files(input_files, outfile):
##     """Concatenate multiple files into a single one"""
##     import os
##     import gzip
    
##     with open(outfile, 'w') as outfile:
##         for fname in input_files:
##             if fname.endswith('.gz'):
##                 if os.path.exists(fname):
##                     infile = gzip.open(fname)
##                 else:
##                     # see if the unzipped file exists
##                     fname = os.path.splitext(fname)[0]
##                     if os.path.exists(fname):
##                         infile = open(fname)
##                     else:
##                         raise ValueError('cannot find %s' % fname)
##             else:
##                 assert os.path.exists(fname), 'cannot find %s' % fname
##                 infile = open(fname)
##             for line in infile:
##                 outfile.write(line)


## def concatenate_multiple_line(group_dict, foldername):
##     """given group_dict, concatenate multiple line fastq files under foldername
##        output is a new dictionary in format <sample:[R1_sample,R2_sample]>
##                    or <sample:[sample]>"""
##     import os

##     current_dir = os.getcwd()
##     new_group_dict = {}
##     for sample_name, files in group_dict.iteritems():
##         if len(files) == 1:
##             outfile = os.path.join(current_dir, sample_name + '.fastq')
##             inputfiles = [os.path.join(foldername, x) for x in files[0]]
##             concatenate_files(inputfiles, outfile)
##             new_group_dict[sample_name] = [outfile]
##         elif len(files) == 2:
##             outfile_left = os.path.join(current_dir, sample_name + '_R1.fastq')
##             outfile_right = os.path.join(current_dir, sample_name + '_R2.fastq')
##             inputfiles_left = [os.path.join(foldername, x) for x in files[0]]
##             inputfiles_right = [os.path.join(foldername, x) for x in files[1]]
##             concatenate_files(inputfiles_left, outfile_left)
##             concatenate_files(inputfiles_right, outfile_right)
##             new_group_dict[sample_name] = [outfile_left, outfile_right]
##         else:
##             raise ValueError(
##                 'the number of files for each sample is not correct')
##     return new_group_dict


## def preprocess_multiple_sample(folder, group_dict, outfile, ref):
##     """folder:the folder where files are stored,
##        group_dict: A dictionary in format <sample:[[R1_samples],[R2_samples]]>
##                    or <sample:[[samples]]>
##        outfile: output file name,
##        ref: reference species, human or mouse"""
##     import os
##     import subprocess
##     from genomicode import config

##     if not os.path.exists(outfile):
##         os.mkdir(outfile)
##     if ref == 'human':
##         ref_file = config.rna_human
##     elif ref == 'mouse':
##         ref_file = config.rna_mouse
##     else:
##         raise ValueError("we cannot handle %s" % ref)
##     new_group_dict = concatenate_multiple_line(group_dict, folder)
##     for sample in new_group_dict:
##         files = new_group_dict[sample]
##         if len(files) == 1:
##             input_file = os.path.join(folder, files[0])
##             command = ['bowtie', '-q', '--phred33-quals', '-n', '2', '-e',
##                        '99999999', '-l', '25', '-p', '8', '-a', '-m', '200',
##                        '-S', ref_file, input_file]
##         elif len(files) == 2:
##             input_file1 = os.path.join(folder, files[0])
##             input_file2 = os.path.join(folder, files[1])
##             command = ['bowtie', '-q', '--phred33-quals', '-n', '2', '-e',
##                        '99999999', '-l', '25', '-I', '1', '-X', '1000', '-p',
##                        '8', '-a', '-m', '200', '-S', ref_file, '-1',
##                        input_file1, '-2', input_file2]
##         else:
##             raise ValueError('number files is not correct')
##         outfilename = os.path.join(outfile, sample + '.sam')
##         f = file(outfilename, 'w')
##         try:
##             process = subprocess.Popen(
##                 command, shell=False, stdout=f, stderr=subprocess.PIPE)
##             process.wait()
##             error_message = process.communicate()[1]
##             if 'error' in error_message:
##                 raise ValueError(error_message)
##         finally:
##             f.close()
