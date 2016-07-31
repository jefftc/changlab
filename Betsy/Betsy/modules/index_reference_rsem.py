from Module import AbstractModule

class Module(AbstractModule):
    def __init__(self):
        AbstractModule.__init__(self)

    def run(
        self, network, antecedents, out_attributes, user_options, num_cores,
        out_path):
        import os
        import shutil
        #from genomicode import config
        from genomicode import filelib
        from genomicode import parallel
        from genomicode import alignlib
        from Betsy import module_utils as mlib

        ref_node, gene_node = antecedents
        # Don't copy the whole path.  Just get the fasta file.
        #ref = alignlib.standardize_reference_genome(
        #    ref_node.identifier, out_path, use_symlinks=True)
        ref = alignlib.create_reference_genome(ref_node.identifier)
        gtf_file = gene_node.identifier
        filelib.assert_exists_nz(gtf_file)

        # Symlink the fasta file into the out path.
        filelib.safe_mkdir(out_path)
        x = os.path.join(out_path, ref.fasta_file)
        os.symlink(ref.fasta_file_full, x)

        # rsem-prepare-reference --bowtie --bowtie2 --gtf gtf02.gtf
        #   <reference.fa> <reference_name>
        # <reference_name>.[1234].ebwt    # Bowtie1.
        # <reference_name>.rev.[12].ebwt
        # <reference_name>.[1234].bt2     # Bowtie2.
        # <reference_name>.rev.[12].bt2
        # <reference_name>.chrlist        # RSEM.
        # <reference_name>.grp
        # <reference_name>.idx.fa
        # <reference_name>.n2g.idx.fa
        # <reference_name>.seq
        # <reference_name>.ti
        # <reference_name>.transcripts.fa
        # chrLength.txt                   # STAR
        # chrNameLength.txt
        # chrName.txt
        # chrStart.txt
        # exonGeTrInfo.tab
        # exonInfo.tab
        # gencode.vM8.annotation.gtf
        # geneInfo.tab
        # Genome
        # genomeParameters.txt
        # SA
        # SAindex
        # sjdbInfo.txt
        # sjdbList.fromGTF.out.tab
        # sjdbList.out.tab
        # transcriptInfo.tab

        rsem_prepare = mlib.get_config("rsem_prepare", which_assert_file=True)
        bowtie = mlib.get_config("bowtie", which_assert_file=True)
        bowtie2 = mlib.get_config("bowtie2", which_assert_file=True)
        STAR = mlib.get_config("STAR", which_assert_file=True)

        # RSEM wants the path that contains the executables.
        bowtie = os.path.split(bowtie)[0]
        bowtie2 = os.path.split(bowtie2)[0]
        STAR = os.path.split(STAR)[0]

        sq = parallel.quote
        cmd = [
            sq(rsem_prepare),
            "--num-threads", num_cores,
            "--bowtie",
            "--bowtie-path", sq(bowtie),
            "--bowtie2",
            "--bowtie2-path", sq(bowtie2),
            "--star",
            "--star-path", sq(STAR),
            "--gtf", sq(gtf_file),
            sq(ref.fasta_file_full),
            ref.name,
            ]
        parallel.sshell(cmd, path=out_path)

        # Copy the GTF file into the output path.
        shutil.copy2(gtf_file, out_path)

        assembly = ref.name
        # Check to make sure index was created successfully.
        x1 = ["%s.%d.ebwt" % (assembly, i+1) for i in range(4)]
        x2 = ["%s.rev.%d.ebwt" % (assembly, i+1) for i in range(2)]
        x3 = ["%s.%d.bt2" % (assembly, i+1) for i in range(4)]
        x4 = ["%s.rev.%d.bt2" % (assembly, i+1) for i in range(2)]
        x5 = [
            "%s.chrlist" % assembly,
            "%s.grp" % assembly,
            "%s.idx.fa" % assembly,
            "%s.n2g.idx.fa" % assembly,
            "%s.seq" % assembly,
            "%s.ti" % assembly,
            "%s.transcripts.fa" % assembly,
            ]
        x6 = [
            "chrLength.txt",
            "chrNameLength.txt",
            "chrName.txt",
            "chrStart.txt",
            "exonGeTrInfo.tab",
            "exonInfo.tab",
            "gencode.vM8.annotation.gtf",
            "geneInfo.tab",
            "Genome",
            "genomeParameters.txt",
            "SA",
            "SAindex",
            "sjdbInfo.txt",
            "sjdbList.fromGTF.out.tab",
            "sjdbList.out.tab",
            "transcriptInfo.tab",
            ]
        x = x1 + x2 + x3 + x4 + x5 + x6
        index_files = [os.path.join(out_path, x) for x in x]
        filelib.assert_exists_nz_many(index_files)


    def name_outfile(self, antecedents, user_options):
        return "reference.rsem"
