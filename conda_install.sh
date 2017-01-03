# Contributed by Brad Chapman (https://github.com/chapmanb);
# chapmanb@fastmail.com.


# Installing dependencies using Miniconda.

wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -b -p ./install

#wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
#bash Miniconda2-latest-MacOSX-x86_64.sh -b -p ./install


# Base system
./install/bin/conda install -c conda-forge -c r -c bioconda \
  numpy matplotlib pygraphviz rpy2 xlwt openpyxl xlrd pil pandas
./install/bin/python setup.py install

# NGS utilities -- missing gtfutils
./install/bin/conda install -c conda-forge -c r -c bioconda \
  samtools bcftools picard vcftools bedtools tophat

# NGS alignment/processing
./install/bin/conda install -c conda-forge -c r -c bioconda \
  trimmomatic fastqc bowtie bowtie2 bwa star

# NGS RNA-seq
./install/bin/conda install -c conda-forge -c r -c bioconda \
  tophat rsem htseq rseqc

# NGS variant and somatic calling -- missing several
./install/bin/conda install -c conda-forge -c r -c bioconda \
  gatk platypus-variant varscan snpeff pindel
