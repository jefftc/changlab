#!/usr/bin/env python

# Some important paths:
# /usr/lib64/python2.7/              library_path
#   os.py
#   site-packages/                   site_packages_path
#     genomicode/
#       Rlib/                        genomicode_rlib_path
#       bin/                         genomicode_bin_path
#     Betsy/
#     arrayio/
#     numpy/core/include/numpy/      numpy_include_path
#       arrayobject.h
# $HOME/                             home_path.  User's home directory.
#   .genomicoderc
#   .betsyrc


# Functions:
# get_home_path
# get_genomicode_rlib_path
# get_genomicode_bin_path
# get_site_packages_path
# get_numpy_include_path
#
# get_original_user
# username2uid
# get_original_uid
#
# read_config_file
# merge_config_files
# setup_dot_betsyrc
# setup_dot_genomicoderc
# autoset_betsyrc_value
# autoset_genomicoderc_value
# 
# check_required_software


# Configuration parameters that point to executable files.  This
# enables this setup script to search for them in defined places,
# e.g. /usr/local/bin/.
CONFIG_BINARY_FILES = [
    # General Applications
    "python",
    "matlab",
    "povray",
    "Rscript",

    # Bioinformatics Software
    "primer3",
    "patser",
    "cluster",
    "blat",

    # Chang Lab
    "slice_matrix.py",
    "annotate_matrix.py",
    "align_matrices.py",
    "preprocess.py",
    "bfrmnorm.py",
    "combatnorm.py",
    "analyze_clinical_outcome.py",
    "find_diffexp_genes.py",
    "analyze_phenotype.py",
    "pybinreg.py",
    "scoresig.py",
    "score_geneset.py",
    "gsea.py",
    "pyrseqc.py",
    "pypeakseq.py",
    "pyspp.py",
    "xls2txt",
    "txt2xls",
    "run_genepattern.py",
    "arrayplot2.py",
    "download_tcga.py",

    # NGS Software
    "samtools",
    "bcftools",
    "vcfutils.pl",
    "vcftools",
    "bgzip",
    "tabix",
    "bedtools",
    "gtfutils",
    "bam2fastx",

    # NGS Alignment
    "bowtie-build",
    "bowtie",
    "bowtie2-build",
    "bowtie2",
    "bwa",

    # RNA-SEQ
    "tophat",
    "rsem-prepare-reference",
    "rsem-calculate-expression",
    "STAR",
    "htseq-count",

    # NGS Variant Calling
    "Platypus.py",
    "bam-somaticsniper",
    "MuSE",
    "table_annovar.pl",

    # [NGS ChIP-Seq]
    "macs14",
    "macs2",
    "PeakSeq",
    ]

# Configuration parameters that are not files.
CONFIG_NOT_FILE = [
    # Gene Pattern
    "gp_server",
    "gp_user",
    "gp_passwd",
    # Gene Map
    "gm_user",
    "gm_passwd",
    "gm_db",
    "gm_host",
    "gm_port",
    # BETSY
    "CACHE_PATH",    # Hack.  This path doesn't have to exist yet.
    ]


def get_home_path():
    import os
    return os.environ["HOME"]


def get_genomicode_rlib_path():
    import os
    
    opj = os.path.join
    return opj(get_site_packages_path(), "genomicode", "Rlib")


def get_genomicode_bin_path():
    import os
    
    opj = os.path.join
    return opj(get_site_packages_path(), "genomicode", "bin")


def get_site_packages_path():
    import os
    x = os.path.split(os.__file__)[0]      # library_path
    x = os.path.join(x, "site-packages")   # site_path
    assert os.path.exists(x)
    return x


def get_numpy_include_path():
    # Figure out the path for numpy include files.
    import os
    import numpy

    # numpy.get_include() returns:
    # .../python2.7/site-packages/numpy/core/include
    # However, headers are in:
    # .../python2.7/site-packages/numpy/core/include/numpy
    # Fix this path.
    numpy_include = numpy.get_include()
    x = os.path.join(numpy_include, "arrayobject.h")
    if not os.path.exists(x):
        numpy_include = os.path.join(numpy_include, "numpy")
    x = os.path.join(numpy_include, "arrayobject.h")
    assert os.path.exists(x), "I could not find arrayobject.h in numpy."
    return numpy_include


def get_original_user():
    # Get the username.  If this is run with sudo, get the username
    # before changing with sudo.
    import os

    assert "USER" in os.environ
    user = os.environ["USER"]
    if "SUDO_USER" in os.environ:
        user = os.environ["SUDO_USER"]
    return user


def username2uid(username):
    import pwd
    return pwd.getpwnam(username).pw_uid


def get_original_uid():
    user = get_original_user()
    uid = username2uid(user)
    return uid
    

def find_bin(name):
    import os

    opj = os.path.join
    SEARCH_PATH = [
        opj(get_home_path(), "bin"),
        "/opt/local/bin",
        "/usr/local/bin",
        "/usr/bin",
        get_genomicode_bin_path(),
        ]
    path, name = os.path.split(name)
    if path:
        # If this is a full path name, then search the full path.
        SEARCH_PATH = [path] + SEARCH_PATH
    for path in SEARCH_PATH:
        filename = os.path.join(path, name)
        if os.path.exists(filename):
            return filename
    # Not found.
    return None


def cluster_find_bin(name):
    import os

    opj = os.path.join
    SEARCH_PATH = [
        opj(get_home_path(), "bin"),
        "/opt/local/bin",
        "/usr/local/bin",
        "/usr/bin",
        get_genomicode_bin_path(),
        ]
    path, name = os.path.split(name)
    if path:
        # If this is a full path name, then search the full path.
        SEARCH_PATH = [path] + SEARCH_PATH
    for path in SEARCH_PATH:
        filename = os.path.join(path, name)
        if not os.path.exists(filename):
            continue
        try:
            get_cluster_version(filename)
        except AssertionError, x:
            # Was not correct.
            continue
        # Found correct one.
        return filename
    # Not found.
    return None
    

def read_config_file(filename):
    # Read a configuration file and return a dictionary of:
    # <section> -> <name> -> <value>
    import ConfigParser
    
    # Read the configuration.
    config = ConfigParser.ConfigParser()
    config.optionxform = str   # use case sensitive option names
    # Will overwrite duplicated keys.
    config.read(filename)

    config_values = {}
    for section in config.sections():
        section_values = config_values.get(section, {})
        for (name, value) in config.items(section):
            section_values[name] = value
        config_values[section] = section_values
    return config_values


def _print_name_value(handle, name, delim, value):
    import os
    
    x = "%s%s%s" % (name, delim, value)
    if name not in CONFIG_NOT_FILE and not os.path.exists(value):
        x = "#%s" % x
    print >>handle, x
    

def merge_config_files(old_file, new_file, value_handler_fn=None):
    # Merge two configuration files.  Use the formatting of the
    # new_file.  Take values from old_file when available.

    # BUG: Comments from the old_file will be lost.
    import os
    import StringIO

    print "Generating configuration file: %s" % old_file
    old_values = {}
    if os.path.exists(old_file):
        old_values = read_config_file(old_file)
    merged_values = {}  # <section> -> <name> -> <value>

    # Make the merged configuration file.
    handle = StringIO.StringIO()
    section = None
    for line in open(new_file):
        # Line can be:
        # 1.  Blank
        # 2.  Comment
        #     Starts with "#" or ";".
        # 3.  Section
        #     [<section>]
        # 4.  Variable definition.
        #     <name>=<value>
        #     <name>:<value>
        
        line_s = line.strip()
        # Case 1.  Print out blank lines.
        if not line_s:
            print >>handle, line_s
            continue
        # Case 2.  Print out comment lines.
        elif line.startswith("#") or line.startswith(";"):
            print >>handle, line_s
            continue
        # Internal comments not supported.  Check after Case 2 to
        # allow:
        # ### Comment.
        assert "#" not in line_s and ";" not in line_s, \
               "Internal comments not supported"
        # Continuations not supported.  Check after Case 2 to allow:
        #    # Comment.
        assert not line.startswith(" "), "Continuations not supported"
        # Case 3.  Section.
        if line_s.startswith("[") and line_s.endswith("]"):
            # New section.
            if section is not None:
                # See if the old_file contains any variables in this
                # section that were not already written out.  If so,
                # then write them out now.
                for name, value in old_values.get(section, {}).iteritems():
                    if name in merged_values.get(section, {}):
                        continue
                    print "Betsy config file %s has unrecognized variable: %s"\
                          % (old_file, name)
                    _print_name_value(handle, name, "=", value)
                    merged_values[name] = value
            print >>handle, line_s
            x = line_s[1:-1]
            section = x.strip()
            assert section not in merged_values, \
                   "Duplicate section: %s" % section
            merged_values[section] = {}
            continue
        # Case 4.  Variable definition.
        # Each variable should be in a section.
        assert section is not None, "Missing section header"
        delim = "="
        if line_s.find("=") < 0:
            delim = ":"
        assert line_s.find(delim) >= 0, "Unknown format: %s" % line_s
        x = line_s.split(delim, 1)
        assert len(x) == 2
        name, value = x

        # If the old value already exists, then don't overwrite it.
        assert section in merged_values
        old_value = old_values.get(section, {}).get(name)
        if old_value is not None:
            merged_values[section][name] = old_value
            _print_name_value(handle, name, delim, old_value)
            continue
        # New value.  See if we can figure out a good default value
        # for it.
        if value_handler_fn is not None:
            value = value_handler_fn(name, value)
        merged_values[section][name] = value
        _print_name_value(handle, name, "=", value)

    # See if the old_file contains any variables in this
    # section that were not already written out.  If so,
    # then write them out now.
    if section is not None:
        for name, value in old_values.get(section, {}).iteritems():
            if name in merged_values.get(section, {}):
                continue
            print "Betsy config file %s has unrecognized variable: %s" % (
                old_file, name)
            _print_name_value(handle, name, "=", value)
            merged_values[name] = value

    # See if the old_file contains any sections that were not already
    # written out.
    for section in old_values:
        if section in merged_values:
            continue
        print "Betsy config file %s has unrecognized section: %s" % (
            old_file, section)
        print >>handle
        print >>handle, "[%s]" % section
        print >>handle
        for name, value in old_values[section].iteritems():
            _print_name_value(handle, name, "=", value)
            merged_values[name] = value
        
    handle.seek(0)
    open(old_file, 'w').write(handle.read())

    # If this script was run with sudo, the config_file will be owned
    # by "root".  Change back to the original user.
    uid = get_original_uid()
    os.chown(old_file, uid, -1)


def setup_dot_betsyrc():
    import os
    
    template_file = os.path.join("Betsy", "dot_betsyrc")
    config_file = os.path.join(get_home_path(), ".betsyrc")
    assert os.path.exists(template_file), "File not found: %s" % template_file

    # If the configuration file already exists, merge them.  This will
    # make sure the configuration includes all the values from
    # dot_betsyrc.
    merge_config_files(config_file, template_file)


def setup_dot_genomicoderc():
    import os
    
    template_file = os.path.join("genomicode", "dot_genomicoderc")
    config_file = os.path.join(get_home_path(), ".genomicoderc")
    assert os.path.exists(template_file), "File not found: %s" % template_file

    # If the configuration file already exists, merge them.  This will
    # make sure the configuration includes all the values from
    # dot_betsyrc.
    merge_config_files(
        config_file, template_file,
        value_handler_fn=autoset_genomicoderc_value)


def autoset_betsyrc_value(name, default_value):
    # See if we can figure out a better value for the default_value.
    import os

    if name == "CACHE_PATH":
        return os.path.join(get_home_path, "betsy.out")
    return default_value


def quote(s, always_quote=False):
    BAD_CHARS = " \\"

    if type(s) in [type(0), type(0.0)]:
        s = str(s)
    needs_quote = False
    if not always_quote:
        for x in BAD_CHARS:
            if x in s:
                needs_quote = True
                break
    if always_quote or needs_quote:
        s = "'" + s.replace("'", "'\\''") + "'"
    return s


def get_cluster_version(cluster_bin):
    # There's two programs called "cluster".  One is cluster 3.0.  The
    # other is from graphviz.
    
    # cluster --version

    # Output from cluster30:
    # Cluster 3.0, command line version (no GUI support),
    # using the C Clustering Library version 1.50.
    # [...]

    # Output from graphviz:
    # option -- unrecognized - ignored
    #Usage: cluster <options> graphfile
    #    -C k - generate no more than k clusters (0)
    #       0 : no limit
    #    -c k - use clustering method k (0)
    #       0 : use modularity
    #       1 : use modularity quality
    #    -o <outfile> - output file (stdout)
    #    -v   - verbose mode
    #    -?   - print usage

    import re
    import subprocess

    cmd = "%s --version" % quote(cluster_bin)
    x = subprocess.check_output(
        cmd, stderr=subprocess.STDOUT, shell=True)
    x = x.strip()
    if x.find("cluster <options> graphfile") >= 0:
        raise AssertionError, "Found graphviz cluster, not Cluster 3.0"
    m = re.search(r"Cluster ([\w\. ]+)", x)
    assert m, "Missing version string"
    return m.group(1)


def autoset_genomicoderc_value(name, default_value):
    # See if we can figure out a better value for the default_value.
    import os
    
    if not default_value.strip():
        return default_value
    path, file_ = os.path.split(default_value)

    # Handle "cluster" separately.  Can be either GraphViz or
    # Cluster3.0.  Make sure we don't automatically set it to
    # GraphViz's cluster.
    if file_ == "cluster":
        filename = cluster_find_bin(default_value)
        if filename:
            default_value = filename

    if file_ in CONFIG_BINARY_FILES:
        filename = find_bin(default_value)
        if filename:
            default_value = filename
    elif file_ == "Rlib":
        x = get_genomicode_rlib_path()
        if os.path.exists(x):
            default_value = x
    elif file_ == "ComBat.R":
        x = get_genomicode_rlib_path()
        x = os.path.join(x, file_)
        if os.path.exists(x):
            default_value = x
            
    return default_value
    

def check_required_software():
    try:
        import numpy
    except ImportError, x:
        raise AssertionError, \
              "Please make sure numpy is available before installing."


def main():
    import glob
    from distutils.core import setup
    from distutils.core import Extension
    from distutils.command.install_scripts import install_scripts
    from distutils.command.install import install

    # Check for required software.
    check_required_software()

    # Make a list of the scripts to install.
    changlab_scripts = glob.glob("scripts/*.py")
    betsy_scripts = glob.glob("Betsy/scripts/*.py")
    assert changlab_scripts
    assert betsy_scripts

    # Make a list of the Rlib/ files.  Set the install path.
    rlib_files = glob.glob("Rlib/*.R")
    rlib_path = get_genomicode_rlib_path()

    # Set up configuration files.
    class my_install(install):
        def run(self):
            install.run(self)
            setup_dot_betsyrc()
            setup_dot_genomicoderc()

    # Install scripts under genomicode.
    class my_install_scripts(install_scripts):
        def run(self):
            self.install_dir = get_genomicode_bin_path()
            install_scripts.run(self)

    setup(
        name="Changlab Software",
        version="1",
        author="Chang Lab",
        author_email="jeffrey.t.chang@uth.tmc.edu",
        description="Software developed in the Jeffrey Chang lab.",
        license="MIT",
        #python_requires=">= 2.7",

        packages=[
            "arrayio", "genomicode", "Betsy", "Betsy.modules", "Betsy.rules"],
        package_dir={
            "Betsy" : "Betsy/Betsy",
            "Betsy.modules" : "Betsy/Betsy/modules",
            "Betsy.rules" : "Betsy/Betsy/rules",
            },
        ext_modules=[
            Extension("genomicode.ciolib", ["genomicode/ciolibmodule.c"]),
            Extension("genomicode.cMatrix", ["genomicode/cMatrixmodule.c"]),
            Extension(
                "genomicode.cjmath",
                ["genomicode/cjmathmodule.c", "genomicode/mathlib.c"]),
            Extension("genomicode.cgeolib", ["genomicode/cgeolibmodule.c"]),
            Extension(
                "genomicode.ccoherencelib",
                ["genomicode/ccoherencelibmodule.c", "genomicode/mathlib.c"]),
            Extension("genomicode.cGeneSet", ["genomicode/cGeneSetmodule.c"]),
            Extension("genomicode.cMarkovModel", ["genomicode/cMarkovModel.c"],
                      include_dirs=[get_numpy_include_path()],
                      ),
            Extension("Betsy.cbie3", ["Betsy/Betsy/cbie3module.c"]),
            ],

        data_files=[
            (rlib_path, rlib_files),
            ],
        scripts=changlab_scripts + betsy_scripts,

        cmdclass = {
            "install_scripts" : my_install_scripts,
            "install" : my_install,
            }
        )
    

if __name__ == '__main__':
    main()
