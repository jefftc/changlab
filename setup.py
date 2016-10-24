import os
from distutils.core import setup
from distutils.core import Extension


# XXX make sure numpy is installed first.
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


# Rlib/
#Betsy/dot_betsyrc
#genomicode/dot_genomicoderc


setup(
    name="Changlab Software",
    version="161023",
    description="Software developed in the Jeffrey Chang lab.",
    author="Chang Lab",
    author_email="jeffrey.t.chang@uth.tmc.edu",
    # XXX Betsy not in right place.
    packages=["arrayio", "genomicode", "Betsy"],
    package_dir={
        "Betsy" : "Betsy/Betsy",
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
                  include_dirs=[numpy_include],
                  ),
        Extension("Betsy.cbie3", ["Betsy/Betsy/cbie3module.c"]),
        ],
    scripts=[
        "scripts/align_matrices.py",
        "scripts/analyze_clinical_outcome.py",
        "scripts/analyze_phenotype.py",
        "scripts/annotate_geneset.py",
        "scripts/annotate_matrix.py",
        "scripts/annotate_nearby_genes.py",
        "scripts/arrayplot.py",
        "scripts/arrayplot2.py",
        "scripts/beeswarmplot.py",
        "scripts/bfrmfactor.py",
        "scripts/bfrmnorm.py",
        "scripts/bfrmproject.py",
        "scripts/boxplot.py",
        "scripts/calc_venn.py",
        "scripts/combatnorm.py",
        "scripts/convert_cel_cc1_to_v3.py",
        "scripts/download_tcga.py",
        "scripts/dwdnorm.py",
        "scripts/find_diffexp_genes.py",
        "scripts/find_subtypes.py",
        "scripts/gmt2matrix.py",
        "scripts/gsea.py",
        "scripts/histoplot.py",
        "scripts/lineplot.py",
        "scripts/pcaplot.py",
        "scripts/preprocess.py",
        "scripts/pybinreg.py",
        "scripts/pypeakseq.py",
        "scripts/pyrseqc.py",
        "scripts/pyspp.py",
        "scripts/run_genepattern.py",
        "scripts/scatterplot.py",
        "scripts/score_geneset.py",
        "scripts/scoresig.py",
        "scripts/slice_annot.py",
        "scripts/slice_genesets.py",
        "scripts/slice_matrix.py",
        "scripts/tfplot.py",
        "Betsy/scripts/betsy_run.py",
        "Betsy/scripts/betsy_manage_cache.py",
        ]
    )

