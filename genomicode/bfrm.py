"""

Functions:
find_bfrm_bin
find_bfrm_normalize
find_bfrm_project

get_default_parameters
init_parameters

read_model
clean_model
read_clean_model

get_affy_row_name
is_affx

"""
import os, sys

def find_bfrm_bin(default_path=None):
    import config
    
    search_paths = [
        default_path,
        config.bfrm, 
        "/home/jchang/projects/bfrm/data/bfrm_20_4/bfrm64",
        ]
    path = None
    for spath in search_paths:
        assert path is None
        if spath is None or not os.path.exists(spath):
            continue
        path = spath
        break
    return path

def find_bfrm_normalize(default_path=None):
    import config
    
    search_paths = [
        default_path,
        config.bfrm_normalize_path,
        "/home/jchang/projects/genepattern/data/BFRM_normalize",
        "/home/jefftc/projects/genepattern/BFRM_normalize",
        "/home/jefftc/projects/genepattern/data/BFRM_normalize",
        ]
    path = None
    for spath in search_paths:
        assert path is None
        if spath is None or not os.path.exists(spath):
            continue
        # Test for some BFRM files.
        files = [
            "bfrm64", "dataframe.jar", "computeCorrected.m", "setup.m",
            ]
        complete = True   # Is this distribution complete.
        for file in files:
            filename = os.path.join(spath, file)
            if not os.path.exists(filename):
                complete = False
                # Don't break here, so we can diagnose missing files.
                #break
        if not complete:
            continue
        path = spath
        break
    return path

def find_bfrm_project(default_path=None):
    import config
    
    search_paths = [
        default_path,
        config.bfrm_project_path,
        "/home/jchang/projects/genepattern/data/BFRM_project",
        "/home/jefftc/projects/genepattern/BFRM_project",
        "/home/jefftc/projects/genepattern/data/BFRM_project",
        ]
    path = None
    for spath in search_paths:
        assert path is None
        if spath is None or not os.path.exists(spath):
            continue
        # Test for some BFRM files.
        files = [
            "bfrm/getFacScores.m",
            ]
        complete = True   # Is this distribution complete.
        for file in files:
            filename = os.path.join(spath, file)
            if not os.path.exists(filename):
                complete = False
                # Don't break here, so we can diagnose missing files.
                #break
        if not complete:
            continue
        path = spath
        break
    return path

def get_default_parameters(bfrm_bin=None):
    import subprocess
    import tempfile
    import stat
    
    bfrm_bin = find_bfrm_bin(bfrm_bin)
    assert bfrm_bin and os.path.exists(bfrm_bin), "I could not find BFRM"

    temp_path = None
    bfrm_file = param_file = None
    cwd = None
    try:
        # Create the parameters in a temporary directory, to make sure
        # we don't overwrite anything.
        temp_path = tempfile.mkdtemp(dir=".")
        temp_path = os.path.realpath(temp_path)
        bfrm_file = os.path.join(temp_path, "bfrm64")
        param_file = os.path.join(temp_path, "default.parameters.txt")
        x = open(bfrm_bin, 'rb').read()
        open(bfrm_file, 'wb').write(x)
        os.chmod(bfrm_file, stat.S_IRUSR|stat.S_IWUSR|stat.S_IXUSR)
        
        # Run BFRM.
        cwd = os.getcwd()
        os.chdir(temp_path)
        cmd = "%s -default" % bfrm_file
        p = subprocess.Popen(
            cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
        w, r = p.stdin, p.stdout
        w.close()
        r.read()

        # Read the parameter file.
        assert os.path.exists(param_file), "BFRM failed."
        parameter_str = open(param_file).read()
    finally:
        if cwd:
            os.chdir(cwd)
        if bfrm_file and os.path.exists(bfrm_file):
            os.unlink(bfrm_file)
        if param_file and os.path.exists(param_file):
            os.unlink(param_file)
        if temp_path:
            os.rmdir(temp_path)
    return parameter_str

def set_params_dataset(parameters, NObservations, NVariables, DataFile):
    x = parameters
    x = _set_param(x, "NObservations", NObservations)
    x = _set_param(x, "NVariables", NVariables)
    x = _set_param(x, "DataFile", DataFile)
    return x

def set_params_latent_factors(parameters, num_factors):
    x = parameters
    x = _set_param(x, "NLatentFactors", num_factors)
    return x

def set_params_control(parameters, NControlVariables, HFile):
    x = parameters
    x = _set_param(x, "NControlVariables", NControlVariables)
    x = _set_param(x, "HFile", HFile)
    return x

def set_params_evol(parameters, evol_file, max_factors, max_genes):
    assert os.path.exists(evol_file), "Missing: %s" % evol_file
    x = open(evol_file).read().split()
    num_genes = len(x)

    x = parameters
    x = _set_param(x, "Evol", 1)
    x = _set_param(x, "EvolVarIn", num_genes)
    x = _set_param(x, "EvolVarInFile", evol_file)
    if max_factors:
        x = _set_param(x, "EvolMaximumFactors", max_factors)
    if max_genes:
        x = _set_param(x, "EvolMaximumVariables", max_genes)
    return x

def _get_param(parameters, name, default=None):
    import re

    m = re.search(
        r"\b%s\b *= *(.*?) *?$" % name, parameters,
        flags=re.IGNORECASE|re.MULTILINE)
    if not m:
        return default
    return m.group(1)

def _set_param(parameters, name, value):
    import re

    value = str(value).strip()

    x = parameters
    assert _get_param(x, name) is not None, "Missing parameter: %s" % name
    # If the name is commented, then remove the comment.
    x = x.replace("#%s" % name, "%s" % name)
    x = re.sub(
        r"(\b%s\b *= *).*? *?$" % name, r"\g<1>%s" % value, x,
        flags=re.IGNORECASE|re.MULTILINE)
    return x

def _open_file_or_archive(path_or_archive, name):
    # Return None if it does not exist.
    import zipfile
    import archive
    
    if path_or_archive.lower().endswith(".zip"):
        s2f = archive.unzip_dict(path_or_archive)
        if name not in s2f:
            return None
        zfile = zipfile.ZipFile(path_or_archive)
        return zfile.open(s2f[name])
    filename = os.path.join(path_or_archive, name)
    if not os.path.exists(filename):
        return None
    return open(filename)

def read_model(root, param_file=None):
    # Read the files from a BFRM analysis at root.  root can either be
    # a path to the BFRM output, or it can be a zipped BFRM model.
    # Does no extra processing on the contents of the files.
    # param_file should be a local filename, relative to root.
    # 
    # Return a dictionary of the contents of the files.  The
    # dictionary contains:
    # B              pxq          mA.txt
    # Bnz            pxq          mBz.txt
    # PostPib        pxq          mPostPib.txt
    # F              qxn          mF.txt
    # Psi            p            mPsi.txt
    # Tau            q            mTau.txt
    # VariablesIn    p            mVariablesIn.txt    Only for evol searches.
    # ExternalProb   (m-p)x(k+1)  mExternalProb.txt   Only in later BFRM versn.
    # NVariables
    # NObservations
    # NDesigns
    # NControlVariables
    # DataFile
    # 
    # p  num_genes_in_model
    # n  num_samples
    # q  num_factors              Includes intercept, control, designs, latent.
    # k  num_factors_no_design    Just intercept, latent, no design or control.
    # m  num_genes_in_dataset
    import arrayio
    
    tdf = arrayio.tab_delimited_format

    assert os.path.exists(root)
    param_file = param_file or "parameters.txt"
        
    opj = os.path.join
    # variable_name, filename, is_matrix
    outfiles = [
        ("A", "mA.txt", 1),
        ("Bnz", "mBz.txt", 1),
        ("PostPib", "mPostPib.txt", 1),
        ("F", "mF.txt", 1),
        ("Psi", "mPsi.txt", 0),
        ("Tau", "mTau.txt", 0),
        ("VariablesIn", "mVariablesIn.txt", 0),
        ("ExternalProb", "mExternalProb.txt", 1),
        ]
    model = {}
    for x in outfiles:
        name, file, is_matrix = x

        handle = _open_file_or_archive(root, file)
        # Some files might not exist.
        if not handle:
            continue
        
        if is_matrix:
            d = tdf.read(handle, 0, 0)
        else:
            # Should be list.
            type_fn = float
            if name == "VariablesIn":
                type_fn = int
            d = map(type_fn, handle.read().split())
        handle.close()
        model[name] = d

    handle = _open_file_or_archive(root, param_file)
    assert handle, "Missing parameter file: %s." % param_file
    params = handle.read()
    handle.close()

    num_designs = _get_param(params, "NDesigns")
    if num_designs is None:
        num_designs = _get_param(params, "NDesignVariables")
    assert num_designs is not None, "I could not find the NDesignVariables."
    model["NDesigns"] = int(num_designs)

    x = _get_param(params, "NControlVariables", 0)
    model["NControlVariables"] = int(x)
    model["NVariables"] = int(_get_param(params, "NVariables"))
    model["NObservations"] = int(_get_param(params, "NObservations"))
    model["DataFile"] = _get_param(params, "DataFile")

    # Figure out the dimensions of this model.
    p = len(model["Psi"])
    n = model["NObservations"]
    q = len(model["Tau"])
    m = model["NVariables"]
    k = q - model["NControlVariables"] - model["NDesigns"]

    # Check the dimensions of the matrices.
    assert model["A"].dim() == (p, q)
    assert model["Bnz"].dim() == (p, q)
    assert model["PostPib"].dim() == (p, q)
    assert model["F"].dim() == (q, n)
    assert len(model["Psi"]) == p
    assert len(model["Tau"]) == q
    if "VariablesIn" in model:
        assert len(model["VariablesIn"]) == p
    if "ExternalProb" in model:
        assert model["ExternalProb"].dim() == (m-p, k+1)

    return model
        
def clean_model(model, factor_cutoff=None):
    # Process the model for simpler handling.  Return a dictionary of
    # the cleaned up model.  The dictionary contains:
    # A              pxk
    # Bnz            pxk
    # PostPib        pxk
    # factors        pxk    1/0, based on factor_cutoff
    # ExternalProb   mxk    Parallel to data set.  (May be missing.)
    # F              kxn
    # Psi            p
    # Tau            k
    # VariablesIn    p      Indexes (0-based) of genes, relative to data set.
    # 
    # GENE_O         p      For sorting original model to this order.  0-based
    # FACTOR_O       k      For sorting original model to this order.  0-based
    # 
    # p  num_genes_in_model
    # n  num_samples
    # k  num_factors              Just latent factors.
    # m  num_genes_in_dataset
    #
    # Matrix variables, A, F, PostPib, factors, VariablesIn,
    # etc. are all sorted according to FACTOR_O and GENE_O.  FACTOR_O
    # and GENE_O are provided to help convert the original data set to
    # the same order as these variables.
    #
    # Changes from the original model:
    # o Remove the designs and controls.
    # o Make VariablesIn, if it doesn't exist.
    # o Make VariablesIn 0-based, instead of 1-based.
    # o Sort the genes and factors.
    from genomicode import jmath
    from genomicode import Matrix

    factor_cutoff = factor_cutoff or 0.99

    # Figure out the dimensions of this model.
    p = len(model["Psi"])
    n = model["NObservations"]
    q = len(model["Tau"])
    m = model["NVariables"]
    k = q - model["NControlVariables"] - model["NDesigns"]
                                
    # Get rid of the design variables.
    # Boundary case: k == 0, otherwise won't slice correctly.
    if k > 0:
        A = model["A"].matrix(None, (-k, None))
        Bnz = model["Bnz"].matrix(None, (-k, None))
        PostPib = model["PostPib"].matrix(None, (-k, None))
        F = model["F"].matrix((-k, None), None)
        Tau = model["Tau"][-k:]
    else:
        A = model["A"].matrix(None, [])
        Bnz = model["Bnz"].matrix(None, [])
        PostPib = model["PostPib"].matrix(None, [])
        F = model["F"].matrix([], None)
        Tau = []

    # Make the factors variable by thresholding PostPib on factor_cutoff.
    factors = PostPib.matrix()
    for x in factors._X:
        for i in range(len(x)):
            if x[i] >= factor_cutoff:
                x[i] = 1
            else:
                x[i] = 0

    # If VariablesIn doesn't exist, make it.
    VariablesIn = model.get("VariablesIn")
    if not VariablesIn:
        # First, make a 1-based index.
        VariablesIn = list(range(1, p+1))
    # Convert VariablesIn to 0-based index.
    VariablesIn = [x-1 for x in VariablesIn]

    # Make the ExternalProb variable.
    ExternalProb = None
    if "ExternalProb" in model:
        # First column is the (1-based) index of the genes.  These
        # should not overlap with VariablesIn.
        seen = {}.fromkeys(VariablesIn)
        for index in model["ExternalProb"].value(None, 0):
            index = int(index)-1
            assert index not in seen
            seen[index] = 1
        assert len(seen) == m
        
        prob = [[0]*k for i in range(m)]

        # Set the probabilities from PostPib.
        for i, index in enumerate(VariablesIn):
            prob[index] = PostPib.value(i, None)
        for x in model["ExternalProb"]._X:
            index, x = int(x[0])-1, x[1:]
            assert len(x) == k
            prob[index] = x
        #for index in model["mVariablesIn"]:
        #    print "%d\tmVariablesIn" % index
        #for index in model["mExternalProb"].value(None, 0):
        #    index = int(index)
        #    print "%d\tmExternalProb" % index
        #sys.exit(0)
        ExternalProb = Matrix.InMemoryMatrix(prob)

    # Order the factors based on decreasing number of genes.
    if factors._X:
        sums = jmath.mysum(factors._X, byrow=0)
        FACTOR_O = jmath.order(sums, decreasing=1)
    else:
        # No factors.
        FACTOR_O = []

    # Order the genes based on decreasing number of factors.  Earlier
    # factors should get much higher weights.
    X = factors.slice(None, FACTOR_O)
    if X:
        weights = [2**x for x in reversed(range(k))]
        sums = [None] * p
        for i in range(p):
            sums[i] = sum([x1*x2 for (x1, x2) in zip(X[i], weights)])
        GENE_O = jmath.order(sums, decreasing=1)
    else:
        # No factors.
        # Can't specify genes, or the indexes will be out of range of
        # an empty matrix.
        #GENE_O = list(range(p))
        GENE_O = []
    
    cmod = {}
    cmod["A"] = A.matrix(GENE_O, FACTOR_O)
    cmod["Bnz"] = Bnz.matrix(GENE_O, FACTOR_O)
    cmod["PostPib"] = PostPib.matrix(GENE_O, FACTOR_O)
    cmod["factors"] = factors.matrix(GENE_O, FACTOR_O)
    if ExternalProb:
        cmod["ExternalProb"] = ExternalProb.matrix(None, FACTOR_O)
    cmod["F"] = F.matrix(FACTOR_O, None)
    cmod["Psi"] = [model["Psi"][i] for i in GENE_O]
    cmod["Tau"] = [Tau[i] for i in FACTOR_O]
    cmod["VariablesIn"] = [VariablesIn[i] for i in GENE_O]
    cmod["GENE_O"] = GENE_O
    cmod["FACTOR_O"] = FACTOR_O

    # Check the dimensions of the matrices.
    if k:
        # If no factors, then the number of rows of the matrices won't
        # be preserved.
        assert cmod["A"].dim() == (p, k), "%s %s" % (A.dim(), (p, k))
        assert cmod["Bnz"].dim() == (p, k)
        assert cmod["PostPib"].dim() == (p, k)
        assert cmod["factors"].dim() == (p, k)
        if "ExternalProb" in cmod:
            assert cmod["ExternalProb"].dim() == (m, k)
        assert cmod["F"].dim() == (k, n)
        # Length of all these will be 0.
        assert len(cmod["Psi"]) == p
        assert len(cmod["Tau"]) == k
        assert len(cmod["VariablesIn"]) == p
        assert len(cmod["GENE_O"]) == p
        assert len(cmod["FACTOR_O"]) == k
    
    return cmod

def read_clean_model(root, param_file=None, factor_cutoff=None):
    model = read_model(root, param_file=param_file)
    model = clean_model(model, factor_cutoff=factor_cutoff)
    return model

def get_affy_row_name(MATRIX):
    affx_name = None
    for name in MATRIX.row_names():
        ids = MATRIX.row_names(name)
        affx_ids = [x for x in ids if is_affx(x)]
        if not affx_ids:
            continue
        affx_name = name
        break
    assert affx_name, "I could not find any affymetrix IDs."
    return affx_name
    
def is_affx(id):
    id = id.upper()
    return id.startswith("AFFX-") and id.endswith("_AT")

