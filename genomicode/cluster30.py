"""

Functions:
cluster30
cluster_hierarchical
cluster_kmeans

cut_dendrogram          Cut a tree into clusters.

"""

# _cluster
# 
# _choose_gene_id
# _choose_gene_label
#
# _find_cluster_files
# _read_cluster_files
# _write_cluster_files    DEPRECATED?
# _unlink_cluster_files
# 
# _guess_filestem
# _convert_to_pcl
# _exists_nz


DIST2ID = {
    "uncent-cor"  : 1,  "pearson"    : 2,  "abs-uncent-cor" : 3,
    "abs-pearson" : 4,  "spearman"   : 5,  "kendall"        : 6,
    "euclidean"   : 7,  "city-block" : 8,
    }
METHOD2ID = {
    "complete" : "m", "single" : "s", "centroid" : "c", "average" : "a",
    }


class ClusterData:
    def __init__(
        self, matrix, gene_tree, array_tree,
        gene_tree_cluster, array_tree_cluster, gene_cluster, array_cluster):
        # matrix   Matrix object
        # gtr      gene_tree            hierarchical
        # atr      array_tree           hierarchical
        # gtc      gene_tree_cluster    hierarchical
        # atc      array_tree_cluster   hierarchical
        # kgg      gene_cluster         kmeans
        # kag      array_cluster        kmeans
        # These can be None, if missing.
        self.matrix = matrix
        self.gene_tree = gene_tree
        self.array_tree = array_tree
        self.gene_tree_cluster = gene_tree_cluster
        self.array_tree_cluster = array_tree_cluster
        self.gene_cluster = gene_cluster
        self.array_cluster = array_cluster


def cluster30_file(
    filename, cluster_genes, cluster_arrays, algorithm, 
    distance=None, method=None, kmeans_k=None, som_rows=None, som_cols=None,
    cluster_bin=None, jobname=None):
    # Run cluster3.0.  Return the command used to run it.
    import config
    import parallel
    import filelib
    
    assert algorithm in ["hierarchical", "kmeans", "som"]
    distance = distance or "pearson"
    assert distance in DIST2ID, "Unknown distance: %s" % distance
    if algorithm == "hierarchical":
        method = method or "average"
    if method:
        assert method in METHOD2ID, "Unknown method: %s" % method
    assert not (method and algorithm in ["kmeans", "som"])
    assert not (kmeans_k and algorithm != "kmeans")
    assert not (som_rows and algorithm != "som")
    assert not (som_cols and algorithm != "som")
    assert cluster_genes or cluster_arrays
    
    assert not kmeans_k or (kmeans_k > 1 and kmeans_k <= 100)
    assert not som_rows or (som_rows >= 1 and som_rows <= 100)
    assert not som_cols or (som_cols >= 1 and som_cols <= 100)

    cluster = cluster_bin or config.cluster or "cluster"
    cluster = filelib.which_assert(cluster)

    array_metric = 0  # no clustering
    gene_metric = 0
    if cluster_arrays:
        array_metric = DIST2ID[distance]
    if cluster_genes:
        gene_metric = DIST2ID[distance]
        
    sq = parallel.quote
    cmd = [
        sq(cluster),
        "-f", sq(filename),
        ]
    if jobname:
        cmd += ["-u", jobname]
    if method is not None:
        cmd += ["-m", METHOD2ID[method]]
    if kmeans_k is not None:
        cmd += ["-k", kmeans_k]
    if algorithm == "som":
        cmd += ["-s"]
    if som_rows:
        cmd += ["-y", som_rows]
    if som_cols:
        cmd += ["-x", som_cols]
    cmd += [
        "-e", array_metric,
        "-g", gene_metric,
        ]
    cmd = " ".join(map(str, cmd))
    parallel.sshell(cmd)
    return cmd
    #return _read_cluster_files(jobname)


def cluster30(MATRIX, cluster_genes, cluster_arrays,
              algorithm, distance=None, method=None, kmeans_k=None,
              cluster_bin=None):
    # Cluster a Matrix and return a ClusterData object.
    # 
    # Arguments:
    # cluster_genes   boolean, whether to cluster the genes
    # cluster_arrays  boolean, whether to cluster the arrays
    # algorithm       "hierarchical", "kmeans"
    # distance        See code (hierarchical or kmeans).
    # method          See code (hierarchical only)
    # kmeans_k        Number of clusters (kmeans only)

    # UNUSED?
    # gene_k          Number of gene clusters (hierarchical only).
    # array_k         Number of array clusters (hierarchical only).
    
    assert algorithm in ["hierarchical", "kmeans"]
    distance = distance or "pearson"
    assert distance in DIST2ID, "Unknown distance: %s" % distance
    if algorithm == "hierarchical":
        method = method or "average"
    if method:
        assert method in METHOD2ID, "Unknown method: %s" % method
    assert not (algorithm == "kmeans" and method)
    assert not (algorithm == "hierarchical" and kmeans_k)
    assert cluster_genes or cluster_arrays
    if algorithm == "kmeans":
        assert kmeans_k and kmeans_k > 1
    else:
        assert not kmeans_k
    
    args = []
    id_ = DIST2ID[distance]
    if cluster_genes:
        args.append("-g %s" % id_)
    else:
        args.append("-g 0")
    if cluster_arrays:
        args.append("-e %s" % id_)
    else:
        args.append("-e 0")

    if method:
        id_ = METHOD2ID[method]
        args.append("-m %s" % id_)

    if algorithm == "kmeans":
        assert kmeans_k
        args.append("-k %d" % kmeans_k)

    cluster_data = _cluster(MATRIX, cluster=cluster_bin, *args)

    ## # Cluster the hierarchical trees, if necessary.
    ## gene_tree_cluster = array_tree_cluster = None
    ## # If I haven't reclustered the data, then the old tree is still
    ## # valid.
    ## if not cluster_genes:
    ##     gene_tree_cluster = cluster_data.gene_tree_cluster
    ## if not cluster_arrays:
    ##     array_tree_cluster = cluster_data.array_tree_cluster
    ## if cluster_data.gene_tree and gene_k:
    ##     assert gene_k <= MATRIX.nrow(), "more gene clusters than genes"
    ##     gene_tree_cluster = clusterio.cut_dendrogram(
    ##         cluster_data.gene_tree, gene_k)
    ## if cluster_data.array_tree and array_k:
    ##     assert array_k <= MATRIX.ncol(), "more array clusters than arrays"
    ##     array_tree_cluster = clusterio.cut_dendrogram(
    ##         cluster_data.array_tree, array_k)
    ## cluster_data.gene_tree_cluster = gene_tree_cluster
    ## cluster_data.array_tree_cluster = array_tree_cluster
    return cluster_data


## def cluster_matrix(
##     MATRIX, cluster, cluster_data,
##     cluster_genes, cluster_arrays, algorithm, distance, method,
##     gene_k, array_k, kmeans_k):
##     from genomicode import clusterio

##     assert algorithm in ["hierarchical", "kmeans"]

##     dist2id = {
##         "uncent-cor"  : 1,  "pearson"    : 2,  "abs-uncent-cor" : 3,
##         "abs-pearson" : 4,  "spearman"   : 5,  "kendall"        : 6,
##         "euclidean"   : 7,  "city-block" : 8,
##         }
##     method2id = {
##         "complete" : "m", "single" : "s", "centroid" : "c", "average" : "a",
##         }

##     # Skip if all conditions are true:
##     # - not clustering genes
##     # - not clustering arrays
##     # - not cutting gene tree (and gene tree already exists)
##     # - not cutting array tree (and array tree already exists)
##     if (not cluster_genes and not cluster_arrays and
##         not (gene_k and cluster_data.gene_tree) and 
##         not (array_k and cluster_data.array_tree)):
##         return MATRIX, cluster_data

##     # If not clustering and just re-cutting the tree, then don't
##     # bother regenerating the clusters.
##     if cluster_genes or cluster_arrays:
##         args = []

##         id_ = dist2id[distance]
##         if cluster_genes:
##             args.append("-g %s" % id_)
##         else:
##             args.append("-g 0")
##         if cluster_arrays:
##             args.append("-e %s" % id_)
##         else:
##             args.append("-e 0")

##         id_ = method2id[method]
##         args.append("-m %s" % id_)

##         if algorithm == "kmeans":
##             args.append("-k %d" % kmeans_k)
            
##         filestem = _cluster(MATRIX, cluster=cluster, *args)
##         files = find_data_files(filestem)
##         #print filestem, files
##         assert "cdt" in files, "No cdt file produced."
##         MATRIX, cluster_data = read_data_set(filestem, cluster_data)
##         _cleanup_cluster(filestem)
    
##     # Cluster the hierarchical trees, if necessary.
##     gene_tree_cluster = array_tree_cluster = None
##     # If I haven't reclustered the data, then the old tree is still
##     # valid.
##     if not cluster_genes:
##         gene_tree_cluster = cluster_data.gene_tree_cluster
##     if not cluster_arrays:
##         array_tree_cluster = cluster_data.array_tree_cluster
##     if cluster_data.gene_tree and gene_k:
##         assert gene_k <= MATRIX.nrow(), "more gene clusters than genes"
##         gene_tree_cluster = clusterio.cut_dendrogram(
##             cluster_data.gene_tree, gene_k)
##     if cluster_data.array_tree and array_k:
##         assert array_k <= MATRIX.ncol(), "more array clusters than arrays"
##         array_tree_cluster = clusterio.cut_dendrogram(
##             cluster_data.array_tree, array_k)
##     cluster_data.gene_tree_cluster = gene_tree_cluster
##     cluster_data.array_tree_cluster = array_tree_cluster

##     return MATRIX, cluster_data


def cluster_hierarchical(MATRIX, cluster_genes, cluster_arrays,
                         distance=None, method=None, cluster_bin=None):
    return cluster30(
        MATRIX, cluster_genes, cluster_arrays, "hierarchical",
        distance=distance, method=method, cluster_bin=cluster_bin)


def cluster_kmeans(MATRIX, cluster_genes, cluster_arrays,
                   distance=None, kmeans_k=None, cluster_bin=None):
    return cluster30(
        MATRIX, cluster_genes, cluster_arrays, "kmeans",
        distance=distance, kmeans_k=kmeans_k, cluster_bin=cluster_bin)


def cut_dendrogram(tree, k):
    # Use a BFS to cut the tree into k clusters.  Return a dictionary
    # of id -> cluster.
    assert len(tree)
    
    tree_cluster = {}  # id -> cluster
    depth = 1          # The root is at depth 1.
    next_cluster = 0

    # Do a BFS through the tree, assigning clusters.
    stack = []         # node id, dist, cluster_num
    
    # Assume the root node is the last element of the tree and add it
    # to the stack.
    stack.append((-len(tree), tree[-1][2], None))
    while stack:
        # Sort the nodes from highest (closest to root) to lowest.
        # Small numbers are high nodes, big numbers are low nodes.
        stack = sorted(stack, cmp=lambda x, y: cmp(x[1], y[1]))

        # To do a BFS, get the highest node.
        id, dist, num = stack.pop(0)

        # If deep enough, and there's no cluster, create a new one.

        # If this node is not already in a cluster, and I'm deep
        # enough down the tree (for the number of clusters wanted), or
        # this node is a leaf (which is automatically deep enough),
        # then assign it to a cluster.
        if num is None and (depth >= k or id >= 0):
            num = next_cluster
            next_cluster += 1
            assert next_cluster <= k, "more clusters than items"

        # Assign the cluster.
        tree_cluster[id] = num

        if id >= 0:
            # If this is a leaf, then do nothing more.
            continue

        # Add the left and right nodes to the stack.
        left_id, right_id, dist = tree[-id-1]
        left_num = right_num = num
        left_dist = right_dist = 0
        if left_id < 0:
            # If this is an internal node, then set the distance.  The
            # default distance for leaves are 0, so they'll be
            # processed before internal nodes.
            left_dist = tree[-left_id-1][2]
        if right_id < 0:
            right_dist = tree[-right_id-1][2]
        stack.append((left_id, left_dist, left_num))
        stack.append((right_id, right_dist, right_num))

        # I've gone down one node, so increase the depth.
        depth += 1

    return tree_cluster

##     # The clustering reordered the items in the expression data set.
##     # Order cluster so that it matches the order in the expression
##     # data.
##     # item_order  item_ids      oindex   cindex
##     # GENE12X     200648_s_at     12       0
##     # GENE5X      200604_s_at      5       1
##     # GENE20X     200722_s_at     20       2
##     # [...]
##     #
##     # item2cluster is based on oindex.  Need to return the clusters in
##     # cindex order.

##     # Calculate the original indexes.
##     oindex = [_parse_gene_or_node(x) for x in item_order]

##     assert len(item2cluster) == len(tree)+1
##     cluster = [None] * (len(item2cluster))
##     for oi, n in item2cluster.iteritems():
##         ci = oindex.index(oi)
##         cluster[ci] = item_ids[oi], n

##     return cluster


def _cluster(MATRIX, *args, **params):
    # Return a ClusterData object.
    #
    # Parameters:
    # cluster   path for "cluster" binary (optional)
    import os
    import tempfile
    import subprocess
    import arrayio
    import filelib
    import config

    path = "."
    cluster = params.get("cluster") or config.cluster or "cluster"

    filestem = pcl_file = None
    try:
        x, filestem = tempfile.mkstemp(dir=path); os.close(x)
        filelib.safe_unlink(filestem)

        # Write the data set in PCL format.
        # This implementation requires a matrix in PCL format.
        MATRIX = _convert_to_pcl(MATRIX)
        pcl_file = filestem + ".pcl"
        arrayio.pcl_format.write(MATRIX, open(pcl_file, 'w'))

        args = list(args)
        args.append('-f "%s"' % pcl_file)

        cmd = "%s %s" % (cluster, " ".join(args))
        #print cmd
        #w, r = os.popen4(cmd)
        p = subprocess.Popen(
            cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, close_fds=True)
        w, r = p.stdin, p.stdout
        w.close()
        output = r.read()
        #print output
        if output.find("cluster: command not found") >= 0:
            raise AssertionError, "cluster: command not found"
        elif output.find("cluster: No such file or directory") >= 0:
            raise AssertionError, output.strip()
        elif output.find("command not found") >= 0:
            raise AssertionError, "%s: command not found" % cluster
        elif output.find("Error reading file") >= 0:
            raise AssertionError, "%s\n%s" % (cmd, output)
        elif output.find("cluster <options> graphfile") >= 0:
            raise AssertionError, "ran Graphviz cluster, not Cluster 3.0"

        files = _find_cluster_files(filestem)
        #print filestem, files
        assert "cdt" in files, "No cdt file produced."
        cluster_data = _read_cluster_files(filestem)
    finally:
        if pcl_file:
            os.unlink(pcl_file)
        if filestem:
            _unlink_cluster_files(filestem)

    return cluster_data


def _choose_gene_id(MATRIX):
    # Given a user-specified matrix, try to pick a good unique ID for
    # the genes.
    import arrayio
    
    headers = MATRIX.row_names()

    # Prioritize some potential ones.  Don't use the standard headers,
    # e.g. arrayio.ROW_ID, so that we can preserve the user's header.
    IDS = ["Probe.Set.ID", "Probe Set ID", "Probe ID"]
    for id_ in IDS:
        if id_ in headers:
            return id_

    # If no known headers are found, then choose a standard one.
    IDS = [arrayio.AFFY_PROBESET_ID, arrayio.GENE_ID, arrayio.ROW_ID]
    for id_ in IDS:
        if id_ in headers:
            return id_

    # If no standard ones are found, then just arbitrarily use the
    # first column that is not missing any values.
    for header in headers:
        names = MATRIX.row_names(header)
        missing = [x for x in names if not x.strip()]
        if not missing:
            return header
    
    raise AssertionError, "I could not find an ID for the matrix."


def _choose_gene_label(MATRIX):
    import arrayio

    names = MATRIX.row_names()

    # Prioritize some potential ones.
    IDS = [
        arrayio.GENE_SYMBOL, "Gene.Symbol", "Gene Symbol", "Symbol",
        #arrayio.GENE_DESCRIPTION, "Description",
        "DESCRIPTION",       # For GCT files.  Use the pretty name.
        "NAME",
        arrayio.GENE_ID, "LocusLink",
        arrayio.AFFY_PROBESET_ID, "Probe.Set.ID",
        arrayio.ROW_ID
        ]
    # Exception: If the GCT files have generic descriptions,
    # e.g. DESC0001, then use the name field instead.
    if "DESCRIPTION" in names:
        desc = MATRIX.row_names("DESCRIPTION")
        if desc[0].startswith("DESC"):
            i = IDS.index("DESCRIPTION")
            IDS.pop(i)
    for id_ in IDS:
        if id_ in names:
            return id_
    if names:
        return names[0]
    raise AssertionError, "I could not find an ID for the matrix."


def _find_cluster_files(file_or_stem):
    # Return a dictionary of extension -> filename.
    import os

    # Needs to get the stem, because it can be hard to find the file
    # names based on the CDT file, because the stem is different.
    # Example files:
    # <stem>.gtc
    # <stem>_K_A5.kag
    # <stem>_K_G5.kgg
    # <stem>_K_G5_A5.cdt    Hard to cut into a stem.
    # <stem>_SOM_G2-2_A2-2.txt
    # <stem>_SOM_G2-2_A2-2.anf
    # <stem>_SOM_G2-2_A2-2.gnf
    fullstem = _guess_filestem(file_or_stem)

    # Bug: should do a case insensitive search.
    path, stem = os.path.split(fullstem)
    if not path:
        path = "."
    EXTENSIONS = [
        "nrm", "pcl", "cdt", "gtr", "atr", "gtc", "atc", "kgg", "kag",
        "txt", "anf", "gnf"]
    ext2file = {}
    for file_ in os.listdir(path):
        if not file_.startswith(stem):
            continue
        f, e = os.path.splitext(file_)
        if e.startswith("."):
            e = e[1:]
        if e not in EXTENSIONS:
            continue
        
        # Do some checking to make sure file looks reasonable.
        recognize_file = False
        if f == stem:
            recognize_file = True
        elif f.startswith("%s_K_A" % stem):
            recognize_file = True
        elif f.startswith("%s_K_G" % stem):
            recognize_file = True
        elif f.startswith("%s_SOM_" % stem):
            recognize_file = True
        if not recognize_file:
            continue
        
        ext2file[e] = os.path.join(path, file_)
    return ext2file
    

def _read_cluster_files(file_or_stem):
    # Return a ClusterData object.
    import os
    import arrayio
    from genomicode import parselib
    from genomicode import clusterio

    files = _find_cluster_files(file_or_stem)
    #print "FOUND", files; sys.exit(0)
    
    filename = file_or_stem
    if not os.path.exists(filename):
        # If this file does not exist, then look for a CDT, NRM, or
        # PCL file (in that order).
        DATA_EXTS = ["cdt", "nrm", "pcl"]
        DATA_EXTS = [x for x in DATA_EXTS if x in files]
        assert DATA_EXTS, "I could not find the expression data file."
        ext = DATA_EXTS[0]
        filename = files[ext]
    MATRIX = arrayio.read(filename)

    # If no gene IDs were provided, then just make some up.
    if not MATRIX.row_names():
        header = "GENE.ID"
        MATRIX._row_order.append(header)
        x = ["R%s" % x for x in parselib.pretty_range(0, MATRIX.nrow())]
        MATRIX._row_names[header] = x
        synonyms = {}
        synonyms[arrayio.ROW_ID] = header
        #MATRIX = Matrix.add_synonyms(MATRIX, synonyms)
        MATRIX._synonyms.update(synonyms)
    if not MATRIX.col_names():
        header = arrayio.tdf.SAMPLE_NAME
        MATRIX._col_order.append(header)
        x = ["C%s" % x for x in parselib.pretty_range(0, MATRIX.ncol())]
        MATRIX._col_names[header] = x
        synonyms = {}
        synonyms[arrayio.COL_ID] = header
        #MATRIX = Matrix.add_synonyms(MATRIX, synonyms)
        MATRIX._synonyms.update(synonyms)
        

    # Read the clustering files.
    formats = [
        ("gtr", clusterio.read_gtr_file),
        ("atr", clusterio.read_atr_file),
        ("gtc", clusterio.read_gtc_file),
        ("atc", clusterio.read_atc_file),
        ("kgg", clusterio.read_kgg_file),
        ("kag", clusterio.read_kag_file),
        ]
    data = {}  # ext -> output
    for ext, read_fn in formats:
        if ext not in files:
            continue
        data[ext] = read_fn(files[ext])

    cluster_data = ClusterData(
        MATRIX,
        data.get("gtr"),
        data.get("atr"),
        data.get("gtc"),
        data.get("atc"),
        data.get("kgg"),
        data.get("kag"),
        )
    return cluster_data


def _write_cluster_files(MATRIX, SCALED, cluster_data, jobname):
    from arrayio import tab_delimited_format
    from genomicode import clusterio

    raise NotImplementedError

    matrix_file = "%s.cdt" % jobname
    tab_delimited_format.write(MATRIX, open(matrix_file, 'w'))
    scaled_file = "%s_s.cdt" % jobname
    tab_delimited_format.write(SCALED, open(scaled_file, 'w'))

    cd = cluster_data
    formats = [
        ("gtr", clusterio.write_gtr_file, cd.gene_tree),
        ("atr", clusterio.write_atr_file, cd.array_tree),
        ("gtc", clusterio.write_gtc_file, cd.gene_tree_cluster),
        ("atc", clusterio.write_atc_file, cd.array_tree_cluster),
        ("kgg", clusterio.write_kgg_file, cd.gene_cluster),
        ("kag", clusterio.write_kag_file, cd.array_cluster),
        ]

    for ext, write_fn, data in formats:
        if not data:
            continue
        outfile = "%s.%s" % (jobname, ext)
        write_fn(data, open(outfile, 'w'))


def _unlink_cluster_files(filestem):
    # Just remove all the files with the filestem.
    import os
    import filelib

    path, filestem = os.path.split(filestem)
    for file_ in os.listdir(path):
        if not file_.startswith(filestem):
            continue
        filename = os.path.join(path, file_)
        filelib.safe_unlink(filename)


def _guess_filestem(file_or_job):
    # Examples:
    # file_or_job           stem
    # test                  test
    # GSE5451.l2.mas5.gtr   GSE5451.l2.mas5
    # GSE1456.mas5.gz       GSE1456
    # GSE1456.mas5          GSE1456
    # out.dat               out
    # out.pcl               out
    # out.txt               out
    # /home/jchang/out.txt  /home/jchang/out
    # 
    # Rule:
    # - If the file doesn't exist, then use the whole thing as the
    #   stem.
    # - If there's a .gz, then chop it off.
    # - If there's an extension, then chop it off.
    import os
    
    if not _exists_nz(file_or_job):
        return file_or_job
    
    stem = file_or_job

    # Chop off the .gz at the end.
    COMPRESSION_EXTS = [".gz", ".bz2", ".zip"]
    s, e = os.path.splitext(stem)
    for ext in COMPRESSION_EXTS:
        if e.lower() == ext:
            stem = s
            break

    # Chop off one more extension, if it exists.
    stem, e = os.path.splitext(stem)

    return stem


def _convert_to_pcl(MATRIX):
    # Convert the matrix to PCL format.
    # Row names   <ID>  NAME
    # Col names   
    import arrayio
    from genomicode import Matrix

    # Select from the row names an ID and a NAME.
    id_name = _choose_gene_id(MATRIX)
    name_name = _choose_gene_label(MATRIX)

    # Make sure there aren't any blank gene IDs, or cluster will
    # complain.  Also, make sure they are unique.
    seen = {}
    for id_ in MATRIX.row_names(id_name):
        id_ = id_.strip()
        assert id_, "Missing gene IDs (header %s)." % id_name
        assert id_ not in seen, "Duplicate gene ID %s." % id_
        seen[id_] = 1

    # Should not use "GID" as column name for PCL file.  When
    # clustering, cluster will add another "GID" column, and then
    # there will be two columns called "GID".  Rename this to
    # something else, if necessary.
    pretty_id_name = id_name
    if pretty_id_name == "GID":
        pretty_id_name = "GID.OLD"
    if pretty_id_name == "NAME":
        # GCT files uses "NAME" for ID, which conflicts with PCL definition.
        pretty_id_name = "ID.NAME"
    pretty_name_name = "NAME"

    SAMPLE_NAME = arrayio.tab_delimited_format.SAMPLE_NAME 
    row_order = [pretty_id_name, pretty_name_name]
    col_order = [SAMPLE_NAME]
    row_names = {}
    col_names = {}
    synonyms = {}
    
    row_names[pretty_id_name] = MATRIX.row_names(id_name)
    row_names[pretty_name_name] = MATRIX.row_names(name_name)
    col_names[SAMPLE_NAME] = MATRIX.col_names(arrayio.COL_ID)
    synonyms[arrayio.ROW_ID] = pretty_id_name
    synonyms[arrayio.COL_ID] = SAMPLE_NAME

    pcl_matrix = Matrix.InMemoryMatrix(
        MATRIX.slice(), row_names=row_names, col_names=col_names,
        row_order=row_order, col_order=col_order, synonyms=synonyms)
    #pcl_matrix = Matrix.add_synonyms(x, synonyms)
    assert arrayio.pcl_format.is_matrix(pcl_matrix)
    return pcl_matrix


def _exists_nz(filename):
    import os
    import stat
    
    if not os.path.exists(filename):
        return None
    if os.stat(filename)[stat.ST_SIZE] > 0:
        return filename
    return None
