"""

Functions:
is_valid_gene_id
select_valid_gene_id
select_valid_gene_id_I

Classes:
AbstractGeneSet          Abstract base class.
NoHeaderFileGeneSet      No column headers, just specify a column.
HeaderFileGeneSet        GeneSet from a file with column headers.
RandomGeneSet            Draw random genes from a Matrix.
GMTGeneSet               GMT format file.

GeneSetHomologConverter  Change geneset to another organism.

"""
import os, sys

class AbstractGeneSet:
    def __init__(self, gene_id_name):
        self.gene_id_name = gene_id_name

    def __iter__(self):
        for gene in self.get_genes():
            yield gene

    # Implement in a derived class.
    def get_genes(self):
        # Return the genes in the geneset.
        # Should only return valid gene IDs.
        raise NotImplementedError
    
    # Optional.  I provide a default implementation.  Can implement an
    # optimized version if desired.
    def __len__(self):
        return len(self.get_genes())
    def __getitem__(self, index):
        return self.get_genes()[index]
    def get_indexes(self, dataset):
        # Return the indexes of the genes.
        # Although geneset implements the iterator protocol, it's
        # faster to call get_genes than to iterate over a list.
        x = dataset._index(row=self.get_genes(), row_header=self.gene_id_name)
        I_row, I_col = x
        assert I_col is None
        return I_row

class GeneSet(AbstractGeneSet):
    def __init__(self, gene_id_name, genes):
        AbstractGeneSet.__init__(self, gene_id_name)
        self.genes = genes[:]
    def get_genes(self):
        return self.genes
    
class NoHeaderFileGeneSet(AbstractGeneSet):
    def __init__(self, gene_id_name, filename, column_num):
        AbstractGeneSet.__init__(self, gene_id_name)
        # Save the parameters and load the geneset when necessary.
        self.filename = filename
        self.column_num = column_num
        self.genes = None
    def _get_genes(self):
        import iolib
        from filelib import openfh

        i = self.column_num
        data = openfh(self.filename).read()
        x = [cols[i] for cols in iolib.split_tdf(data)]
        x = iolib.strip_each(x)
        x = {}.fromkeys(x).keys()
        x = select_valid_gene_id(x)
        return x
    def get_genes(self):
        if self.genes is None:
            self.genes = self._get_genes()
        return self.genes

class HeaderFileGeneSet(NoHeaderFileGeneSet):
    def __init__(self, gene_id_name, filename, column_header):
        from filelib import openfh
        
        handle = openfh(filename)
        header = handle.readline().rstrip("\r\n").split("\t")
        assert column_header in header, (
            "Invalid column name: %s" % column_header)
        colnum = header.index(column_header)
        NoHeaderFileGeneSet.__init__(self, gene_id_name, handle, colnum)

class RandomGeneSet(AbstractGeneSet):
    def __init__(self, gene_id_name, dataset, num_probes):
        # num_probes is the total number of probes on the chip that
        # should be matched with each geneset.
        assert num_probes > 0, "Need to match at least 1 probe"
        AbstractGeneSet.__init__(self, gene_id_name)
        self.dataset = dataset
        self.num_probes = num_probes
        self.fixed = None

        # Index the IDs in the data set for faster selection.
        ids = dataset.row_names(gene_id_name)
        I_valid = select_valid_gene_id_I(ids)
        assert self.num_probes <= len(I_valid), \
               "Not enough gene IDs for sampling [%d/%d]." % (
            len(I_valid), self.num_probes)

        id2indexes = {}      # id -> list of indexes
        all_singles = True
        for i in I_valid:
            id = ids[i]
            if id not in id2indexes:
                id2indexes[id] = [i]
            else:
                id2indexes[id] += [i]
                all_singles = False
                
        # id, list of indexes, number of indexes
        population = [None] * len(id2indexes)
        if all_singles:
            for i, (id, indexes) in enumerate(id2indexes.iteritems()):
                population[i] = id, indexes[0], 1
        else:
            for i, (id, indexes) in enumerate(id2indexes.iteritems()):
                population[i] = id, indexes, len(indexes)
        self.gene_ids = ids
        self.population = population
        self.all_singles = all_singles
        
    def _make_gene_set(self, num_probes):
        # Make a random geneset.  Return a list of the indexes into
        # the original data set.
        if num_probes is None:
            num_probes = self.num_probes
        
        while 1:
            # Choose a random set of genes.  Since genes may have many
            # probe sets, this may result in too many probe sets.
            # Using random.shuffle on ids is very slow.  random.sample
            # generates many calls to random and set.add.
            selected = sample(self.population, self.num_probes)
            
            # If every gene points to exactly 1 probe, then the number
            # of genes is the same as the number of probes.
            if self.all_singles:
                num_indexes = len(selected)
            else:
                num_indexes = sum([x[2] for x in selected])

            # If there are too many genes, then back off.
            while num_indexes > self.num_probes:
                x = selected.pop(0)
                num_indexes -= x[2]
                
            # If I have the right number of indexes, then stop.
            # Otherwise, I have too many, so try again.
            if num_indexes == self.num_probes:
                break
            # Gene set fail!  Try again.  Does not happen very often.

        # Extract the indexes from the selected genes.
        if self.all_singles:
            indexes = [x[1] for x in selected]
        else:
            indexes = [None] * self.num_probes
            i = 0
            for x in selected:
                for y in x[1]:
                    indexes[i] = y
                    i += 1
        return indexes
    def fix(self):
        # Fix so that get_genes always returns the same gene set.
        self.fixed = self._make_gene_set()
    def get_indexes(self, dataset, num_probes=None):
        # Return a list of indexes.
        assert self.dataset is dataset
        num_probes = num_probes or self.num_probes
        if self.fixed is not None:
            indexes = self.fixed
        else:
            indexes = self._make_gene_set(num_probes)
        assert len(indexes) == num_probes
        return indexes
    def get_genes(self, num_probes=None):
        # Return a list of random genes.
        indexes = self.get_indexes(self.dataset, num_probes=num_probes)
        geneset = [self.gene_ids[i] for i in indexes]
        return geneset

class GMTGeneSet(AbstractGeneSet):
    def __init__(self, gene_id_name, filename, geneset_name, *more_genesets):
        AbstractGeneSet.__init__(self, gene_id_name)
        # Save the parameters and load the geneset when necessary.
        self.filename = filename
        self.geneset_names = [geneset_name] + list(more_genesets)
        self.genes = None
    def _get_genes(self):
        import iolib
        from filelib import openfh

        data = openfh(self.filename).read()
        # <name> <comment> <gene1> <gene2> ... <genen>
        genes = []
        for cols in iolib.split_tdf(data):
            if cols[0] not in self.geneset_names:
                continue
            genes.extend(cols[2:])
        if not genes:
            raise AssertionError, "I could not find gene set: %s" % \
                  ",".join(self.geneset_names)
        x = genes
        x = iolib.strip_each(x)
        x = {}.fromkeys(x).keys()
        x = select_valid_gene_id(x)
        return x
    def get_genes(self):
        if self.genes is None:
            self.genes = self._get_genes()
        return self.genes

class GeneSetHomologConverter(AbstractGeneSet):
    def __init__(self, geneset, converter):
        # converter is a dict of old_gene_id -> new_gene_id.
        self.geneset = geneset
        self.converter = converter
        AbstractGeneSet.__init__(self, geneset.gene_id_name)
    def get_genes(self):
        x = self.geneset.get_genes()
        x = [self.converter.get(x) for x in x]
        x = {}.fromkeys(x).keys()
        x = select_valid_gene_id(x)
        return x
    
def is_valid_gene_id(id):
    if id is None:
        return 0
    id = id.strip()
    if id in ["", "0", "---"]:
        return 0
    if id.find("///") >= 0:
        return 0
    return 1

def select_valid_gene_id(L):
    I = select_valid_gene_id_I(L)
    return [L[i] for i in I]

def select_valid_gene_id_I(L):
    return [i for i, x in enumerate(L) if is_valid_gene_id(x)]

def sample(L, n):
    # Choose n random objects from L.
    import random
    return random.sample(L, n)

try:
    #raise ImportError
    import cGeneSet
except ImportError:
    pass
else:
    this_module = sys.modules[__name__]
    for name in cGeneSet.__dict__.keys():
        if name.startswith("__"):
            continue
        this_module.__dict__[name] = cGeneSet.__dict__[name]
