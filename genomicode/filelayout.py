"""

Describes the (real or imagined) layout of some files.

Classes:
FileObject
PathObject

Functions:
File   Generate a FileObject.
Path   Generate a PathObject.
is_file
is_path
get_paths
get_missing_names

"""

import os

class Node:
    def __init__(self, id):
        # id is a user-defined ID.
        self.id = id
        self.name = None
        self.parent = None
    def filename(self):
        assert self.name is not None
        name = self.name
        if self.parent is not None:
            name = os.path.join(self.parent.filename(), name)
        return name
##     def filename(self, id):
##         if id != self.id:
##             return None
##         name = self.name
##         if name is None:
##             name = "[%s]" % id
##         return name
    def __str__(self):
        x = self.name
        if x is None:
            x = "[%s]" % self.id
        return x
    def __getattr__(self, id):
        if self.id == id:
            return self.filename()
        raise AttributeError, id
##     def find_node(self, id):
##         if self.id == id:
##             return self
##         return None
##     def get_ids(self):
##         return [self.id]
        
class FileObject(Node):
    def __init__(self, id):
        Node.__init__(self, id)
    def __call__(self, name):
        # Should only call once when creating the file.
        assert name is not None
        assert self.name is None
        self.name = name
        return self

class PathObject(Node):
    def __init__(self, id):
        Node.__init__(self, id)
        self.subnodes = []
    def __del__(self):
        # Break the subnode-parent cycle, so can be garbage collected.
        for n in self.subnodes:
            assert n.parent == self
            n.parent = None
        self.subnodes = []
    def __call__(self, name, *subnodes):
        # Should only call once when creating this path.
        assert name is not None
        assert self.name is None
        self.name = name
        for n in subnodes:
            assert n.parent is None, "Node has multiple parents."
            n.parent = self
            self.subnodes.append(n)
        return self
##     def filename(self, id):
##         name = self.name
##         if name is None:
##             name = "[%s]" % self.id
##         if self.id == id:
##             return name
##         for node in self.subnodes:
##             x = node.filename(id)
##             if x is None:
##                 continue
##             return os.path.join(name, x)
##         return None
    def __str__(self):
        x = "\n".join([str(x) for x in self.subnodes])
        lines = x.split("\n")
        lines = ["  %s" % x for x in lines]
        x = "\n".join(lines)
        x = "%s/\n%s" % (Node.__str__(self), x)
        return x
    def __getattr__(self, id):
        if self.id == id:
            return self.filename()
        for node in self.subnodes:
            try:
                x = getattr(node, id)
            except AttributeError, x:
                continue
            return x
        raise AttributeError, id

##     def find_node(self, id):
##         if self.id == id:
##             return self
##         for node in self.subnodes:
##             x = node.find_node(id)
##             if x is not None:
##                 return x
##         return None
##     def find_node(self, id):
##         if self.id == id:
##             return self
##         for node in self.subnodes:
##             x = node.find_node(id)
##             if x is not None:
##                 return x
##         return None
##     def pop(self, id):
##         # Remove a node in place and return it.
##         if self.id == id:
##             return self
##         for i in range(len(self.subnodes)):
##             if self.subnodes[i].id == id:
##                 return self.subnodes.pop(i)
##     def push(self, parent_id, node):
##         # Return a boolean indicating whether node is successfully placed.
##         if self.id == parent_id:
##             self.subnodes.append(node)
##             return True
##         for node_ in self.subnodes:
##             if node_.push(parent_id, node):
##                 return True
##         return False
##     def move(self, id, parent_id):
##         # Return a boolean indicating whether node is successfully moved.
##         if not self.find_node(parent_id):
##             return False
##         node = self.pop(id)
##         if node is None:
##             return False
##         return self.push(parent_id, node)
##     def get_ids(self):
##         ids = [self.id]
##         for n in self.subnodes:
##             ids.extend(n.get_ids())
##         return ids
##     def __iter__(self):
##         yield self.id
##         for n in self.subnodes:
##             for id in n:
##                 yield id

class FileMaker:
    def __getattr__(self, id):
        return FileObject(id)
File = FileMaker()

class PathMaker:
    def __getattr__(self, id):
        return PathObject(id)
Path = PathMaker()

def walk(layout):
    # yield dirpath, dirnames, filenames
    stack = [layout]
    while stack:
        node = stack.pop(0)
        if isinstance(node, FileObject):
            continue
        # Push all the subnodes on the stack.
        stack.extend(node.subnodes)
        
        dirpath = node.filename()
        dirnames, filenames = [], []
        for n in node.subnodes:
            fn = n.filename()
            if isinstance(n, FileObject):
                if fn not in filenames:
                    filenames.append(fn)
            elif isinstance(n, PathObject):
                if fn not in dirnames:
                    dirnames.append(fn)
        yield dirpath, dirnames, filenames
    
