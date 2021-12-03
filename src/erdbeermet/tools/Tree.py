# -*- coding: utf-8 -*-

from erdbeermet.visualize.RecognitionVis import Visualizer
from erdbeermet.tools.FileIO import write_recognition


__author__ = 'David Schaller'


class TreeNode:
    """Tree node class for type R matrix recognition tree.
    
    Components for class 'Tree'. Contains a list of children as well as a
    reference to the parent node. 
    The root of the tree corresponds to the full list of items and the full
    distance matrix, and every other node corresponds the remaining distance
    matrix after the application of a specific R-step (x, y: z)alpha.
    
    Attributes
    ----------
    n : int
        Number of items remaining after removal of z.
    V : list
        List of items after removal of z.
    D : 2-dimensional numpy array
        Distance matrix after removal of z and update of the distances.
    R_step : tuple
        x, y, z, and alpha representing the R-step (x, y: z)alpha that was
        applied last.
    valid_ways : int
        Total number of recognition paths leading to a success in the subtree
        below this node.
    info : str
        Info string why the recognition failed after the application of
        this R-step (if this is the case).
    """
    
    def __init__(self, n, V, D=None, R_step=None):
        
        self.parent = None
        self.children = []
        
        self.n = n
        self.V = V
        self.D = D
        self.R_step = R_step
        
        self.valid_ways = 0
        self.info = ''
        
        
    def __str__(self):
        
        token = f'<<n={self.n}'
        
        if self.R_step:
            token += '|({},{}:{}){:.4f}'.format(*self.R_step)
        
        return token + '>>'
        
        
    def add_child(self, child):
        
        self.children.append(child)
        child.parent = self
        

class Tree:
    """Tree for type R matrix recognition.
    
    The root of the tree corresponds to the full list of items and the full
    distance matrix, and every other node corresponds the remaining distance
    matrix after the application of a specific R-step (x, y: z)alpha.
    
    Attributes
    ----------
    root : TreeNode
        The root corresponds to the full distance matrix.
    """
    
    def __init__(self, root):
        
        self.root = root
    
    
    def preorder(self):
        """Generator for preorder traversal of the tree."""
        
        def _preorder(node):
            yield node
            for child in node.children:
                yield from _preorder(child)
        
        if self.root:
            yield from _preorder(self.root)
        else:
            yield from []
    
    
    def postorder(self):
        """Generator for postorder traversal of the tree."""
        
        def _postorder(node):
            for child in node.children:
                yield from _postorder(child)
            yield node
        
        if self.root:
            yield from _postorder(self.root)
        else:
            yield from []
            
    
    def inner_vertices(self):
        """Generator for inner vertices in preorder."""
        
        def _inner_vertices(node):
            if node.children:
                yield node
                for child in node.children:
                    yield from _inner_vertices(child)
        
        if self.root:
            yield from _inner_vertices(self.root)
        else:
            yield from []
            
    
    def edges(self):
        """Generator for all edges of the tree."""
        
        def _edges(node):
            for child in node.children:
                yield (node, child)
                yield from _edges(child)
        
        if self.root:
            yield from _edges(self.root)
        else:
            yield from []
            
    
    def inner_edges(self):
        """Generator for all inner edges of the tree."""
        
        def _inner_edges(node):
            for child in node.children:
                if child.children:
                    yield (node, child)
                    yield from _inner_edges(child)
        
        if self.root:
            yield from _inner_edges(self.root)
        else:
            yield from []


    def to_newick(self, node=None):
        """Recursive Tree --> Newick (str) function."""
        
        def _to_newick(node):
            if not node.children:
                return str(node)
            else:
                s = ''
                for child in node.children:
                    s += _to_newick(child) + ','
                return "({}){}".format(s[:-1], node)
        
        if self.root:
            return _to_newick(self.root) + ';'
        else:
            return ';'
        
    
    def visualize(self, decimal_prec=4, save_as=None):
        
        Visualizer(self, decimal_prec=decimal_prec, save_as=save_as)
    
    
    def write_to_file(self, filename):
        
        write_recognition(filename, self)
    
    
    def _assert_integrity(self):
        
        for v in self.preorder():
            for child in v.children:
                if child is v:
                    raise RuntimeError('loop at {}'.format(v))
                if child.parent is not v:
                    raise RuntimeError('Tree invalid for '\
                                       '{} and {}'.format(v, child))
        
        return True