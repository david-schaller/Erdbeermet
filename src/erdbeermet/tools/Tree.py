# -*- coding: utf-8 -*-

from erdbeermet.visualize.RecognitionVis import Visualizer


__author__ = 'David Schaller'


class TreeNode:
    """Class 'TreeNode'.
    
    Components for class 'Tree'. Contains a list of children as well as a
    reference to the parent node.
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
        
    
    def visualize(self, decimal_prec=4, save_as=None, use_latex=False):
        
        Visualizer(self, decimal_prec=decimal_prec, save_as=save_as,
                   use_latex=use_latex)
    
    
    def _assert_integrity(self):
        
        for v in self.preorder():
            for child in v.children:
                if child is v:
                    raise RuntimeError('loop at {}'.format(v))
                if child.parent is not v:
                    raise RuntimeError('Tree invalid for '\
                                       '{} and {}'.format(v, child))
        
        return True