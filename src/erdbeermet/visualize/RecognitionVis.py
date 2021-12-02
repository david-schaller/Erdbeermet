# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('font',**{'family':'serif','serif':['Palatino']})


__author__ = 'David Schaller'


class Visualizer:
    
    info_dict = {'no pseodometric': r'no pseodometric',
                 'negative delta/dxy': r'negative $\delta$/$d_{xy}$',
                 'no candidate': r'no candidate',
                 'spikes too short': r'spikes too short'}
    
    def __init__(self, tree, decimal_prec=4, save_as=None):
        
        self.tree = tree
        self.decimal_prec = decimal_prec
        
        self.edge_length = 0.5
        self.symbolsize = 0.03
        self.symbollw = 0.04
        self.leafs_per_vertical_unit = 10
        self.fs = 9
        self.symbol_zorder = 3       
        
        # print(tree.to_newick())
        
        self.distance_dict = {}
        self.colors = {}
        self.leaf_counter = 0
        self.node_positions = {}
        
        self.draw()
        
        if save_as:
            self.fig.savefig(save_as)
        
    
    def draw(self):
        
        self.fig, self.ax = plt.subplots()
        self.ax.set_aspect('equal')
        self.ax.invert_yaxis()
        
        self.initial_traversal()
        self.assign_positions()
        self.draw_edges()
        self.draw_nodes()
        
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.spines['bottom'].set_visible(False)
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['left'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        
        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        self.fig.set_size_inches(5*abs(xmax-xmin), 5*abs(ymax-ymin)+0.4)
        plt.tight_layout()
        plt.show()
        
    
    def initial_traversal(self):
        
        xmax = 0.0
        
        for v in self.tree.preorder():
            if not v.parent:
                self.distance_dict[v] = self.edge_length
            else:
                self.distance_dict[v] = (self.distance_dict[v.parent] +
                                         self.edge_length)
                if self.distance_dict[v] > xmax:
                    xmax = self.distance_dict[v]

            if not v.children:
                self.leaf_counter += 1
                
        self.ax.set_xlim(-0.1,xmax+0.5)
    
    
    def assign_positions(self):
        
        ymax = (self.leaf_counter-1)/self.leafs_per_vertical_unit
        self.ax.set_ylim(ymax+self.symbolsize*0.6, -self.symbolsize*0.6)
        
        yposition = 0
        for v in self.tree.postorder():
            if not v.children:
                self.node_positions[v] = (self.distance_dict[v],
                                          yposition)
                yposition += 1/self.leafs_per_vertical_unit
            else:
                ymean = (self.node_positions[v.children[0]][1] +
                         self.node_positions[v.children[-1]][1])/2
                self.node_positions[v] = (self.distance_dict[v],
                                          ymean)
    
    
    def draw_edges(self):
        
        for v in self.tree.preorder():
            if v.parent:
                self.ax.plot([self.node_positions[v.parent][0],
                              self.node_positions[v][0]],
                             [self.node_positions[v][1],
                              self.node_positions[v][1]],
                     color='black',
                     linestyle='-', linewidth=1)
            else:
                self.ax.plot([self.node_positions[v][0]-self.edge_length,
                              self.node_positions[v][0]],
                             [self.node_positions[v][1],
                              self.node_positions[v][1]],
                     color='black',
                     linestyle='-', linewidth=1)
            if v.children:
                self.ax.plot([self.node_positions[v][0],
                              self.node_positions[v][0]],
                             [self.node_positions[v.children[0]][1],
                              self.node_positions[v.children[-1]][1]],
                     color='black',
                     linestyle='-', linewidth=1)
    
    
    def draw_nodes(self):
        
        for v in self.tree.preorder():
            
            x, y = self.node_positions[v]
            
            color = 'red'
            if v.valid_ways:
                color = 'lightgreen'
            
            self.draw_circle(x, y, color=color)
            self.write_V_and_R_step(v)
            
            if not v.children:
                self.write_abort_info(v)
                
    
    def draw_circle(self, x, y, color='white'):
        
        fill = mpatches.Circle((x, y), self.symbolsize/2,
                               color=color, fill=True,
                               zorder=self.symbol_zorder)
        self.ax.add_patch(fill)
        outer = mpatches.Circle((x, y), self.symbolsize/2,
                                color='black', fill=False,
                                lw=self.symbolsize/self.symbollw,
                                zorder=self.symbol_zorder)
        self.ax.add_patch(outer)
    
    
    def write_V_and_R_step(self, v):
        
        x, y = self.node_positions[v]
        
        self.ax.text(x-self.symbolsize/2-0.005, y-0.002,
                     #f'n={v.n}, {v.V}',
                     f'{v.V}',
                     horizontalalignment='right',
                     verticalalignment='bottom',
                     fontsize=self.fs)
        
        if v.R_step:
            str_templ = '({},{}:{})  \u03B1={:.' + str(self.decimal_prec) + 'f}'
            
            self.ax.text(x-self.edge_length+self.symbolsize/2+0.005, y-0.002,
                         str_templ.format(*v.R_step),
                         horizontalalignment='left',
                         verticalalignment='bottom',
                         fontsize=self.fs)
    
    
    def write_abort_info(self, v):
        
        x, y = self.node_positions[v]
        
        text = v.info
        
        self.ax.text(x+self.symbolsize/2+0.02, y, text,
                     horizontalalignment='left',
                     verticalalignment='center',
                     fontsize=self.fs)
