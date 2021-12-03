# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx

import erdbeermet.tools.FileIO as FileIO


class Scenario:
    
    def __init__(self, N, history, circular):
        
        self.N = N
        self.history = history
        self.circular = circular
        
        self._build_matrix()
        
    
    def distances(self):
        
        return self.D
    
    
    def get_history(self):
        
        return self.history
    
    
    def get_circular_ordering(self):
        
        if not self.circular:
            return False
        
        visited = {0}
        ordering = [0]
        while True:
            succ = next(self.L.successors(ordering[-1]))
            if succ in visited:
                break
            else:
                ordering.append(succ)
                visited.add(succ)
        
        return ordering
    
    
    def write_history(self, filename):
        
        FileIO.write_history(filename, self.history)
    
    
    def print_history(self):
        
        for x, y, z, alpha, delta in self.history:
            
            delta_str = '[' + ','.join(str(d) for d in delta) + ']'
            print(f"({x}, {y}: {z}) {alpha}; {delta_str}")
    
    
    def _build_matrix(self):
        
        self.D = np.zeros((self.N, self.N))
        D = self.D
    
        if self.circular:
            self.L = nx.DiGraph()
            L = self.L
            L.add_node(0)
            L.add_edge(0, 0)
            
        for x, y, z, alpha, delta in self.history:
            
            # simple duplication event
            if (x == y or
                (x is None) or (y is None) or
                alpha == 1.0 or alpha == 0.0):
                
                if x is None or alpha == 0.0:
                    x = y
                    
                D[x, z] = 0.0
                D[z, x] = 0.0
                
                for u in range(z):
                    D[u, z] = D[x, z]
                    D[z, u] = D[x, z]
                    
                if self.circular:
                    L.remove_edge(x, y)
                    L.add_edge(x, z)
                    L.add_edge(z, y)
                    
            # recombination event      
            else:
                if self.circular:
                    if not L.has_edge(x, y):
                        raise RuntimeError(f'{x} and {y} are not neighbors')
                        
                    L.remove_edge(x, y)
                    L.add_edge(x, z)
                    L.add_edge(z, y)
                
                for u in range(z):
                    if u != x and u != y:
                        
                        d = alpha * D[x, u] + (1 - alpha) * D[y, u]
                        D[u, z] = d
                        D[z, u] = d
                        
                D[z, x] = (1 - alpha) * D[x, y]
                D[x, z] = (1 - alpha) * D[x, y]
                D[z, y] = alpha * D[x, y]
                D[y, z] = alpha * D[x, y]
            
            # distance increment, i.e., independent evolution after event
            if len(delta) != z + 1:
                raise RuntimeError(f'invalid length of delta array for z={z}')
                        
            for p in range(z):
                for q in range(p+1, z+1):
                    D[p, q] += delta[p] + delta[q]
                    D[q, p] = D[p, q]


def random_history(N, branching_prob=0.0, circular=False, clocklike=False):
    
    history = []
    
    if circular:
        L = nx.DiGraph()
        L.add_node(0)
        L.add_edge(0, 0)
    
    for z in range(1, N):
        
        # simple duplication event
        if z == 1 or np.random.random() < branching_prob:
            
            x = np.random.randint(z)
            y, alpha = x, 1.0           # for the history
                
            if circular:
                y = next(L.successors(x))
                L.remove_edge(x, y)
                L.add_edge(x, z)
                L.add_edge(z, y)
                
        # recombination event      
        else:
            if not circular:
                x, y = np.random.choice(z, size=2, replace=False)
            else:
                x = np.random.randint(z)
                y = next(L.successors(x))
                    
                L.remove_edge(x, y)
                L.add_edge(x, z)
                L.add_edge(z, y)
                
            alpha = np.random.random()
        
        # distance increment, i.e., independent evolution after event
        if not clocklike:
            delta = np.random.exponential(scale=1/N, size=z+1)
        else:
            delta = np.random.exponential(scale=1/N) * np.ones((z+1,))
                
        history.append( (x, y, z, alpha, delta) )
    
    return history


def simulate(N, branching_prob=0.0, circular=False, clocklike=False):
    
    if circular and branching_prob > 0.0:
        raise ValueError('pure duplication events are not allowed for '\
                         'circular type R metrics')
    
    return Scenario(N,
                    random_history(N, branching_prob=branching_prob,
                                   circular=circular,
                                   clocklike=clocklike),
                    circular)


def scenario_from_history(history, circular=False, stop_after=False):
        
    if stop_after is False:
        N = len(history) + 1
    elif stop_after <= len(history) + 1:
        N = stop_after
    else:
        raise RuntimeError(f'not enough events to simulate {stop_after} items')
    
    return Scenario(N, history[:N-1], circular)


def load(filename, circular=False, stop_after=False):
    
    return scenario_from_history(FileIO.parse_history(filename),
                                 circular=circular,
                                 stop_after=stop_after)


def R_metric_on_4(p, q, a, dx=0, dy=0, dz=0, du=0):
    
    xy = p + q + dx + dy
    xz = (1 - a) * (p + q) + dx + dz
    xu = p + du + dx
    yz = a * (p + q) + dy + dz
    yu = q + du + dy
    zu = a * p + (1 - a) * q + dz + du
    
    return np.array([xy, xz, xu, yz, yu, zu])

    