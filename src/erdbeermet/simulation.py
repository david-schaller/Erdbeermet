# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx

import erdbeermet.tools.FileIO as FileIO


class Scenario:
    
    def __init__(self, N, linear):
        
        self.N = N
        self.linear = linear
        
        self.history = []
        self.D = np.zeros((self.N, self.N))
        
    
    def distances(self):
        
        return self.D
    
    
    def get_history(self):
        
        return self.history
    
    
    def get_linear_ordering(self):
        
        return self.linear_ordering
    
    
    def write_history(self, filename):
        
        FileIO.write_history(filename, self.history)
    
    
    def print_history(self):
        
        if hasattr(self, 'stop_after') and self.stop_after is not False:
            events = self.stop_after - 1
        else:
            events = len(self.history)
        
        for i in range(events):
            x, y, z, alpha, delta = self.history[i]
            
            delta_str = '[' + ','.join(str(d) for d in delta) + ']'
            print(f"({x}, {y}: {z}) {alpha}; {delta_str}")
        
        if hasattr(self, 'stop_after') and self.stop_after is not False:
            print(f'--- stopped after {self.stop_after} ---')


def simulate(N, branching_prob=0.0, linear=False, clocklike=False):
    
    if linear and branching_prob > 0.0:
        raise ValueError('pure duplication events are not allowed for '\
                         'linear type R metrics')
    
    scenario = Scenario(N, linear)
    scenario.branching_prob = branching_prob
    scenario.clocklike = clocklike
    
    D = scenario.D
    
    if linear:
        L = nx.DiGraph()
        L.add_node(0)
    
    for z in range(1, N):
        
        # simple duplication event
        if z == 1 or np.random.random() < branching_prob:
            
            x = np.random.randint(z)
            y, alpha = x, 1.0           # for the history
            D[x, z] = 0.0
            D[z, x] = 0.0
            
            for u in range(z):
                D[u, z] = D[x, z]
                D[z, u] = D[x, z]
                
            if linear:
                L.add_edge(x, z)
                
        # recombination event      
        else:
            if not linear:
                x, y = np.random.choice(z, size=2, replace=False)
            else:
                x = np.random.randint(z)
                
                if L.out_degree(x) > 0:
                    y = list(L.successors(x))[0]
                else:
                    y = x
                    x = list(L.predecessors(y))[0]
                    
                L.remove_edge(x, y)
                L.add_edge(x, z)
                L.add_edge(z, y)
                
            alpha = np.random.random()
            
            for u in range(z):
                if u != x and u != y:
                    
                    d = alpha * D[x, u] + (1 - alpha) * D[y, u]
                    D[u, z] = d
                    D[z, u] = d
                    
            D[z, x] = (1 - alpha) * D[x, y]
            D[x, z] = (1 - alpha) * D[x, y]
            D[z, y] = alpha * D[x, y]
            D[y, z] = alpha * D[x, y]
        
        # distance increment, i.e. independent evolution after event
        if not clocklike:
            delta = np.random.exponential(scale=1/N, size=z+1)
        else:
            delta = np.random.exponential(scale=1/N) * np.ones((z+1,))
        
        for p in range(z):
            for q in range(p+1, z+1):
                
                D[p, q] += delta[p] + delta[q]
                D[q, p] = D[p, q]
                
        scenario.history.append( (x, y, z, alpha, delta) )
    
    scenario.linear_ordering = [u for u in nx.topological_sort(L)] if linear \
                               else None
    
    return scenario
    

def scenario_from_history(history, linear=False, stop_after=False):
        
    if stop_after is False:
        N = len(history) + 1
    elif stop_after <= len(history) + 1:
        N = stop_after
    else:
        raise RuntimeError(f'not enough events to simulate {stop_after} items')
    
    scenario = Scenario(N, linear)
    
    D = scenario.D
    
    if linear:
        L = nx.DiGraph()
        L.add_node(0)
        
    for i in range(N-1):
        x, y, z, alpha, delta = history[i]
        
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
                
            if linear:
                L.add_edge(x, z)
                
        # recombination event      
        else:
            if linear:
                if y not in L.successors(x):
                    raise RuntimeError(f"'{x}' and '{y}' are not neighbors!")
                    
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
        
        # distance increment, i.e. independent evolution after event
        if len(delta) != z + 1:
            raise RuntimeError(f'invalid length of delta array for z={z}')
                    
        for p in range(z):
            for q in range(p+1, z+1):
                
                D[p, q] += delta[p] + delta[q]
                D[q, p] = D[p, q]
        
        scenario.history.append( (x, y, z, alpha, delta) )
                
    scenario.linear_ordering = [u for u in nx.topological_sort(L)] if linear \
                               else None
    
    return scenario


def load(filename, linear=False, stop_after=False):
    
    return scenario_from_history(FileIO.parse_history(filename),
                                 linear=linear,
                                 stop_after=stop_after)


def R_metric_on_4(p, q, a, dx=0, dy=0, dz=0, du=0):
    
    xy = p + q + dx + dy
    xz = (1 - a) * (p + q) + dx + dz
    xu = p + du + dx
    yz = a * (p + q) + dy + dz
    yu = q + du + dy
    zu = a * p + (1 - a) * q + dz + du
    
    return np.array([xy, xz, xu, yz, yu, zu])


if __name__ == '__main__':
    
    scenario = simulate(6)
    print(scenario.D)
    scenario.print_history()
    
    # write history to file
    import os
    result_dir = 'testfiles'
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    scenario.write_history(os.path.join(result_dir, 'history'))
    
    scenario_linear = simulate(6, linear=True)
    print(scenario_linear.D)
    scenario_linear.print_history()
    print(scenario_linear.get_linear_ordering())
    
    scenario_linear2 = scenario_from_history(scenario_linear.history,
                                             linear=True)
    print(scenario_linear2.D)
    print(scenario_linear2.get_linear_ordering())
    