# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx


def R_metric_on_4(p, q, a, dx=0, dy=0, dz=0, du=0):
    
    xy = p + q + dx + dy
    xz = (1 - a) * (p + q) + dx + dz
    xu = p + du + dx
    yz = a * (p + q) + dy + dz
    yu = q + du + dy
    zu = a * p + (1 - a) * q + dz + du
    
    return np.array([xy, xz, xu, yz, yu, zu])


class SimulatorInterface:
    
    
    def __init__(self, N, history, linear):
        
        self.N = N
        self.history = history
        self.linear = linear
        
    
    def distances(self):
        
        return self.D
    
    
    def get_history(self):
        
        return self.history
    
    
    def get_linear_ordering(self):
        
        return self.linear_ordering
    
    
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
    


class Simulator(SimulatorInterface):
    
    def __init__(self, N, branching_prob=0.0, linear=False, clocklike=False):
        
        self.N = N
        
        if linear and branching_prob > 0.0:
            raise ValueError("Pure duplication events are not allowed for "\
                             "linear type R metrics!")
        
        self.branching_prob = branching_prob
        self.clocklike = clocklike
        
        history = []
        super().__init__(N, history, linear)
        
        self._run()
        

    def _run(self):
        
        D = np.zeros((self.N, self.N))
        
        if self.linear:
            L = nx.DiGraph()
            L.add_node(0)
        
        for z in range(1, self.N):
            
            # simple duplication event
            if z == 1 or np.random.random() < self.branching_prob:
                
                x = np.random.randint(z)
                y, alpha = x, 1.0           # for the history
                D[x, z] = 0.0
                D[z, x] = 0.0
                
                for u in range(z):
                    D[u, z] = D[x, z]
                    D[z, u] = D[x, z]
                    
                if self.linear:
                    L.add_edge(x, z)
                    
            # recombination event      
            else:
                if not self.linear:
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
            if not self.clocklike:
                delta = np.random.exponential(scale=1/self.N, size=z+1)
            else:
                delta = np.random.exponential(scale=1/self.N) * np.ones((z+1,))
            
            for p in range(z):
                for q in range(p+1, z+1):
                    
                    D[p, q] += delta[p] + delta[q]
                    D[q, p] = D[p, q]
                    
            self.history.append( (x, y, z, alpha, delta) )
        
        self.D = D
        
        if self.linear:
            self.linear_ordering = [u for u in nx.topological_sort(L)]
        else:
            self.linear_ordering = None
            

class MetricFromEvents(SimulatorInterface):
    
    
    def __init__(self, history, linear=False, stop_after=False):
        
        self.stop_after = stop_after
        
        if self.stop_after is False:
            N = len(history) + 1
        elif self.stop_after <= len(history) + 1:
            N = self.stop_after
        else:
            raise RuntimeError(f'not enough events to simulate {stop_after} items')
        super().__init__(N, history, linear)
        
        self._run()
    
    
    def _run(self):
        
        D = np.zeros((self.N, self.N))
        
        if self.linear:
            L = nx.DiGraph()
            L.add_node(0)
            
        for i in range(self.N-1):
            x, y, z, alpha, delta = self.history[i]
            
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
                    
                if self.linear:
                    L.add_edge(x, z)
                    
            # recombination event      
            else:
                if self.linear:
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
                raise RuntimeError(f"Invalid length of delta array for z={z}!")
                        
            for p in range(z):
                for q in range(p+1, z+1):
                    
                    D[p, q] += delta[p] + delta[q]
                    D[q, p] = D[p, q]
        
        self.D = D
        
        if self.linear:
            self.linear_ordering = [u for u in nx.topological_sort(L)]
        else:
            self.linear_ordering = None


if __name__ == '__main__':
    
    p = 0.01
    q = 1.0
    a = 0.3
    dx = 0.5
    dy = 0.5
    dz = 0.5
    du = 0.5
    
#    from erdbeermet.Box4 import Box4
    
#    b = R_metric_on_4(p, q, a, dx=dx, dy=dy, dz=dz, du=du)
#    print(b)
#    box = Box4(b)
#    print(box._diagonal_mode, box.solutions)
#    box.plot()
    
    sim = Simulator(6)
    print(sim.D)
    sim.print_history()
    
    siml = Simulator(6, linear=True)
    print(siml.D)
    siml.print_history()
    print(siml.get_linear_ordering())
    
    siml2 = MetricFromEvents(siml.history, linear=True)
    print(siml2.D)
    print(siml2.get_linear_ordering())
    
    import os
    import FileIO
    
    result_dir = 'testfiles'
    
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    
    FileIO.write_history(os.path.join(result_dir, 'history'), sim.history)