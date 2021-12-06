# -*- coding: utf-8 -*-

import numpy as np

import erdbeermet.tools.FileIO as FileIO


class Scenario:
    """Scenario for the generation of a type R matrix.
    
    Attributes
    ----------
    N : int
        Number of items.
    history : list of tuples
        The history of merge and branching events.
    circular : bool
        Indicates whether the scenario has a circular type R matrix
    D : 2-dimensional numpy array
        The distance matrix.
    
    See Also
    --------
    simulate()
    scenario_from_history()
    load()
    """
    
    def __init__(self, history):
        """Constructor for Scenario class.
        
        Parameters
        ----------
        history : list of tuples
            The history of merge and branching events.
        """
        
        self.N = len(history) + 1
        self.history = history
        
        self._build_matrix()
        
    
    def distances(self):
        """Distance matrix of the scenario.
        
        Returns
        -------
        2-dimensional numpy array
        """
        
        return self.D
    
    
    def get_history(self):
        """History of merge and branching events.
        
        Returns
        -------
        list of tuples
        """
        
        return self.history
    
    
    def get_circular_order(self):
        """Circular order of the items.
        
        Returns
        -------
        list or bool
            A list representing the circular order (cut between item 0 and
            its predecessor); or False if the scenario is not circular.
        """
        
        if not self.circular:
            return False
        
        visited = {0}
        order = [0]
        while True:
            succ = self._circ_order[order[-1]]
            if succ in visited:
                break
            else:
                order.append(succ)
                visited.add(succ)
        
        return order
    
    
    def write_history(self, filename):
        """Write the event history into a file.
        
        Parameters
        ----------
        filename : str
            Path and filename.
        """
        
        FileIO.write_history(filename, self.history)
    
    
    def print_history(self):
        """Print the event history."""
        
        for x, y, z, alpha, delta in self.history:
            
            delta_str = '[' + ','.join(str(d) for d in delta) + ']'
            print(f"({x}, {y}: {z}) {alpha}; {delta_str}")
    
    
    def _build_matrix(self):
        """Generate the distance matrix and determine whether it is circular.
        """
        
        self.D = np.zeros((self.N, self.N))
        D = self.D
        
        # initialize circular as True and set to False if non-neighbor merge
        # event is encountered
        self.circular = True
        self._circ_order = {0: 0}
            
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
                    old_succ = self._circ_order[x]
                    self._circ_order[x] = z
                    self._circ_order[z] = old_succ
                    
            # recombination event      
            else:
                if self.circular:
                    if self._circ_order[x] == y:
                        self._circ_order[x] = z
                        self._circ_order[z] = y
                    elif self._circ_order[y] == x:
                        self._circ_order[y] = z
                        self._circ_order[z] = x
                    else:
                        self.circular = False
                
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
    """Generate a random history of merge and branching events.
    
    Parameters
    ----------
    N : int
        Number of items.
    branching_prob : float, optional
        Probability that an event is a pure branching event. The default is
        0.0, i.e., pure branching events are disabled.
    circular : bool, optional
        If set to True, the resulting history is guaranteed to produce a
        circular type R matrix. The default is False.
    clocklike : bool, optional
        If set to True, the distance increment is equal for all items within
        each iteration (comprising a merge or branching event and the distance
        increments) and only varies between iteration. The default is False,
        in which case the increments are also drawn independently for the
        items within an iteration.
        
    Returns
    -------
    list of tuples
        Each tuple corresponds to an iteration comprising a merge or branching
        event and the distance increments.
    """
    
    history = []
    
    if circular:
        successors = {0: 0}
    
    for z in range(1, N):
        
        # simple duplication event
        if z == 1 or np.random.random() < branching_prob:
            
            x = np.random.randint(z)
            y, alpha = x, 1.0           # for the history
                
            if circular:
                y = successors[x]
                successors[x] = z
                successors[z] = y
                
        # recombination event      
        else:
            if not circular:
                x, y = np.random.choice(z, size=2, replace=False)
            else:
                x = np.random.randint(z)
                y = successors[x]
                successors[x] = z
                successors[z] = y
                
            alpha = np.random.random()
        
        # distance increment, i.e., independent evolution after event
        if not clocklike:
            delta = np.random.exponential(scale=1/N, size=z+1)
        else:
            delta = np.random.exponential(scale=1/N) * np.ones((z+1,))
                
        history.append( (x, y, z, alpha, delta) )
    
    return history


def simulate(N, branching_prob=0.0, circular=False, clocklike=False):
    """Simulate a random type R matrix.
    
    Parameters
    ----------
    N : int
        Number of items.
    branching_prob : float, optional
        Probability that an event is a pure branching event. The default is
        0.0, i.e., pure branching events are disabled.
    circular : bool, optional
        If set to True, the resulting distance matrix is guaranteed to be a
        circular type R matrix. The default is False.
    clocklike : bool, optional
        If set to True, the distance increment is equal for all items within
        each iteration (comprising a merge or branching event and the distance
        increments) and only varies between iteration. The default is False,
        in which case the increments are also drawn independently for the
        items within an iteration.
        
    Returns
    -------
    Scenario
        Comprises the history of merge and branching events as well as the
        distance matrix.
    """
    
    # if circular and branching_prob > 0.0:
    #     raise ValueError('pure duplication events are not allowed for '\
    #                      'circular type R metrics')
    
    return Scenario(random_history(N, branching_prob=branching_prob,
                                   circular=circular,
                                   clocklike=clocklike))


def scenario_from_history(history, stop_after=False):
    """Generate a type R matrix from a list of merge and branching events.
    
    Parameters
    ----------
    history : list of tuples
        The history of merge and branching events.
        
    Returns
    -------
    Scenario
        Comprises the history of merge and branching events as well as the
        distance matrix.
    """
        
    if stop_after is False:
        N = len(history) + 1
    elif stop_after <= len(history) + 1:
        N = stop_after
    else:
        raise RuntimeError(f'not enough events to simulate {stop_after} items')
    
    return Scenario(history[:N-1])


def load(filename, stop_after=False):
    """Generate an event history from a file and generate the type R matrix.
    
    Parameters
    ----------
    filename : str
        Path and filename.
        
    Returns
    -------
    Scenario
        Comprises the history of merge and branching events as well as the
        distance matrix.
    """
    
    return scenario_from_history(FileIO.parse_history(filename),
                                 stop_after=stop_after)


def R_metric_on_4(p, q, a, dx=0, dy=0, dz=0, du=0):
    
    xy = p + q + dx + dy
    xz = (1 - a) * (p + q) + dx + dz
    xu = p + du + dx
    yz = a * (p + q) + dy + dz
    yu = q + du + dy
    zu = a * p + (1 - a) * q + dz + du
    
    return np.array([xy, xz, xu, yz, yu, zu])

    