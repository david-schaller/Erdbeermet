# -*- coding: utf-8 -*-

import re


def write_history(filename, history):
    
    with open(filename, 'w') as f:
        
        start = True
        
        for x, y, z, alpha, delta in history:
            delta_str = '[' + ','.join(str(d) for d in delta) + ']'
            if start:
                f.write(f"({x}, {y}: {z}) {alpha}; {delta_str}")
                start = False
            else:
                f.write(f"\n({x}, {y}: {z}) {alpha}; {delta_str}")
                

def _split_floats(floats):
    
    return [float(item) for item in floats.split(',')]


def parse_history(filename):
    
    event_regex = re.compile(r"\((\d+)\,\s*(\d+)\:\s*(\d+)\)\;?\s*(\d+\.?\d*e?-?\d+)\;\s*\[(?P<delta>(\s*\d+\.?\d*e?-?\d+,?)+)\]")
    
    with open(filename, 'r') as f:
        
        lines = f.readlines()
      
    history = []
    for line in lines:
        
        match = event_regex.match(line.strip())
        
        if match:
            
            x = int(match.group(1))
            y = int(match.group(2))
            z = int(match.group(3))
            alpha = float(match.group(4))
            delta = _split_floats(match.group('delta'))
            
        history.append((x, y, z, alpha, delta))
            
    return history


def _write_matrix(f, V, D):
    
    for i in range(len(V)):
        f.write(f'\n{V[i]}  ')
        for j in range(len(V)):
            f.write('{: 12.8f}'.format(D[i,j]))


def write_recognition(filename, tree, matrices=True):
    
    with open(filename, 'w') as f:
        
        start = True
        for v in tree.preorder():
            
            if not start:
                f.write('\n')
                f.write(80 * '-')
                f.write('\n')
            else:
                start = False
                
            f.write(f'n={v.n}\n')
            if v.R_step is not None:
                f.write('(result of R-step: ({},{}:{}){:.8f})\n'.format(*v.R_step))
            f.write(f'V={v.V}\n')
            f.write(f'total successes of this branch: {v.valid_ways}\n')
            
            if matrices and v.D is not None:
                f.write(f'Matrix on {v.n} elements:\n')
                _write_matrix(f, v.V, v.D)
                f.write('\n')
                
            if not v.valid_ways:
                f.write(f'reason of abort: {v.info}\n')
        
        