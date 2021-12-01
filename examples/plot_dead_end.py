# -*- coding: utf-8 -*-

import numpy as np

from erdbeermet.simulation import Simulator, MetricFromEvents
from erdbeermet.recognition import recognize
import erdbeermet.tools.FileIO as FileIO
from erdbeermet.Box4 import Box4


# --- change filename here ---
history_file = './example_histories/eid0003_n6_history'


history = FileIO.parse_history(history_file)
sim = MetricFromEvents(history)
sim.print_history()
print(sim.D)

print('-------------------- Recognition --------------------')
rec_tree = recognize(sim.D, first_candidate_only=False)

print('---------')
print('# Successes:', rec_tree.successes)


for v in rec_tree.preorder():
    
    if v.n == 4 and v.info == 'spikes too short':
        V, D = v.V, v.D
        print(f'plotting box for {V} .........')
        print('matrix:')
        print(D)

        box = Box4(D, labels=V)
        print(box._diagonal_mode)
        print(box.solutions)
        box.plot()
