# -*- coding: utf-8 -*-

from erdbeermet.simulation import load
from erdbeermet.recognition import recognize
from erdbeermet.visualize.BoxGraphVis import plot_box_graph


# --- change filename here ---
history_file = './example_histories/eid0003_n6_history'


scenario = load(history_file)
scenario.print_history()
print(scenario.D)

print('-------------------- Recognition --------------------')
rec_tree = recognize(scenario.D, first_candidate_only=False)

print('---------')
print('# Successes:', rec_tree.successes)


# visualize recognition
rec_tree.visualize()

# visualize metric on remaining 4 elements
for v in rec_tree.preorder():
    
    if v.n == 4 and v.info == 'spikes too short':
        V, D = v.V, v.D
        print(f'plotting box for {V} .........')
        print('matrix:')
        print(D)

        plot_box_graph(D, labels=V)
        break
