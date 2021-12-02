# -*- coding: utf-8 -*-


from erdbeermet.simulation import load
from erdbeermet.visualize.BoxGraphVis import plot_box_graph

# --- change filename here ---
history_file = './example_histories/eid0003_n6_history'

# load a history until the fourth item is created (i.e. three event)
scenario = load(history_file, stop_after=4)
scenario.print_history()
print(scenario.D)

# plot the 4x4 distance matrix as a box graph
plot_box_graph(scenario.D, labels=range(4))