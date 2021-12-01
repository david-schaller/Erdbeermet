# -*- coding: utf-8 -*-

import numpy as np

from erdbeermet.simulation import Simulator, MetricFromEvents
from erdbeermet.recognition import recognize
import erdbeermet.FileIO as FileIO
from erdbeermet.Box4 import Box4

# --- change filename here ---
history_file = './example_HISTORIES/eid0003_n6_history'
history = FileIO.parse_history(history_file)

# load a history until the fourth item is created (i.e. three event)
sim = MetricFromEvents(history, stop_after=4)
sim.print_history()
print(sim.D)

# plot the 4x4 distance matrix as a box graph
box = Box4(sim.D, labels=range(4))
box.plot()