# -*- coding: utf-8 -*-

import numpy as np

from erdbeermet.simulation import Simulator, MetricFromEvents
from erdbeermet.recognition import recognize
import erdbeermet.FileIO as FileIO
from erdbeermet.Box4 import Box4

# --- change filename here ---
history_file = '../examples_ids/eid0003_n6_history'
history = FileIO.parse_history(history_file)

sim = MetricFromEvents(history, stop_after=2)
sim.print_history()
print(sim.D)

sim = MetricFromEvents(history, stop_after=3)
sim.print_history()
print(sim.D)

sim = MetricFromEvents(history, stop_after=4)
sim.print_history()
print(sim.D)

box = Box4(sim.D, labels=range(4))
print(box._diagonal_mode)
print(box.solutions)
print(box.first_solution())
box.plot()

sim = MetricFromEvents(history, stop_after=5)
sim.print_history()
print(sim.D)

sim = MetricFromEvents(history)
sim.print_history()
print(sim.D)