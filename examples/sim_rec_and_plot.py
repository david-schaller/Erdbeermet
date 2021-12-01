# -*- coding: utf-8 -*-


from erdbeermet.simulation import Simulator
from erdbeermet.recognition import recognize


sim = Simulator(7)
sim.print_history()
print(sim.D)

print('-------------------- Recognition --------------------')
rec_tree = recognize(sim.D)
print(rec_tree.to_newick())

rec_tree.visualize(decimal_prec=5, save_as=False)
