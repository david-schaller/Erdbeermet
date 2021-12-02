# -*- coding: utf-8 -*-


from erdbeermet.simulation import simulate
from erdbeermet.recognition import recognize


scenario = simulate(7)
scenario.print_history()
print(scenario.D)

print('-------------------- Recognition --------------------')
rec_tree = recognize(scenario.D)
print(rec_tree.to_newick())

rec_tree.visualize(decimal_prec=5, save_as=False)
