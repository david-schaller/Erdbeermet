# -*- coding: utf-8 -*-


__author__ = 'David Schaller'


from erdbeermet.simulation import Simulator, MetricFromEvents
from erdbeermet.recognition import recognize
import erdbeermet.tools.FileIO as FileIO


# --- change filename here ---
history_file = './example_histories/eid0003_n6_history'


history = FileIO.parse_history(history_file)
sim = MetricFromEvents(history)
sim.print_history()
print(sim.D)

print('-------------------- Recognition --------------------')
# it makes more sense to set 'first_candidate_only=False'
# even though True works (not the full tree)
rec_tree = recognize(sim.D, first_candidate_only=False)

print('---------')
print('# Successes:', rec_tree.successes)
      

# ------------------------------------------------------------------------------
   
# call visualize() function of the Tree instance 'rec_tree'
rec_tree.visualize()

# ------------------------------------------------------------------------------

# save the plot as pdf by specififying the parameter:
#      --> save_as='valid/path/filename_ending_with.pdf'

#rec_tree.visualize(save_as='testfile_tree.pdf')