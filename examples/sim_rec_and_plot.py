# -*- coding: utf-8 -*-

import os


from erdbeermet.simulation import simulate
from erdbeermet.recognition import recognize



print('-------------------- Recognition --------------------\n')
scenario = simulate(7)
scenario.print_history()
print(scenario.D)



print('\n-------------------- Recognition --------------------\n')
recognition_tree = recognize(scenario.D, print_info=True)
print('\n')
print(recognition_tree.to_newick())
print('\n# successes:', recognition_tree.successes)


# write the recognition into a text file
result_dir = 'example_histories'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
recognition_tree.write_to_file(os.path.join(result_dir,
                                            'testfile_recognition.txt'))

# visulize the recognition
recognition_tree.visualize(decimal_prec=5, save_as=False)

# save the plot as pdf by specififying the parameter:
#      --> save_as='valid/path/filename_ending_with.pdf'

# recognition_tree.visualize(save_as='testfile_tree.svg')