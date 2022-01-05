# -*- coding: utf-8 -*-

from erdbeermet.simulation import simulate
from erdbeermet.recognition import recognize


counter = 1

while True:
    
    print(counter)
    
    scenario = simulate(6)
    
    recognition_tree = recognize(scenario.D, print_info=False,
                                 first_candidate_only=True)
    
    if recognition_tree.valid_ways == 0:
        print('\n')
        scenario.print_history()
        
        print('\n')
        print(scenario.D)
    
        print('\n')
        print(recognition_tree.to_newick())
        print('\n# successes:', recognition_tree.successes)
    
        # visulize the recognition
        recognition_tree.visualize(decimal_prec=5, save_as=False)
        
        break
    
    counter += 1
