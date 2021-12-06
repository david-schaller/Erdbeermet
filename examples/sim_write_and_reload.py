# -*- coding: utf-8 -*-


from erdbeermet.simulation import simulate, load
import os

result_dir = 'example_histories'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)


# ------------------ non-circular R matrix ------------------

scenario = simulate(6)
print(scenario.D)
scenario.print_history()

# write history to file
scenario.write_history(os.path.join(result_dir, 'testfile_history'))

# reload history to file
scenario_reloaded = load(os.path.join(result_dir, 'testfile_history'))

print(scenario_reloaded.D)
scenario_reloaded.print_history()


# ------------------ circular R matrix ------------------

scenario_circular = simulate(6, circular=True)
print(scenario_circular.D)
scenario_circular.print_history()
print(scenario_circular.get_circular_order())

# write history to file
scenario_circular.write_history(os.path.join(result_dir, 'testfile_history_circular'))

# reload history to file
scenario_reloaded2 = load(os.path.join(result_dir, 'testfile_history_circular'))

print(scenario_reloaded2.D)
print(scenario_reloaded2.get_circular_order())
