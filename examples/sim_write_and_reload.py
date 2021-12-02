# -*- coding: utf-8 -*-


from erdbeermet.simulation import simulate, load
import os

result_dir = 'example_histories'
if not os.path.exists(result_dir):
    os.makedirs(result_dir)


# ------------------ non-linear R matrix ------------------

scenario = simulate(6)
print(scenario.D)
scenario.print_history()

# write history to file
scenario.write_history(os.path.join(result_dir, 'testfile_history'))

# reload history to file
scenario_reloaded = load(os.path.join(result_dir, 'testfile_history'))

print(scenario_reloaded.D)
scenario_reloaded.print_history()


# ------------------ linear R matrix ------------------

scenario_linear = simulate(6, linear=True)
print(scenario_linear.D)
scenario_linear.print_history()
print(scenario_linear.get_linear_ordering())

# write history to file
scenario.write_history(os.path.join(result_dir, 'testfile_history_linear'))

# reload history to file
scenario_reloaded2 = load(os.path.join(result_dir, 'testfile_history_linear'),
                          linear=True)

print(scenario_reloaded2.D)
print(scenario_reloaded2.get_linear_ordering())
