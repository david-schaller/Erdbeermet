# Erdbeermet

A Python library for generating, visualizing, manipulating and recognizing type R (pseudo)metrics
(German: **Er**zeugung, **D**arstellung, **Be**arbeitung und **E**rkennung von **R**-(Pseudo)**Met**riken).

## Installation

Download or clone the repo, go to the root folder of package and install it using the command:

    python setup.py install

#### Dependencies

The package requires Python 3.7 or higher.

* [Numpy](https://numpy.org)
* [Scipy](http://www.scipy.org/install.html)
* [Matplotlib](https://matplotlib.org/)

## Usage and description

### Simulation

#### R-steps / event histories

The elements of the `history` are the "R-steps" comprising a merge event of `x` and `y` that creates `z` with a parameter `alpha` (if `alpha` is 0 or 1, we have a pure branching event), and additionally a list of distance increments ("`delta`s") for the currently existing items.
A history written to file has one line per R-step of the form `(x, y: z) alpha; [deltas]`.
<details>
<summary>Example file: (Click to expand)</summary>

    (0, 0: 1) 1.0; [0.26,0.11]
    (1, 0: 2) 0.27; [0.004,0.18,0.2]
    (1, 2: 3) 0.61; [0.49,0.08,0.37,0.17]
    (0, 3: 4) 0.82; [0.21,0.06,0.03,0.42,0.02]
    (1, 4: 5) 0.81; [0.02,0.004,0.08,0.11,0.01,0.04]

</details>

#### The class `Scenario`

The module `erdbeermet.simulation` contains functions for simulating R matrices, writing simulated scenarios to file, and reloading them.
The class `Scenario` acts as a wrapper for histories of merge and branching events and the corresponding R matrix.
<details>
<summary>The class has the following attributes: (Click to expand)</summary>

| Attribute | Type | Description |
| --- | --- | --- |
| `N` | `int` | the number of items that were simulated |
| `history` | `list` of `tuple`s | the history of merge and branching events |
| `circular` | `bool` | indicates whether the scenario has a circular type R matrix |
| `D` | `N`x`N` `numpy` array | the distance matrix |

</details>

<details>
<summary>and the following functions: (Click to expand)</summary>

| Function | Parameter/return type | Description |
| --- | --- | --- |
| `distances()` | returns `N`x`N` `numpy` array | getter for the distance matrix |
| `get_history()` | returns `list` of `tuple`s | getter for the event history |
| `get_circular_ordering()` | returns `list` of `int`s | list representing the circular ordering (cut between item 0 and its predecessor); or `False` if the scenario is not circular |
| `write_history(filename)` | parameter of type `str` | write the event history into a file |
| `print_history()` |  | print the event history |

</details>

#### Generation of `Scenario`s

Instances of `Scenario` can be generated using the function `simulate()` in the module `erdbeermet.simulation`.

<details>
<summary>Parameters of this function (Click to expand)</summary>

| Parameter (with default values) | Type | Description |
| --- | --- | --- |
| `N` | `int` | number of items to be generated |
| `branching_prob=0.0` | `float` | probability that an event is a pure branching event; the default is 0.0, i.e., pure branching events are disabled |
| `circular=False` | `bool` | if set to True, the resulting distance matrix is guaranteed to be a circular type R matrix (only "neighbors" can be involves in merge events) |
| `clocklike=False` | `bool` | if set to True, the distance increment is equal for all items within each iteration (comprising a merge or branching event and the distance increments) and only varies between iteration; the default is False, in which case the increments are also drawn independently for the items within an iteration |

</details>

Simulated scenarios can be saved to a file (in form of their event history) using their function `write_history(filename)`.

    # simulate a scenario with six items
    scenario = simulate(6, branching_prob=0.3, circular=True, clocklike=False)

    # print the individual R-steps
    scenario.print_history()

    # write history to file
    scenario.write_history('path/to/history.txt')

    # reload history to file
    scenario_reloaded = load('path/to/history.txt')

Alternatively, the function `load(filename, stop_after=False)` returns an instance of `Scenario` after reading an event history from an earlier simulated scenario from a file.
The parameter `stop_after` can be set to an `int` x>0 to only include the R-steps until the x'th item is created, i.e., x-1 R-steps are executed.



### Recognition

(description)



## References

What are R pseudometrics/matrices?

* **Prohaska, S.J., Berkemer, S.J., Gärtner, F., Gatter, T., Retzlaff, N., The Students of the Graphs and Biological Networks Lab 2017, Höner zu Siederdissen, C., Stadler, P.F. (2017) Expansion of gene clusters, circular orders, and the shortest Hamiltonian path problem. Journal of Mathematical Biology. doi: 10.1007/s00285-017-1197-3.**
