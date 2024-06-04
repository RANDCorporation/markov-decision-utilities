# Markov Decision Utilities

Julia implementation for some Markov chain processes, including defining and implementing transition matrix from a sparse text file using variables, finding steady states and absorption probabilities, and simulating random walks.


# About This Package

This package was created as part of the research described in [Developing Models and Metrics to Assess the Impacts of Complexity in Operational Settings (RR-A1596-1)](https://www.rand.org/t/RRA1596-1), in which authors assessed mathematical strategies for quantifying complexity in wartime environments, and provide a framework that will help make complexity a concrete consideration in operational settings.

This project, _Complexity and Adversary Decisions: Developing Models and Metrics to Assess the Impacts of Complexity in Operational Settings_, was conducted between August 2021 and May 2022. The research was sponsored by Dr. Mark Linderman, Senior Scientist for Command and Control, Air Force Research Lab Information Directorate (AFRL RI) and conducted within the Force Modernization and Employment Program, RAND Project AIR FORCE. 

For more information about RAND Project AIR FORCE, visit [the rand.org website](https://www.rand.org/paf).

For questions about this package, contact [Matthew Sargent](mailto:msargent@rand.org)


#  0. Quick Start

It is easy to begin with the default settings located in the `complexity_modeling.config` file (using the `cuban_missile_crisis` example dataset). The do a basic run based on configuration defaults, open the Julia REPL, navigate to the `./julia` directory, and run the following lines of code:

```
# load the ./julia directory to the load path to find the modules
!(pwd() in LOAD_PATH) ? push!(LOAD_PATH, pwd()) : nothing

# load the MDU module
using MarkovDecisionUtilities

# run all scenarios
run_package_simulation!()
```

This will create an output dataset in the `./out` subdirectory. For more information on `run_package_simulation!`, see section 3 below. 

---


#  1. Entering and Interacting with Data: the `MarkovDataset` object

The MDU ingests data using an object called **`MarkovDataset`** (defined in `DataStructures.jl`). The `MarkovDataset` object requires two initialiation arguments, which are described below. Additionally, relevant data tables and associate fields (as well as field data type requirements) are detailed below.


##  `dir_dataset`

`dir_dataset` gives the directory containing the dataset to run. The base file path of the directory is assumed to be the name of the dataset in the MDU. In the MDU scripts, the directory `./ref/datasets` contains subdirectories with datasets. Each dataset must contain three files (described below):


###  Scenario Attribute Table (attribute_scenario.csv)

This file contains basic information about each scenario. It requires three fields:
    - `scenario_id` (type `Int64`): key field giving the scenario id, an integer. This index is used throughout the MDU to refer to scenarios.
    - `scenario_name` (type `String`): name of the scenario, used for tracking purposes.
    - `default_starting_state` (type `Int64`): a default starting state for the scenario to use in simulations.


###  State Attribute Table (attribute_state.csv)

`attribute_state.csv` gives information about **all possible states across all scenarios**. Different scenarios do not re-index states. They all operate on the same set of available decision states. Not every state must be used in every scenario; however, any states that are defined in a scenario must be included here, and states do not change across scenarios. This file has $2(s + 1)$ fields, where $s$ is the number of scenarios.

This file requires the following 2 fields (independent of number of scenarios) with entried defined for each state:
    - `state` (type `Int64`): the state number, an integer index that is used to define transition matrices in `transition_matrices_by_scenario.csv`
    - `name` (type `String`): name of the state, used for tracking (can be merged back into data for reporting purposes)

Additionally, `attribute_state.csv` requires 2 additional fields *for each scenario defined*:
    - `time_#` (type `Int64`): the number of times to spend in this state under scenario number `#`. The number should be entered as an integer; for example, for scenario 0, `time_0` gives the time spent in the state (1 means one step; 10 means ten steps), while `time_12` would give the times in each state for scenario 12.
    - `outcome_value_#` (type `Int64` or `Float64`): value used to score ending a simulation in this state. This is used to calcualte expected outcomes and simulated outcomes. For example, values < 0 may indicate undesirable outcomes in the decision space, while positive values may indicate more desirable outcomes. These values are set by the analyst and research team and depend on the context of the decision environment and analytical goal.


###  Transition Matrices by Scenario (transition_matrices_by_scenario.csv)

`transition_matrices_by_scenario.csv` defines the transition matrices for each scenario using a sparse array schema. Transitions are added as elements in a stochastic matrix (**NOTE:** the orientation of the matrix is set in `complexity_modeling.config`. This example uses the default, a row-stochastic matrix.); i.e., to specify a transition from $1 \to 2$, ensure that there is a row where $i = 1$ and $j = 2$. The transition probability associated with scenario $\#$ is entered in a row with the name `transition_expression_scenario_#`.

Users should keep the following in mind when entering transition matrices:
    - Transitions that aren't entered are assumed to occur with probability 0
    - A transition link must be entered as a row in this file if it exists in **any** scenario.
    - Links that do not exist in any scenario do not have to be entered in the file.
    - In row-stochastic orientation, all entries in a row must sum to 1 (i.e., for a given scenario, all entries with the same state index `i` must sum to 1)

`transition_matrices_by_scenario.csv` requires $2 + s$ fields. Two fields are required independent of the number of scenarios:
    - `i` (type `Int64`): the source or origin state (transition out)--must be defined in `attribute_state.csv`
    - `j` (type `Int64`): the target or end state (transition in)--must be defined in `attribute_state.csv`

Finally, the following field must be entered for each scenario defined in `attribute_scenario.csv`:
    - `transition_expression_scenario_#` (type `Float64` or `String`): Most users should only enter Float64 expressions. This field gives transition probabilities for scenario $\#$. However, it is possible to enter symbolic expressions, such as `"1 - A"` (where `A` is the symbolic variable) in this field, allowing users to evaluate a range of matrices. However, if doing so, an associated variable must be specified in the configuration file that gives a default value for `A` (for example `variable_A: 0.15`). The definition of a default symbolic variable value allows the model to run a number of scenarios without exploring over the symbolic variable. MarkovModel functions allow for users to set specific values of symbolic variables for individual runs.


##  `default_transition_scenario`

The `default_transition_scenario` gives the default scenario (integer) to run, when running a single simulation run, in the absence of a specified scenario. This scenario must be defined in `attribute_scenario.csv`.

---


#  2. Once a dataset is initialized, we can choose a scenario index and initialize a `MarkovModel` object

- the dataset .get_markov_data() function will validate scenarios and check if the specification is valid
- the `MarkovModel` structure performs:
    - checks on the transition matrix;
    - transformations of the transition matrix, such as canonical form;
    - analytical solutions (where possible); and
    - simulations (at the base level).


##  `MarkovModel` can be used to return the canonical form of the matrix (as well as the fundamental matrix)

The `explore_markov_model.ipynb` notebook contains evaluation cells that users can use to return a dictionary giving canonical form A as well as applicable submatrices (Q/R/I) + states; empty arguments returns the default transition matrix. See `?MarkovModel` for more information on `get_canonical`.


##  Next, we can use the `MarkovModel` to perform a random walk based on the starting state.

The default starting state is specified in the scenario attribute table and can be accessed as `dataset.attribute_scenario.field_maps["scenario_id_to_default_starting_state"]`.

Use `do_walk` to do a simple random walk in discrete time with two starting arguments: (`starting_state`, `n_time_steps`)

`do_walk` returns an ordered tuple with the following elements:

- discrete walk starting at state starting_state for n_time_stpes
- all state_changes
- the raw trials (if speciied -- this is an optional argument),
- the set of unique transitions
- the time to absorption  

See `?MarkovModel` for more information on `do_walk`


##  To perform an ensemble of simulated walks, use `MarkovModel.simulate_walks()`

`simulate_walks` generally takes three arguments: (`n_runs`, `n_time_steps`, `starting_state`)

In `explore_markov_model.ipynb`, users can run 10000 walks for 1000 periods using initial state `state_0` (defined above). The result of `MarkovModel.simulate_walks()` is a tuple of two data frames:
- `df_out`: a data frame long by `:run_id` (the index for simulation run) and `time_period` (the index for the number of time periods). The data frame includes two key output fields:
    - `state`: the state of the random walk at the given `time_period` for run `run_id`
    - `new_state`: a binary (0 or 1) indicating whether the state is a change from the last state (used for calculating metrics about state changes)

- `df_abs`: A summary data frame capturing information on absorption. For each `run_id`, it includes information on:
    - `absorption_time_period`: the time period that the walk reached an absorption state
    - `absorbed`: binary indicating whether or not the walk was absorbed
    - `n_state_changes`: the number of times the walk changed states.


####  Using the output data, information on entropies can be calculated.

See `?MarkovModel` for more information on `simulate_walks`.


##  Next, internal tools allow the comparison of some simulated observations to analytical solutions

- Simulated soluations are generated using `MarkovModel.simulate_walks()`; times to absorption and the probabilities of ending in an absorption state given a starting state can be calculated from simulations (convergence probabilities0
- Analytical solutions are also available

See commenting and code in the iJupyter notebook for a comparison of simulated convergence probabilities to those found in the analytical solution.

Note that the analytical solution is $B = FR$, where $F = (I - Q)^{-1}$ is the fundamental matrix and $R \in \mathbb{R}^{t \times a}$ is the upper-right submatrix of the canonical form, such that $R_{ij}$ gives the proabability of a transient state $i$ converging to absorption state $j$ (in this context, $i$ and $j$ are states in $A$, the canonical form of the transition matrix $Q$—not the original transition matrix $Q$).

---


#  3. Finally, the `MarkovDecisionUtilities` package provides a single command can be used to execute simulations of all scenarios included in the configuration dataset

`MarkovDecisionUtilities.run_package_simulation!()` will run all scenarios from the configuration dataset using defaults from the configuration. These include:

- Default number of runs: `config.default_num_runs`
- Default number of time steps: `config.default_num_time_steps`
- Default starting state: for each scenario, the default starting state comes from the scenario attribute table, which can be seen using `dataset.attribute_scenario.table`(applies to this Jupyter notebook, since `dataset` was initialized using the configuration data set specified in `config.dataset`)

- Running `run_package_simulation!()` will create a unique session key based off the date and time to track unique runs. Outputs are sent to a session key subdirectory created in `dir_out`. The directory contains analytical solutions, summary statistics, simulated walks, and all input data. These are stored in several files:

    - `run_package_simulation!()` writes `attribute_scenario.csv`, `attribute_state.cs`, `transition_matrices_by_scenario.csv`, and `complexity_modeling.config` to the session key subdirectory for archival purposes.
    - Analytical solutions by scenario (number `#`) are written to individual files called `analytical_solutions_$(dataset.field_scenario_key)_#.csv`. The files are long by starting state (field `starting_state`) and contain fields for the following info:
        1. the expectd time in each transient state $S$, given as `expected_time_in_state_$S$`;
        2. the expectd number of visits to state $S$, given as `expected_number_of_vists_to_state_$S$`; and
        3. the probability of absorption in state $S$, given as `probability_of_absorption_in_state_$S$`.
    - `expected_outcome_by_edge.csv`: The expected outcome (i.e., the simulation parameter outcome) across each edge (transitinon probability, $i \to j$) by scenario, calculated across all runs.
    - `outcome_values_summary.csv`: summary outcome values calculated for each scenario (but indexed by scenario name). This table is useful for a quick, easily-readable comparison of scenarios.
    - `outcome_values.csv`: selected outcomes for each `scenario_id` and `run_id`, including:
        - `state_period_TIMEPERIODMAX`: the state at the final time period (`TIMEPERIODMAX`);
        - `outcome_value`: the outcome value at final state or absorption;
        - `absorption_time_period`: the time period when the random walk is absorbed;
        - `absorbed`: whether or not the walk was absorbed;
        - `n_state_changes`: the number of state changes in the random walk;
        - `entropy_total_by_time`: total entropy over time;
        - `entropy_total_by_steps`: total entropy over all unique state changes;
        - `entropy_mean_by_time`: mean entropy over time;
        - `entropy_mean_by_steps`: mean entropy over all unique state changes;
    - `simulated_walks.csv`: random walks by scenario and run
    - `state_entropies.csv`: the input state entropies for each state and scenario.

For more information, see the docstring at `?run_package_simulation!` (this can be tried in the iJupyter notebook or from the Julia REPL).

## BibTeX Citation

If you re-use the Markov Decision Utilities code in a scientific publication, we would appreciate it if you use the following citations:

```
@book{Sargent2024,
    url       = {https://rand.org/t/RRA1596-1},
    year      = {2024},
    publisher = {RAND},
    pages     = {67},
    author    = {Matthew Sargent, Edward Parker, James Syme, and Vikram Kilambi},
    title     = {Developing Models and Metrics to Assess the Impacts of Complexity in Operational Settings}
}

@misc{Sargent2024,
    author    = {Matthew Sargent, Edward Parker, James Syme, and Vikram Kilambi},
    year      = {2024},
    publisher = {RAND},
    title     = {Markov Decision Utilities},
    url       = {https://github.com/RANDCorporation/markov-decision-utilities}
}
```


