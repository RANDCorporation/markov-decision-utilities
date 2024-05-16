module MarkovDecisionUtilities

# load
using CSV
using Dates
using DataFrames
using Statistics
# load local modules
using DataStructures
using MarkovModels
using MDUSetups
using SupportFunctions

# export functions
export get_average_outcomes_by_edge
export run_all_markov_scenarios
export run_package_simulation!



"""
Get average outcomes by edge (each transition probability pair i,j) in the
    Markov Chain.

# Constructs

```
get_average_outcomes_by_edge(
    all_state_transitions::Set{Tuple{Int64,Int64}},
    vec_set_transitions::Array{Set{Tuple{Int64,Int64}},1},
    df_outcome_values::DataFrame
)
```

# Function Arguments

- `all_state_transitions`:
- `vec_set_transitions`:

"""
function get_average_outcomes_by_edge(
    all_state_transitions::Set{Tuple{Int64,Int64}},
    vec_set_transitions::Array{Set{Tuple{Int64,Int64}},1},
    df_outcome_values::DataFrame
)
    #
    # assume that df_outcome_values[:, :run_id] and vec_set_transitions are correctly ordered
    #
    # check fields
    field_run = :run_id
    field_outcome = :outcome_value
    check_fields!(df_outcome_values, [field_run, field_outcome])

    # check shapes
    all_runs = Set(df_outcome_values[:, field_run])
    if length(all_runs) != length(vec_set_transitions)
        error("error in get_average_outcomes_by_edge: the length of df_outcomes and vec_set_transitions should be the same.")
    elseif length(all_runs) != size(df_outcome_values)[1]
        error("error in get_average_outcomes_by_edge: multiple values of $(field_run) found in df_outcomes.")
    end

    # initialize output
    vec_i = Int64.(zeros(length(all_state_transitions)))
    vec_j = Int64.(zeros(length(all_state_transitions)))
    vec_mean_outcome = Float64.(zeros(length(all_state_transitions)))

    for ast in enumerate(all_state_transitions)
        k, tup = ast
        i, j = tup
        inds = findall(x -> (tup in x), vec_set_transitions)
        mean_outcome = mean(df_outcome_values[inds, field_outcome])
        vec_i[k] = i
        vec_j[k] = j
        vec_mean_outcome[k] = mean_outcome
    end

    return DataFrame(:i => vec_i, :j => vec_j, :expected_outcome_across_edge => vec_mean_outcome)
end



"""
Run all available scenarios contained in the input `MarkovDataset`

# Constructs

```
run_all_markov_scenarios(
        dataset::MarkovDataset,
        n_runs::Int64,
        n_time_periods::Int64,
        dict_starting_states::Dict{Int64, Int64} = dataset.attribute_scenario.field_maps["scenario_id_to_default_starting_state"]
    )
```

"""
function run_all_markov_scenarios(
        dataset::MarkovDataset,
        n_runs::Int64,
        n_time_periods::Int64,
        # define the scenarios to run and the associated starting state
        dict_starting_states::Dict{Int64, Int64} = dataset.attribute_scenario.field_maps["scenario_id_to_default_starting_state"]
    )

    # initialize some variables defined in the loop
    attr_state = 0
    df_abs_exp = []
    df_anconv_exp = []
    df_avbyedge_exp = []
    df_ent_exp = []
    df_exp = []
    df_outcomes_exp = []
    init_q = true

    # initialize some fields
    field_runkey = :run_id
    field_scenkey = dataset.field_scenario_key
    field_scenname = dataset.field_scenario_name
    field_statekey = dataset.field_state_key
    field_ov = dataset.get_table_field(-1, "outcome_value")

    # get scenarios to loop over
    scens_loop = [x for x in dataset.all_scenarios if x in keys(dict_starting_states)]


    if length(scens_loop) == 0
        vals = dataset.print_valid_scenarios()
        error("No valid scenarios were found in dict_starting_states. Valid scenarios are:\n$(vals)")
    end

    for enum in enumerate(scens_loop)

        i, scen = enum
        state_0 = dict_starting_states[scen]

        # get the state and transition associated with the scenario
        attr_state, transition_specification = dataset.get_markov_data(scen, true)
        # set the model and calculate entropies
        modmark = MarkovModel(transition_specification, attr_state, config)
        df_entropies = DataFrame(dataset.field_state_key => Int64.(collect(1:size(modmark.Q_default)[1])), :entropy => modmark.get_state_entropies())
        dict_entropies = build_dict(df_entropies[:, [dataset.field_state_key, :entropy]])

        # get the fundamental matrix and caluculate (1) expected # of visits, (2) expected time in each state, (3) total steps to absorbtion
        dict_canon = modmark.get_canonical()
        mat_fundamental = modmark.get_fundamental_matrix(dict_canon)
        vectimes_residence = get.(
            (modmark.attribute_states.field_maps["state_to_time"], ),
            dict_canon[:states_transient],
            1.0
        )
        mat_expected_time_in_state = mat_fundamental .* transpose(vectimes_residence)
        # convert to dataframes for output
        df_expected_time_in_state = DataFrame(mat_expected_time_in_state, Symbol.("expected_time_in_state_" .* string.(dict_canon[:states_transient])))
        df_fundamental = DataFrame(mat_fundamental, Symbol.("expected_number_of_visits_to_state_" .* string.(dict_canon[:states_transient])))
        df_expected_absorption = DataFrame(
            :expected_number_of_steps_to_absorption => mat_fundamental*ones(size(mat_fundamental)[1]),
            :expected_time_to_absorption => mat_expected_time_in_state*ones(size(mat_fundamental)[1])
        )
        # get probabiility of absorption in each state using
        mat_probability_absorption = modmark.get_probability_of_absorption_state(dict_canon)
        df_probability_absorption = DataFrame(mat_probability_absorption, Symbol.("probability_of_absorption_in_state_" .* string.(dict_canon[:states_absorption])))

        # build total analytical output
        df_analytical_convergence = hcat(
            DataFrame(:starting_state => dict_canon[:states_transient]),
            df_expected_absorption,
            df_expected_time_in_state[:, sort(names(df_expected_time_in_state); by = x -> tryparse(Int64, split(x, "_")[end]))],
            df_fundamental[:, sort(names(df_fundamental); by = x -> tryparse(Int64, split(x, "_")[end]))],
            df_probability_absorption[:, sort(names(df_probability_absorption); by = x -> tryparse(Int64, split(x, "_")[end]))]
        )
        sort!(df_analytical_convergence, [:starting_state])

        # run simulation and get output, state matrix, and transitions by run
        df_out, df_abs, vec_set_transitions = modmark.simulate_walks(n_runs, n_time_periods, state_0; return_type = :DataFrames, dict_entropy = dict_entropies);

        # some data manipulations
        df_out = leftjoin(df_out, df_entropies, on = :state => dataset.field_state_key)
        df_out[:, dataset.field_scenario_key] = Int64.(ones(nrow(df_out))*scen)
        df_abs[:, dataset.field_scenario_key] = Int64.(ones(nrow(df_abs))*scen)
        df_entropies[:, dataset.field_scenario_key] = Int64.(ones(nrow(df_entropies))*scen)
        df_analytical_convergence[:, dataset.field_scenario_key] = Int64.(ones(nrow(df_analytical_convergence))*scen)

        # build outcome values data frame
        df_outcomes = filter(x -> (x[:time_period] == n_time_periods), df_out)
        df_outcomes = leftjoin(df_outcomes, attr_state.table[:, [field_statekey, field_ov]], on = [field_statekey])
        df_outcomes = leftjoin(df_outcomes, df_abs, on = [field_runkey, field_scenkey])

        # add outcomes by edges
        df_average_by_edge = get_average_outcomes_by_edge(
            modmark.get_all_state_transitions(),
            vec_set_transitions,
            df_outcomes
        )
        nms_average_by_edge = names(df_average_by_edge)
        df_average_by_edge[:, dataset.field_scenario_key] = Int64.(scen * ones(size(df_average_by_edge)[1]))
        select!(df_average_by_edge, [[dataset.field_scenario_key]; Symbol.(nms_average_by_edge)])


        # output dataframes - allocate space if initializing, then just replace in future iterations
        if init_q
            df_abs_exp = [df_abs for x in scens_loop]
            df_anconv_exp = [df_analytical_convergence for x in scens_loop]
            df_avbyedge_exp = [df_average_by_edge for x in scens_loop]
            df_ent_exp = [df_entropies for x in scens_loop]
            df_exp = [df_out for x in scens_loop]
            df_outcomes_exp = [df_outcomes for x in scens_loop]
            init_q = false
        else
            df_abs_exp[i] = df_abs
            df_anconv_exp[i] = df_analytical_convergence
            df_avbyedge_exp[i] = df_average_by_edge
            df_ent_exp[i] = df_entropies
            df_exp[i] = df_out
            df_outcomes_exp[i] = df_outcomes
        end

    end

    # concatenate everything
    df_abs_exp = vcat(df_abs_exp...)
    df_avbyedge_exp = vcat(df_avbyedge_exp...)
    df_ent_exp = vcat(df_ent_exp...)
    df_exp = vcat(df_exp...)
    df_outcomes_exp = vcat(df_outcomes_exp...)

    select!(df_outcomes_exp, Not([:time_period, :entropy]))
    rename!(df_outcomes_exp, :state => Symbol("state_period_$(n_time_periods)"))

    # add a summary table
    df_outcomes_summary = combine(
        groupby(
            df_outcomes_exp[:, [field_scenkey, field_ov, :absorbed]], [field_scenkey]
        ),
        field_ov => mean,
        field_ov => var,
        :absorbed => sum
    )
    df_outcomes_summary = leftjoin(df_outcomes_summary, dataset.attribute_scenario.table[:, [field_scenkey, field_scenname]], on = field_scenkey)
    select!(df_outcomes_summary, Not([field_scenkey]))
    fields_ord = [[field_scenname]; Symbol.([x for x in sort(names(df_outcomes_summary)) if (Symbol(x) != field_scenname)])]
    df_outcomes_summary = df_outcomes_summary[:, fields_ord]

    return df_anconv_exp, df_avbyedge_exp, df_ent_exp, df_exp, df_outcomes_exp, df_outcomes_summary
end



##  run the simulation across specified scenarios
"""
Run the simulation across all scenarios defined in the `MarkovDataset` dataset.

# Constructs

```
run_package_simulation!()
```

```
run_package_simulation!(
    dataset::MarkovDataset,
    n_runs::Int64,
    n_time_periods::Int64,
    dir_output::String = dir_out,
    # define the scenarios to run and the associated starting state
    dict_starting_states::Dict{Int64, Int64} = dataset.attribute_scenario.field_maps["scenario_id_to_default_starting_state"]
)
```

# Function Arguments

- `dataset`: `MarkovDataset` to use to setup decision environment
- `n_runs`: number of walks to simulate
- `n_time_periods`: maximum number of time periods to run each walk for
- `dir_output`: output directory to write subdirectory to
- `dict_starting_states`: dictionary mapping each scenario to its starting state
"""
function run_package_simulation!(
        dataset::MarkovDataset,
        n_runs::Int64,
        n_time_periods::Int64,
        dir_output::String = dir_out,
        # define the scenarios to run and the associated starting state
        dict_starting_states::Dict{Int64, Int64} = dataset.attribute_scenario.field_maps["scenario_id_to_default_starting_state"]
    )

    # get the run time for the simulation
    analysis_run_time = Dates.format(now(), "yyyymmdd_HHMMSS")

    # get the directory of the dataset and copy over input
    dir_target = joinpath(dir_output, "$(analysis_run_time)_output_$(dataset.dataset_name)")
    @info("\nInitializting output using $(dataset.dir_dataset) as $(dir_target)...")
    cp(dataset.dir_dataset, dir_target)
    @info("\tDone.\n")

    # copy configuration file
    @info("\nCopying configuration file...")
    cp(config.fp_config, joinpath(dir_target, basename(config.fp_config)))
    @info("\tDone.")

    # run the simulation
    @info("\nStarting simulation of $(n_runs) runs and $(n_time_periods) time periods...\n\n")
    dfs_analytical_conversion_exp, df_expected_outcomes_by_edge, df_entropies_exp, df_sim, df_outcomes, df_outcomes_summary = run_all_markov_scenarios(
        dataset,
        n_runs,
        n_time_periods,
        dict_starting_states
    )
    @info("\nSimulation complete.\n")

    @info("\nExporting files...")
    # export simulation
    CSV.write(joinpath(dir_target, basename(fp_default_entropies)), df_entropies_exp)
    CSV.write(joinpath(dir_target, basename(fp_default_expected_outcome_by_edge)), df_expected_outcomes_by_edge)
    CSV.write(joinpath(dir_target, basename(fp_default_simulation_walks_out)), df_sim)
    CSV.write(joinpath(dir_target, basename(fp_default_outcomes_out)), df_outcomes)
    CSV.write(joinpath(dir_target, basename(fp_default_outcomes_summary_out)), df_outcomes_summary)
    # loop over analytical files--there should be one file for each scenario
    for tup in enumerate(dfs_analytical_conversion_exp)
        i, df_ac_out = tup
        scen_id = Int64(df_ac_out[1, dataset.field_scenario_key])
        CSV.write(
            joinpath(
                dir_target,
                basename(fpt_default_analytical_convergence(scen_id, dataset.field_scenario_key)
                )
            ),
            df_ac_out
        )
    end
    @info("\tDone.\n")

    try
        # set environment variable then compress
        @info("\nAttempting to ompress output...")
        ENV["COPYFILE_DISABLE"] = 1
        fn_tar = "$(basename(dir_target)).tar.gz"
        dir_targ = basename(dir_target)

        # change directory before compression
        dir_cur = pwd()
        cd(dirname(dir_target))
        run(`tar -cvzf $(fn_tar) $(dir_targ)`)
        cd(dir_cur)
        rm(dir_target, recursive = true)

    catch e
        @error ("Error attempting to create tarball: $(e)")
    end


    # notify of completion
    @info("Done.\n\n***\tSIMULATIONS COMPLETE\t***")

    return 0
end

# set the default behavior
function run_package_simulation!()
    dataset_def = MarkovDataset(joinpath(dir_datasets, config.dataset), config.default_transition_scenario);
    run_package_simulation!(dataset_def, config.default_num_runs, config.default_num_time_steps)
    return 0
end

# end the module
end
