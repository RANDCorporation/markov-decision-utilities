#include("support_functions.jl")
module DataStructures
####################
#    STRUCTURES    #
####################

using DataFrames

using SupportFunctions

export AttributeTable
export Configuration
export MarkovDataset

"""
Store information in an attribute table, which includes a key (or a columns of
    unique entries) for which attributes are defined in additional columns.
    Includes several methods and properties:

    * `get_ordered_attribute` (method for retrieving an attribute, or non-key
        column)
    * `key_values` (ordered unique entries in AttributeTable.key)
    * `field_maps` (dictionaries that map between the key and all other
        attribute fields)
    * `table` (the attribute table)


##  Constructs

struct AttributeTable(
    fp_table::Union{String,DataFrame}
    key::Symbol
    fields_to_dict::Union{Vector{Symbol}, Nothing}
)

##  Initialization Arguments (also accessible as properties)

- `fp_table`: file path to a CSV table or a DataFrame containing the attribute
    table to initialize
- `key`: field in `fp_table` used as the table key
- `fields_to_dict`: vector of fields to include in `field_maps`
    * Enter nothing to initialize **all** field maps (default)
    * Enter Vector{Symbol}() to initialize **no** field maps

"""
struct AttributeTable

    fp_table::Union{String,DataFrame}
    key::Symbol
    fields_to_dict::Union{Vector{Symbol}, Nothing}
    field_maps
    get_ordered_attribute
    key_values
    table

    function AttributeTable(
        fp_table::Union{String,DataFrame},
        key::Symbol,
        fields_to_dict::Union{Vector{Symbol}, Nothing} = nothing
    )

        # verify table exists and check keys
        if isa(fp_table, String)
            table = read_csv(check_path(fp_table, false), true)
        elseif isa(fp_table, DataFrame)
            table = fp_table
            fp_table = "NONE"
        end

        fields_all_nonkey = [x for x in Symbol.(names(table)) if x != key]
        if isa(fields_to_dict, Nothing)
            fields_to_dict = fields_all_nonkey
        else
            fields_to_dict = [x for x in fields_to_dict if x != key]
        end
        check_fields!(table, [[key]; fields_to_dict])

        # check key
        if length(Set(table[:, key])) < nrow(table)
            error("Invalid key '$(key)' found in $(fp_table): the key is not unique. Check the table and specify a unique key.")
        end
        # next, create dict maps
        field_maps = Dict()
        for fld in fields_to_dict
            field_fwd = "$(key)_to_$(fld)"
            field_rev = "$(fld)_to_$(key)"

            field_maps[field_fwd] = build_dict(table[:, [key, fld]])
            # check for 1:1 correspondence before adding reverse
            vals_unique = Set(table[:, fld])
            if (length(vals_unique) == nrow(table))
                field_maps[field_rev] = build_dict(table[:, [fld, key]])
            end
        end

        # get key values
        key_values = sort(unique(table[:, key]))

        # order an attribute by the key values
        function get_ordered_attribute(attribute::String)
            fm = "$(key)_to_$(attribute)"
            check_keys!(field_maps, [fm])
            return [field_maps[fm][x] for x in key_values]
        end


        return new(fp_table, key, fields_to_dict, field_maps, get_ordered_attribute, key_values, table)
    end
end


"""
Ingest configuration parameters for Markov Decision Utilities. The configuration
    file includes the following properties:

    * dataset
        The name of dataset, which is specified as a subfolder in
            `./ref/datasets`, to run
    * default_num_runs
        The default number of runs to execute when conducting a random walk
            simulation. Note that relevant simlulation functions can generally
            be run with a number of simluation runs specified as a function
            argument.
    * default_num_time_steps
        The default number of time steps to run when conducting a random walk
            simulation. Note that relevant simlulation functions can generally
            be run with a number time steps specified as a function argument.
    * default_transition_scenario
        The default scenario index (integer) to call if unspecified in function
            arguments. Transition scenarios are 0-indexed, and it is generally
            recommended to treat 0 as a "baseline".
    * dict_parameters_variable
        A dictionary that contains default *variable* values. The MDU allows
            users to asses symbolic transition matrices as part of scenarios,
            where those values are entered as string variables (e.g., as
            `A` or `1 - A`) in a transition matrix.
    * get
        Function to retrieve configuration parameters using a string or symbol
        (e.g., `Configuration.get(:default_num_time_steps)`)
    * parameters_required
        Vector of parameters that are required to run the Markov Decision
        Utility
    * parameters_variable
        Vector symbolic variables implied by the Configuration entries
    * transition_correction_threshold
        Threshold used to correct row sums in row stochastic transition matrices
        (or column sums in column stochastic transition matrices); accounts for
        floating point errors that may be present in entries. If sum errors
        exceed this threshold, then the simulation will stop.
    * transition_stochastic_orientation
        Stochastic orientation of input transition matrices. Entries are "row"
        or "column" (default is "row").


##  Constructs

struct Configuration(
    fp_config::String
)

##  Initialization Arguments (also accessible as properties)

- `fp_config`: file path to the configuration file used to initialize MDU

"""
##  CONFIGURATION file
struct Configuration

    fp_config::String
    dataset
    default_num_runs
    default_num_time_steps
    default_transition_scenario
    dict_parameters_variable
    get
    parameters_required
    parameters_variable
    transition_correction_threshold
    transition_stochastic_orientation

    function Configuration(fp_config::String)

        dict_conf = parse_config(fp_config)
        # match string to define variables
        matchstr_var = "variable_"
        # set required and variable configuration parametres
        parameters_required = ["dataset", "default_num_runs", "default_num_time_steps", "default_transition_scenario", "transition_correction_threshold", "transition_stochasitc_orientation"]
        parameters_variable = [replace(x, matchstr_var => "") for x in keys(dict_conf) if !(x in parameters_required) & occursin(matchstr_var, x)]
        # clean the dictionary
        keys_clean = copy(collect(keys(dict_conf)))
        for k in keys_clean
            k_new = replace(k, matchstr_var => "")
            if !(k_new in parameters_variable) & !(k_new in parameters_required)
                pop!(dict_conf, k)
            elseif k_new in parameters_variable
                val = dict_conf[k]
                dict_conf[k_new] = val
                pop!(dict_conf, k)
            end
        end

        # check the keys for required
        check_keys!(dict_conf, parameters_required, "Required configuration parameters ##MISSINGHERE## not found. Check the configuration file at '$(fp_config)'.")

        # retrive configuration values
        function get(key::String)
            if key in keys(dict_conf)
                return dict_conf[key]
            else
                valid_opts = print_valid_values(string.(collect(keys(dict_conf))))
                error("Invalid key $(key). Valid op
                tions are $(valid_opts).")
            end
        end
        function get(key::Symbol)
            return get(string(key))
        end

        ##  ASSIGN REQUIRED CONFIGURATION PARAMETERS AND CHECK

        dataset = string(get("dataset"))

        default_num_runs = Int64(get("default_num_runs"))
        # check values
        if default_num_runs < 0
            error("Invalid configuration parameter value of default_num_runs = $(default_num_runs): default_num_runs should be > 0.")
        end

        default_num_time_steps = Int64(get("default_num_time_steps"))
        # check values
        if default_num_time_steps < 0
            error("Invalid configuration parameter value of default_num_time_steps = $(default_num_time_steps): default_num_time_steps should be > 0.")
        end

        default_transition_scenario = Int64(get("default_transition_scenario"))
        # check values
        if default_transition_scenario < 0
            error("Invalid configuration parameter value of default_transition_scenario = $(default_transition_scenario): default_transition_scenario should be > 0.")
        end

        transition_correction_threshold = Float64(get("transition_correction_threshold"))
        # check values
        if transition_correction_threshold < 0
            error("Invalid configuration parameter value of transition_correction_threshold = $(transition_correction_threshold): transition_correction_threshold should be >= 0.")
        end

        transition_stochastic_orientation = Symbol(get("transition_stochasitc_orientation"))
        # check values
        if !(transition_stochastic_orientation in [:row, :column])
            error("Invalid configuration parameter transition_stochasitc_orientation=$(transition_stochasitc_orientation): transition_stochasitc_orientation should be 'row' or 'column'.")
        end


        # variable parameters
        dict_parameters_variable = Dict(zip(parameters_variable, get.(parameters_variable)))

        return new(
            fp_config,
            dataset,
            default_num_runs,
            default_num_time_steps,
            default_transition_scenario,
            dict_parameters_variable,
            get,
            parameters_required,
            parameters_variable,
            transition_correction_threshold,
            transition_stochastic_orientation
        )
    end
end



##  MARKOV DATASET

"""
# Information

The `MarkovDataset` structure is a key component of the Markov Data Utilities
    package collection. It ingests data with the following structre:

    datasets (also referred to herein as decision environments)
        |_ simulation properties and scenario definition (`attribute_scenario.csv`)
        |_ transition matrices by scenario (`transition_matrices_by_scenario.csv`)
        |_ universe of states (`attribute_state.csv`)

#  Constructs

```
struct MarkovDataset(
    dir_dataset::String,
    default_transition_scenario::Int64
)
```

# Initialization Arguments

- `dir_dataset`: Directory representing the dataset to initialize. The directory
    must contain the files described above.

- `default_transition_scenario`: The default scenario (integer) to run, when
    running a single simulation run, in the absence of a specified scenario
"""
struct MarkovDataset

    dir_dataset::String
    default_transition_scenario::Int64
    all_scenarios
    attribute_scenario
    attribute_state
    dataset_name
    field_i
    field_j
    field_scenario_default_starting_state
    field_scenario_key
    field_scenario_name
    field_state_key
    field_state_name
    fn_attribute_scenario
    fn_attribute_state
    fn_transition_by_scenario
    get_markov_data
    get_table_field
    print_valid_scenarios
    required_files


    function MarkovDataset(dir_dataset::String, default_transition_scenario::Int64)

        ##  DIRECTORY AND FILE CHECKS AND OPERATIONS

        # check existence of the directory
        dir_dataset = check_path(dir_dataset, false)
        dataset_name = basename(dir_dataset)
        # specify some default filenames and required files
        fn_attribute_scenario = "attribute_scenario.csv"
        fn_attribute_state = "attribute_state.csv"
        fn_transition_by_scenario = "transition_matrices_by_scenario.csv"
        required_files = [fn_attribute_scenario, fn_attribute_state, fn_transition_by_scenario]
        # check for existence
        for fn in required_files
            if !ispath(joinpath(dir_dataset, fn))
                error("Path to required file '$(fn)' not found in dataset '$(dataset_name)'")
            end
        end


        ##  FIELD SPECIFICATION AND FUNCTIONS

        # function for getting fields based on scenario
        function get_table_field(scen::Int64, type_in::String)
            dict_types_to_field_prepend = Dict{String, String}(
                "time" => "time",
                "outcome_value" => "outcome_value",
                "transition_expression" => "transition_expression_scenario"
            )
            check_keys!(dict_types_to_field_prepend, [type_in])

            prep = dict_types_to_field_prepend[type_in]
            out = (scen < 0) ? prep : "$(prep)_$(scen)"

            return Symbol(out)
        end
        # dispatch to a default
        function get_table_field(type_in::String)
            return get_table_field(default_transition_scenario, type_in)
        end


        # set some required fields
        field_i = :i
        field_j = :j
        field_scenario_key = :scenario_id
        field_scenario_name = :scenario_name
        field_scenario_default_starting_state = :default_starting_state
        field_state_key = :state
        field_state_name = :name

        # set scenario attribute and check for required fields
        attribute_scenario = AttributeTable(joinpath(dir_dataset, fn_attribute_scenario), field_scenario_key)
        check_fields!(attribute_scenario.table, [field_scenario_key, field_scenario_name, field_scenario_default_starting_state])

        attribute_state = AttributeTable(joinpath(dir_dataset, fn_attribute_state), field_state_key)
        transitions_by_scenario = read_csv(joinpath(dir_dataset, fn_transition_by_scenario), true)
        check_fields!(transitions_by_scenario, [field_i, field_j])

        # set of all available scenarios
        all_scenarios = attribute_scenario.key_values
        if !(default_transition_scenario in all_scenarios)
            error("default_transition_scenario $(default_transition_scenario) not found in the scenario attribute table.")
        end

        # set of all available states
        all_states = attribute_state.key_values

        function print_valid_scenarios()
            dict_map = copy(attribute_scenario.field_maps["scenario_id_to_scenario_name"])
            dict_map = sort(dict_map)
            dict_repl = Dict{String, String}("(" => "\tscenario ", "," => ":", ")" => "")
            str_out = join(zip(keys(dict_map), values(dict_map)), "\n")
            for repl in keys(dict_repl)
                str_out = replace(str_out, repl => dict_repl[repl])
            end

            return str_out
        end
        ##  SET SOME FUNCTIONS

        # get datasets for the MarkovModel
        function get_markov_data(
            scen::Int64,
            clean_missing_states_q::Bool = true
        )

            if !(scen in all_scenarios)
                error("Scenario $(scen) not found in the dataset. Valid scenarios are\n$(print_valid_scenarios())")
            end

            # extract the appropriate subset
            fields_ext_state = get_table_field.(scen, ["time", "outcome_value"])
            fields_rnm_state = get_table_field.(-1, ["time", "outcome_value"])
            # transition fields
            field_ext_trans = get_table_field(scen, "transition_expression")
            field_rnm_trans = :expression#get_table_field(-1, "transition_expression")

            # get attribute table for the stat
            att_state = copy(attribute_state.table[:, [[field_state_key, field_state_name]; fields_ext_state]])
            rename!(att_state, Dict(zip(fields_ext_state, fields_rnm_state)))
            # clean data types
            att_state[!, get_table_field(-1, "time")] = Int64.(coalesce(att_state[!, get_table_field(-1, "time")], 1))
            att_state[!, get_table_field(-1, "outcome_value")] = Float64.(coalesce(att_state[!, get_table_field(-1, "outcome_value")], 0))
            # build the attribute table
            att_state = AttributeTable(att_state, field_state_key)

            # get the transition matrix and drop 0s
            table_trans = transitions_by_scenario[:, [field_i, field_j, field_ext_trans]]
            rename!(table_trans, field_ext_trans => field_rnm_trans)
            table_trans = filter(x -> (x[field_rnm_trans] > 0), table_trans)

            # add dummy rows?
            if clean_missing_states_q
                states_to_add = setdiff(all_states, union(table_trans[:, field_i], table_trans[:, field_j]))
                df_append = DataFrame(field_i => states_to_add, field_j => states_to_add)
                df_append[:, field_rnm_trans] = ones(length(states_to_add))
                table_trans = vcat(table_trans, df_append)
            end

            return att_state, table_trans
        end


        return new(
            dir_dataset,
            default_transition_scenario,
            all_scenarios,
            attribute_scenario,
            attribute_state,
            dataset_name,
            field_i,
            field_j,
            field_scenario_default_starting_state,
            field_scenario_key,
            field_scenario_name,
            field_state_key,
            field_state_name,
            fn_attribute_scenario,
            fn_attribute_state,
            fn_transition_by_scenario,
            get_markov_data,
            get_table_field,
            print_valid_scenarios,
            required_files
        )
    end
end

# end the module
end
