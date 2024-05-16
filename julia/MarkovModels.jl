
"""
Implement Markov transition matrices and associated simulation methods.


"""
module MarkovModels

using DataFrames
using LinearAlgebra
using Random
using SparseArrays

using DataStructures
using SupportFunctions

export MarkovModel
"""
# Information

The `MarkovModel` structure includes a range of methods for interacting with
    transition matrices and random walks. The `MarkovModel` structure generally
    includes the following capabilities:

    - checks on the transition matrix;
    - transformations of the transition matrix, such as canonical form;
    - analytical solutions (where possible); and
    - simulations (at the base level).

# Constructs

```
MarkovModel(
    model_specification::DataFrame
    attribute_states::AttributeTable
    configuration::Configuration
)
```


# Initialization Arguments

See `DataStructures.jl` for information on AttributeTable and Configuration
    structures.

- `model_specification`: DataFrame used to input transition matrices. The
    DataFrame should include the following fields:
    - `i`
    - `j`
    - `expression`
- `attribute_states`: AttributeTable
- `configuration`: Configuration object


# Methods

The following methods and properties are defined in `MarkovModel`

- `check_transition_matrix!`
- `do_walk`
- `find_absorption`
- `find_next_state`
- `get_absorbing_states`
- `get_all_state_transitions`
- `get_canonical`
- `get_fundamental_matrix`
- `get_probability_of_absorption_state`
- `get_state_attribute`
- `get_state_entropies`
- `Q_default` (property: default transition matrix defined from the input
    DataFrame `model_specification`)
- `simulate_absorptions`
- `simulate_walks`
- `swap_inds!`
- `swap_states!`
- `transition_matrix_from_data_frame`

Information on these methods is provided below.


## `check_states!`

Check the validity of input states against specified attribute_states.

###  Constructs

```
check_states!(
	states::Vector{Int64}
)
```


## `check_transition_matrix!`

Check the transition matrix `mat` with stochastic orientation `stochastic_orientation`
    and return a valid matrix with stochastic orientation `return_orientation`.
    Use `correction_thresh` to adjust invalid probability transition vectors.

###  Constructs

```
check_transition_matrix!(
    mat::SparseMatrixCSC{Float64, Int64},
    correction_thresh::Float64 = configuration.transition_correction_thresh,
    stochastic_orientation::Symbol = configuration.transition_stochastic_orientation,
    return_orientation::Symbol = :row
)
```


## `format_probs_as_columns`

###  Constructs
```
format_probs_as_columns(
	vec_ext_inds::Array{Tuple{Int64, Int64}, 1}
)
```


## `get_absorbing_states`

Get the absorbing states associated with Sparse matrix `mat`

###  Constructs

```
get_absorbing_states(
    mat::SparseMatrixCSC{Float64,Int64}
)
```


## `get_all_state_transitions`

Get the set of all state transitions (``i \\to j``) contained in `mat`

###  Constructs

```
get_all_state_transitions(
	mat::SparseMatrixCSC{Float64,Int64}
)
```

- use default matrix Q_default
```
get_all_state_transitions()
```


## `get_canonical`

Return a dictionary of the Matrix `mat` in canonical form. If
    - ``n`` is total number of states in a given transition matrix,
    - ``t`` is the number of transient states, and
    - ``a`` is the number of absorping states,

    then the orignal matrix (if it is an absorping discrete-time Markov Chain)
    can be reduced to canonical form,
    dictionary includes the following elements:

    - `:R`: upper-right submatrix of ``A`` (``t \\times a``)
    - `:A`: input transition matrix `mat` in canonical form (``n \\times n``)
    - `:inds`: indices in ``A`` of original states in `mat`
    - `:Q`: the transition matrix of transient states to transient states
        (dimension ``a \\times a``)
        - **NOTE** this is not the original input matrix
    - `:states_absorption`: all states in the input matrix that are absorbing
    - `:states_transient`: all states in the input matrix that are transient


###  Constructs

```
get_canonical(
    mat::SparseMatrixCSC{Float64,Int64} = transition_matrix_from_data_frame()
)
```

- Note: `transition_matrix_from_data_frame` returns `Q_default` with no arguments


## `get_fundamental_matrix`

Returns the fundamental matrix of ``Q``; requires `dict_canon` as input, which
    is generated using `get_canonical`. The fundamental matrix is

    ``F = (I - Q)^{-1}``

###  Constructs

```
get_fundamental_matrix(
	dict_canon::Dict
)
```


## `get_probability_of_absorption_state`

Returns the probability of absorption in each absorbing state ``j`` given
    starting state ``i``. Requires `dict_canon` as input, which is generated
    using `get_canonical`. The probability of the absorption state is calculated
    directly as

    ``FR``

    where ``F`` is the fundamental matrix (see `get_fundamental_matrix`) and ``R``
    is the ``t \times a`` upper-right submatrix of the canonical-form of the
    matrix `mat` used to generate `dict_canon`.


###  Constructs

```
get_probability_of_absorption_state(
	dict_canon::Dict
)
```


## `get_state_attribute`

Retrieve an attribute associated with a state

###  Constructs
```
get_state_attribute(
	state::Int64,
    attribute::String
)
```


## `get_state_entropies`

Get entropies associated with each state based on an input matrix `matrix`

###  Constructs

```
get_state_entropies(
	matrix::SparseMatrixCSC{Float64, Int64}
)
```

```
get_state_entropies()
```

- **Note** no argument will use the `MarkovModel` default matrix Q_default


## `find_absorption`

Return the ``t \\times a`` matrix ``B_{ij}`` that gives the probability that a
    walk starting in transient state ``i`` (`1 \\leq j \\leq a``, where``t`` is
    the number transient states) will end up in state ``j``
    (``1 \\leq j \\leq a``, where ``a`` is the number of absorping states)

###  Constructs

```
find_absorption(
	extraction_transition_states::Vector{Int64} = Vector{Int64}(),
	dict_variable_parameters::Dict{String, Float64} = configuration.dict_parameters_variable,
	configuration::Configuration = configuration,
	df_modspec::DataFrame = model_specification,
	extract_probabilites::Array{Tuple{Int64, Int64}, 1} = Array{Tuple{Int64, Int64}, 1}()
)
```


## `simulate_absorptions`

Simulate absorbptions over a DataFrame

###  Constructs

```
simulate_absorptions(
	df_runs::DataFrame,
	dict_cols_to_params::Dict{Symbol, String},
	extraction_transition_states::Array{Int64,1} = Array{Int64,1}(),
	extract_transition_matrix_indices::Array{Tuple{Int64, Int64}, 1} = Array{Tuple{Int64, Int64}, 1}()
)
```


## `swap_inds!`

Swap indices in sparse array indices over an arbitrary number of dimensions
    (vec_dim can be 1, 2, 3... etc vecs). This function supports `get_canonical`.

###  Constructs

```
swap_inds!(
	ind_1::Int64,
	ind_2::Int64,
	vecs_dim::Vector{Int64}...
)
```


## `swap_states!`

Swap states in sparse array indices over an arbitrary number of dimensions
    (vec_dim can be 1, 2, 3... etc vecs). This function supports `get_canonical`.

###  Constructs

```
swap_states!(
	state_1::Int64,
	state_2::Int64,
	vecs_dim::Vector{Int64}...
)
```


## `transition_matrix_from_data_frame!`

Read a data frame and convert it to a transition matrix. The data frame must
    have the following columns:
    - ``i``: row index in the matrix
    - ``j``: columns index in the matrix
    - ``expression``: the expression in the matrix

###  Constructs

```
transition_matrix_from_data_frame(
	dict_variables::Dict{String, Float64},
	df_specification::DataFrame = model_specification,
	transition_correction_threshold::Float64 = configuration.transition_correction_threshold,
	stochastic_orientation::Symbol = configuration.transition_stochastic_orientation,
	return_orientation::Symbol = :row,
	check_states::Bool = true
)
```

```
transition_matrix_from_data_frame(
	df_specification::DataFrame = model_specification,
	configuration::Configuration = configuration,
	return_orientation::Symbol = :row,
	check_states::Bool = true
)
```


## `do_walk`

Do a random walk on a transition matrix `mat_transition` (default is
    `Q_default`) starting from state `state_0`, and walk for `num_steps`.

    Set `save_trials` to `true` to return the random trials used to conduct the
    random walk.

###  Constructs

```
do_walk(
	state_0::Int64,
	num_steps::Int64,
	save_trials::Bool = false,
	mat_transition::SparseMatrixCSC{Float64,Int64} = Q_default,
	vec_state_times::Vector{Int64} = attribute_states.get_ordered_attribute(state_time_attribute)
)
```


## `find_next_state`

Based on a vector of transitions (e.g., a row in a row-stochastic matrix) and a
    random trial, return the next state.

###  Constructs

using a random trial, find the next state
```
find_next_state(
	vec_transitions::SparseVector{Float64,Int64},
	trial::Float64
)
```


## `simulate_walks`

Run a simulation of multiple walks on a transition matrix.

###  Constructs

```
simulate_walks(
    n_runs::Int64,
    n_time_periods::Int64,
    state_0::Int64;
    save_trials::Bool = false,
    return_type::Symbol = :DataFrame,
    mat_transition::SparseMatrixCSC{Float64,Int64} = Q_default,
    vec_state_times::Vector{Int64} = attribute_states.get_ordered_attribute(state_time_attribute),
    dict_entropy::Dict{Int64, Float64} = Dict{Int64, Float64}(),
    default_entropy::Float64 = 0.0
)
```

###  Function Arguements

- `n_runs`: number of walks to conduct
- `n_time_periods`: maximum number of time periods to walk
- `state_0`: starting state


###  Keyword Arguements

- `save_trials`: save off the random trials used to conduct the walk?
- `mat_transition`: specify a different transition matrix. Default is
    `Q_default`
- `vec_state_times`: vector of times spent in each state (ordered by
    `attribute_state`)
- `dict_entropy`: include a dictionary of entropies by state to characterize
    the entropy faced by a decisionmaker along the walk.
- `default_entropy`: default entropy if no entropy is specified


"""
struct MarkovModel

    model_specification::DataFrame
    attribute_states::AttributeTable
    configuration::Configuration
    check_states!
    check_transition_matrix!
    do_walk
    find_absorption
    find_next_state
    get_absorbing_states
    get_all_state_transitions
    get_canonical
    get_fundamental_matrix
    get_probability_of_absorption_state
    get_state_attribute
    get_state_entropies
    Q_default
    simulate_absorptions
    simulate_walks
    swap_inds!
    swap_states!
    transition_matrix_from_data_frame

    function MarkovModel(
        model_specification::DataFrame,
        attribute_states::AttributeTable,
        configuration::Configuration
    )

        ##  initialize some shared fields/parameters
        state_time_attribute = "time"
        state_outcome_value_attribute = "outcome_value"
        check_fields!(attribute_states.table, Symbol.([state_time_attribute, state_outcome_value_attribute]))

        ##  function to check states
        function check_states!(states::Vector{Int64})
            if !issubset(states, attribute_states.key_values)
                invalid = setdiff(states, attribute_states.key_values)
                error("Invalid values $(print_valid_values(invalid)) found by check_states!.")
            end
        end


        ##  function to check a transition matrix
        function check_transition_matrix!(
            mat::SparseMatrixCSC{Float64, Int64},
            correction_thresh::Float64 = configuration.transition_correction_thresh,
            stochastic_orientation::Symbol = configuration.transition_stochastic_orientation,
            return_orientation::Symbol = :row
        )
            #=
                correction_thresh is used to correct row sums if the devition is small. sets acceptable deviation from 1 for r
            =#
            if minimum(mat) < 0
                error("Negative values detected in the transition matrix: the minimum value is $(minimum(mat)).")
            end

            # get sums
            if stochastic_orientation == :row
                sums = sum.(eachrow(mat))
            elseif stochastic_orientation == :column
                sums = sum.(eachrow(transpose(mat)))
            else
                error("Invalid specification of stochastic_orientation '$(stochastic_orientation)': valid specifications are symbols :row or :column.")
            end

            # check deviation
            absdev = abs.(sums .- 1)
            if maximum(absdev) > correction_thresh
                w = findall(x -> (x > correction_thresh), absdev)
                print_inds = print_valid_values(w)
                error("Invalid transition probabilities: $(stochastic_orientation)s $w do not sum to 1 (and are beyond the threshold for correction).")
            elseif maximum(absdev) > 0
                # dim represents the dimension to sum over for normalization; 1 = rows, 2 = columns
                dim = (stochastic_orientation == :row) ? 2 : 1
                mat = normalize_sparse_matrix(mat, dim)
                if (return_orientation != (stochastic_orientation == :row))
                    mat = transpose(mat)
                end
            end

            return 0
        end

        function format_probs_as_columns(
            vec_ext_inds::Array{Tuple{Int64, Int64}, 1}
        )
            cols_out = ["" for x in 1:length(vec_ext_inds)]
            for i in 1:length(vec_ext_inds)
                cols_out[i] = "q_$(vec_ext_inds[i][1])-$(vec_ext_inds[i][2])"
            end
            return cols_out
        end

        ##  get the absorbing states of an MC
        function get_absorbing_states(
            mat::SparseMatrixCSC{Float64,Int64}
        )
            # row indices
            nz = findnz(mat)
            w = nz[1] .== nz[2]
            rows = nz[1][w]
            w2 = nz[3][w] .== 1
            #rows = findnz(mat)[1]
            #f(ind) = count(x -> (x == ind), rows)
            # row indices
            #sparse_inds_abs = findall(x -> (x == 1), f.(rows))

            if length(w2) > 0
                return sort(rows[w2])
            else
                return []
            end
        end


        ##  quickly return all state transitions
        function get_all_state_transitions(
            mat::SparseMatrixCSC{Float64,Int64}
        )
            return Set(zip(findnz(mat)[1:2]...))
        end

        ##  function to build cannonical form of MCMC
        function get_canonical(
            mat::SparseMatrixCSC{Float64,Int64} = transition_matrix_from_data_frame()
        )

            states_absorbing = get_absorbing_states(mat)
            mat_out = mat
            # get an index vector that tracks states
            vec_indices = collect(1:size(mat)[1])

            # only do this if there are absorbing states
            if length(states_absorbing) > 0

                # get matrix indices from the sparse array
                v_i, v_j, vals = findnz(mat)
                n_abs = length(states_absorbing)
                n = size(mat)[1]
                n_q = n - n_abs
                i = 0

                # swap states (in ascending order)
                while i < n_abs
                    #s1 = states_absorbing[i]
                    s1 = states_absorbing[n_abs - i]
                    s2 = n - i
                    swap_states!(s1, s2, v_i, v_j)
                    swap_inds!(s1, s2, vec_indices)
                    i += 1
                end

                # set the matrix anew
                mat_out = sparse(v_i, v_j, vals, n, n)
                # get Q and R
                Q = mat_out[1:n_q,1:n_q]
                R = mat_out[1:n_q,(n_q + 1):n]

                return Dict(:A => mat_out, :Q => Q, :R => R, :inds => vec_indices, :states_transient => vec_indices[1:n_q], :states_absorption => vec_indices[(n_q + 1):n])
            else
                return Dict(:A => mat_out)
            end
        end


        ##  get the fundamental matrix
        function get_fundamental_matrix(
            dict_canon::Dict
        )
            check_keys!(dict_canon, [:Q])
            return inv(Matrix(I - dict_canon[:Q]))
        end


        ##  get the matrix denoting the probability of starting in i and ending in state j
        function get_probability_of_absorption_state(
            dict_canon::Dict
        )
            check_keys!(dict_canon, [:Q, :R])
            mat_fundamental = get_fundamental_matrix(dict_canon)
            return mat_fundamental*dict_canon[:R]
        end


        ##  function to get attributes of states
        function get_state_attribute(
            state::Int64,
            attribute::String
        )
            #
            fm = "$(key)_to_$(attribute)"
            check_keys!(attribute_states.field_maps, [fm])
            return attribute_states.field_maps[fm][state]
        end


        ## calculate the state entropy of
        function get_state_entropies(
            matrix::SparseMatrixCSC{Float64, Int64}
        )
            # update the matrix with non-zero log entries
            fnz = findnz(matrix)
            matrix_log = sparse(fnz[1], fnz[2], log.(fnz[3]))
            return -1 .* sum(matrix .* matrix_log, dims = 2)[:, 1]
        end


        ##  find the absorbtiona as a function of variables denoted in the transition matrix
        function find_absorption(
            extraction_transition_states::Vector{Int64} = Vector{Int64}(),
            dict_variable_parameters::Dict{String, Float64} = configuration.dict_parameters_variable,
            configuration::Configuration = configuration,
            df_modspec::DataFrame = model_specification,
            extract_probabilites::Array{Tuple{Int64, Int64}, 1} = Array{Tuple{Int64, Int64}, 1}()
        )

            # get the row-oriented transition matrix
            mat = transition_matrix_from_data_frame(
                dict_variable_parameters,
                df_modspec,
                configuration.transition_correction_threshold,
                configuration.transition_stochastic_orientation
            )

            # get the canonical form and the absorbtion state
            dict_canon = get_canonical(mat)
            absorbs = get_probability_of_absorption_state(dict_canon)

            # get information on indexing
            n_transition = size(dict_canon[:Q])[1]
            transition_states = dict_canon[:inds][1:n_transition]
            absorbing_states = dict_canon[:inds][(n_transition + 1):end]

            # use extraction state to return only the probabilities for a list of states - default to all
            if length(extraction_transition_states) > 0
                inds = findall(x -> (x in extraction_transition_states), transition_states)
                absorbs = absorbs[inds, :]
                transition_states = transition_states[inds, :]
            end

            # extract raw probabilities?
            if length(extract_probabilites) > 0
                probs_out = -1 * ones(length(extract_probabilites))
                for prob_ind in 1:length(extract_probabilites)
                    mat_ind = extract_probabilites[prob_ind]
                    if (minimum(mat_ind) > 0) & (maximum(mat_ind) <= size(mat)[1])
                        probs_out[prob_ind] = mat[mat_ind[1], mat_ind[2]]
                    end
                end
            else
                probs_out = Array{Float64,1}()
            end

            return Dict(:B => absorbs, :states_absorb => absorbing_states, :states_transition => transition_states, :extraction_probabilities => probs_out)

        end


        ## simulate absorbptions over a dataframe
        function simulate_absorptions(
            df_runs::DataFrame,
            dict_cols_to_params::Dict{Symbol, String},
            extraction_transition_states::Array{Int64,1} = Array{Int64,1}(),
            extract_transition_matrix_indices::Array{Tuple{Int64, Int64}, 1} = Array{Tuple{Int64, Int64}, 1}()
        )

            check_fields!(df_runs, Symbol.(collect(keys(dict_cols_to_params))))

            # initialize output dataframe
            init_out_q = true
            df_out = DataFrame()
            absorbing_states = []
            transition_states = []
            ext_probs = []
            n_ep = 0
            n_t = 0

            for i in 1:nrow(df_runs)
                # get values
                dict_vals = Dict(zip([dict_cols_to_params[x] for x in keys(dict_cols_to_params)], collect(df_runs[i, Symbol.(collect(keys(dict_cols_to_params)))])))
                dict_abs = find_absorption(extraction_transition_states, dict_vals, config, model_specification, extract_transition_matrix_indices)
                # initialie data output structure
                if init_out_q
                    absorbing_states = dict_abs[:states_absorb]
                    n_t = length(dict_abs[:states_transition])
                    n_ep = length(dict_abs[:extraction_probabilities])
                    df_out = zeros(n_t*nrow(df_runs), length(dict_abs[:states_absorb]))
                    transition_states = Int64.(ones(size(df_out)[1]))
                    ext_probs = zeros(size(df_out)[1], n_ep)
                    init_out_q = false
                end

                # checkt that states line up—i.e., that we have not introduced any new absorbing states
                if dict_abs[:states_absorb] == absorbing_states
                    rng = ((i - 1)*n_t + 1):(i*n_t)
                    df_out[rng, :] = dict_abs[:B]
                    transition_states[rng] = dict_abs[:states_transition]
                    if n_ep > 0
                        ext_probs[rng, :] = transpose(hcat([dict_abs[:extraction_probabilities] for x in rng]...))
                    end
                else
                    error("Error in simulate_absorptions: absorbing states different from iteration $(i - 1) to $(i). Check the model definition file.")
                end
            end

            df_out = DataFrame(df_out, Symbol.(absorbing_states))
            df_out[:, :state_i] = transition_states
            if n_ep > 0
                cols = format_probs_as_columns(extract_transition_matrix_indices)
                df_ep = DataFrame(ext_probs, Symbol.(cols))
                df_out = hcat(df_out, df_ep)
            end

            return df_out
        end


        ##  function to swap indices in sparse array indices over an arbitrary number of dimensions (vec_dim can be 1, 2, 3... etc vecs)
        function swap_inds!(
            ind_1::Int64,
            ind_2::Int64,
            vecs_dim::Vector{Int64}...
        )
            #
            for vec in vecs_dim
                v1 = vec[ind_1]
                v2 = vec[ind_2]
                # reassign positions for state 1 to state 2 and vic-versa
                vec[ind_1] = v2
                vec[ind_2] = v1
            end
        end


        ##  function to swap states in sparse array indices over an arbitrary number of dimensions (vec_dim can be 1, 2, 3... etc vecs)
        function swap_states!(
            state_1::Int64,
            state_2::Int64,
            vecs_dim::Vector{Int64}...
        )
            #
            for vec in vecs_dim
                # get positions of state_1 and state_2 in vec_i
                all_1 = findall(x -> (x == state_1), vec)
                all_2 = findall(x -> (x == state_2), vec)
                # reassign positions for state 1 to state 2 and vic-versa
                vec[all_1] .= state_2
                vec[all_2] .= state_1
            end
        end


        ##  read a data frame and convert it to a transition matrix
        function transition_matrix_from_data_frame(
            dict_variables::Dict{String, Float64},
            df_specification::DataFrame = model_specification,
            transition_correction_threshold::Float64 = configuration.transition_correction_threshold,
            stochastic_orientation::Symbol = configuration.transition_stochastic_orientation,
            return_orientation::Symbol = :row,
            check_states::Bool = true
        )
            #=
                transition_correction_threshold gives the threshold for correcting sums in a stochastic matrix,
                    as small deviations may be acceptable or due to numerical error
                stochastic_orientation gives the orientation of the input matrix
                return_orientation gives the orientation of the cleaned matrix
            =#

            ##  CHECKS

            # check fields and get df copy
            check_fields!(df_specification, [:i, :j, :expression])
            df_to_transition = copy(df_specification)
            # check restrcted var values
            if transition_correction_threshold < 0
                error("Invalid transition_correction_threshold value $(transition_correction_threshold): the transition correction threshold should be positive")
            end

            ##  CLEAN MATRIX AND REPLACE VARIBALES

            df_to_transition[!, :i] = try_parse_int.(df_to_transition[!, :i])
            df_to_transition[!, :j] = try_parse_int.(df_to_transition[!, :j])

            # check for destination states that are not specified as sources
            missing_source_states = Int64.(setdiff(df_to_transition[!, :j], df_to_transition[!, :i]))
            if length(missing_source_states) > 0
                if check_states
                    df_append = DataFrame(:i => missing_source_states, :j => missing_source_states, :expression => ones(length(missing_source_states)))
                    df_to_transition = vcat(df_to_transition, df_append)
                else
                    vals = print_valid_values(missing_source_states)
                    error("Invalid destination states: states $(vals) were found as destinations, but are not specified as sources. Set check_states = true in transition_matrix_from_data_frame to assume these are absorbing states.")
                end
            end

            # replace the expression
            sym = :expression
            df_to_transition[!, sym] = repl_str_from_dict(string.(df_to_transition[!, sym]), dict_variables)
            vec_floats = try_parse_float.(eval.(Meta.parse.(df_to_transition[!, sym])))
            if any(missing .=== vec_floats)
                inds = findall(x -> (x === missing), vec_floats)
                vals = print_valid_values(df_to_transition[inds, sym])
                error("Invalid variable specification in df_specification: unable to replace variables $(vals) in column $(string(sym)).")
            else
                df_to_transition[!, sym] = vec_floats
            end

            # drop invalid rows.
            df_to_transition = filter(x -> ((x[:i] != 0) & (x[:j] != 0)), df_to_transition)


            ##  BUILD SPARSE ARRAY

            # get dimension, then build probabilities + sparse matrix
            n = maximum([maximum(df_to_transition[:, :i]), maximum(df_to_transition[:, :j])])
            #vec_probs = df_to_transition[:, :b] + prob_p*df_to_transition[:, :m]
            vec_probs = df_to_transition[:, :expression]
            mat_return = sparse(df_to_transition[:, :i], df_to_transition[:, :j], vec_probs, n, n)

            # check sums
            check_transition_matrix!(mat_return, transition_correction_threshold, stochastic_orientation)

            return mat_return
        end

        """
            get a transition matrix from a data frame(empty uses default variable values set in config; user can set following arguments:
            dict_variables::Dict{String, Float64}, - e.g., {"A": 0.18, "p" => 0.135}
            df_specification::DataFrame = model_specification,
            transition_correction_threshold::Float64 = configuration.transition_correction_threshold,
            stochastic_orientation::Symbol = configuration.transition_stochastic_orientation,
            return_orientation::Symbol = :row
        """

        # use mutliple dispatch to allow configuration as input
        function transition_matrix_from_data_frame(
            df_specification::DataFrame = model_specification,
            configuration::Configuration = configuration,
            return_orientation::Symbol = :row,
            check_states::Bool = true
        )

            dict_variables = configuration.dict_parameters_variable
            stochastic_orientation = configuration.transition_stochastic_orientation
            transition_correction_threshold = configuration.transition_correction_threshold

            return transition_matrix_from_data_frame(dict_variables, df_specification, transition_correction_threshold, stochastic_orientation, return_orientation, check_states)
        end


        ##  FUNCTIONS FOR SIMULATING A RANDOM WALK INFORMED BY THE MATRIX

        # define the default data frame
        Q_default = transition_matrix_from_data_frame()
        # execute a walk
        function do_walk(
            state_0::Int64,
            num_steps::Int64,
            save_trials::Bool = false,
            mat_transition::SparseMatrixCSC{Float64,Int64} = Q_default,
            vec_state_times::Vector{Int64} = attribute_states.get_ordered_attribute(state_time_attribute)
        )

            # get some components
            states_abs = get_absorbing_states(mat_transition)

            # intiialize output and iterator
            set_transitions = [(0, 0) for x in 1:num_steps]
            vec_out = Int64.(zeros(num_steps))
            vec_state_changes = Int64.(zeros(num_steps))
            vec_trials = save_trials ? zeros(num_steps) : Vector{Float64}()
            cur_state = state_0
            prev_state = state_0
            t_period_absorb = -1

            i = 1

            while i <= num_steps
                # have to spend at least one time period in the state; assume integer
                if (cur_state in states_abs)
                    t = num_steps - i + 1
                    t_period_absorb = i
                else
                    t = max(vec_state_times[cur_state], 1)
                end

                # update number of changes (only do at beginning to ensure we don't count a state change that puts i > t) and the states -  only care about transient states
                vec_state_changes[i] = (i > 1) ? ((cur_state != prev_state) ? 1 : 0) : 1
                set_transitions[i] = (cur_state != prev_state) ? (prev_state, cur_state) : (0, 0)
                vec_out[i:minimum([i + t - 1, num_steps])] .= cur_state

                # get next state
                trial = Random.rand()
                if save_trials & (i + t - 1 < num_steps)
                    vec_trials[i + t - 1] = trial
                end

                prev_state = cur_state
                cur_state = find_next_state(mat_transition[cur_state, :], trial)

                i += t
            end

            # convert the set of transitions to a set and drop
            set_transitions = setdiff(Set(set_transitions), Set([(0, 0)]))
            # return the output walk + the trials
            return vec_out, vec_state_changes, vec_trials, set_transitions, t_period_absorb
        end


        # using a random trial, find the next state
        function find_next_state(
            vec_transitions::SparseVector{Float64,Int64},
            trial::Float64
        )
            inds, vals = findnz(vec_transitions)
            probs = [0; cumsum(vals)]
            ind_state_next = findall(i -> ((trial <= probs[i]) & (trial > probs[i - 1])), collect(2:length(probs)))[1]
            return inds[ind_state_next]
        end


        # simulate multiple walks, return a dataframe
        function simulate_walks(
            n_runs::Int64,
            n_time_periods::Int64,
            state_0::Int64;
            save_trials::Bool = false,
            return_type::Symbol = :DataFrame,
            mat_transition::SparseMatrixCSC{Float64,Int64} = Q_default,
            vec_state_times::Vector{Int64} = attribute_states.get_ordered_attribute(state_time_attribute),
            dict_entropy::Dict{Int64, Float64} = Dict{Int64, Float64}(),
            default_entropy::Float64 = 0.0
        )

            check_states!([state_0])
            n_runs = maximum([n_runs, 1])
            n_time_periods = maximum([n_time_periods, 1])

            include_entropy_q = length(dict_entropy) > 0
            # initialize vectors that are long by runs and time period
            vec_out_trials = save_trials ? zeros(n_runs*n_time_periods) : Vector{Float64}()
            vec_out = Int64.(zeros(n_runs*n_time_periods, 4))
            # initialize vectors that are long by run
            vec_abs_info = Int64.(zeros(n_runs))
            vec_state_change_info = Int64.(zeros(n_runs))
            all_state_transitions = Set(zip(findnz(mat_transition)[1:2]...))
            vec_out_state_transitions = [all_state_transitions for x in 1:n_runs]
            # array of entropy values (total by time, total by steps, mean by time, mean by steps)
            array_entropy = include_entropy_q ? zeros(n_runs, 4) : 0

            time_periods = collect(1:n_time_periods)

            t0 = time()


            for i in 1:n_runs

                walk, state_changes, trials, set_transitions, time_abs = do_walk(
                    state_0,
                    n_time_periods,
                    save_trials,
                    mat_transition,
                    vec_state_times
                )

                rng = ((i - 1) * n_time_periods + 1):(i*n_time_periods)

                # update output vectors
                vec_out[rng, 1] .= i
                vec_out[rng, 2] = time_periods
                vec_out[rng, 3] = walk
                vec_out[rng, 4] = state_changes

                # update some vectors based on run (single value per run) - note that state_changes contains a 1 in the first entry, so we do not count that as a change
                vec_abs_info[i] = time_abs
                vec_state_change_info[i] = Int64(sum(state_changes) - 1)
                vec_out_state_transitions[i] = set_transitions

                if include_entropy_q
                    vec_entropies = get.((dict_entropy, ), walk, default_entropy)

                    if time_abs != 1
                        rng = (time_abs > 1) ? (1:(time_abs - 1)) : (1:length(vec_entropies))
                        # get time-weighted average entropy until absorption and entropy normed by number of state changes - assumes that initial state is a transient state
                        # column ordering is (1) total by time, (2) total by unique changes (3) mean by time, (4) mean by steps
                        array_entropy[i, 1] = sum(vec_entropies[rng])
                        array_entropy[i, 2] = dot(vec_entropies[rng], state_changes[rng])
                        array_entropy[i, 3] = array_entropy[i, 1]/(time_abs - 1)
                        array_entropy[i, 4] = array_entropy[i, 2]/sum(state_changes[rng])
                    elseif time_abs == 1
                        # should be 0 -- this case only occurs if you start in an absorbing state
                        array_entropy[i, 1] = vec_entropies[1]
                        array_entropy[i, 2] = vec_entropies[1]
                        array_entropy[i, 3] = vec_entropies[1]
                        array_entropy[i, 4] = vec_entropies[1]
                    end
                end

                if save_trials
                    vec_out_trials[rng] = trials
                end
            end


            if return_type == :Arrays

                return vec_out, vec_out_trials, array_entropy, vec_abs_info, vec_state_change_info, vec_out_state_transitions

            elseif return_type in [:DataFrame, :DataFrames]
                # build data frame out
                df_out = DataFrame(vec_out, [:run_id, :time_period, :state, :new_state])
                if save_trials
                    df_out[!, :trials] = vec_out_trials
                end

                if return_type == :DataFrame
                    return df_out
                elseif return_type == :DataFrames
                    # add absorption info
                    df_abs = DataFrame(
                        :run_id => collect(1:n_runs),
                        :absorption_time_period => vec_abs_info,
                        :absorbed => (vec_abs_info .!= -1),
                        :n_state_changes => vec_state_change_info
                    )
                    # add entropy information if needed
                    if include_entropy_q
                        df_abs = hcat(
                            df_abs,
                            DataFrame(array_entropy, [:entropy_total_by_time, :entropy_total_by_steps, :entropy_mean_by_time, :entropy_mean_by_steps])
                        )
                    end

                    return df_out, df_abs, vec_out_state_transitions
                end
            end
        end


        ##  some multiple dispatch definitions

        function get_all_state_transitions()
            return get_all_state_transitions(Q_default)
        end

        function get_state_entropies()
            return get_state_entropies(Q_default)
        end

        return new(
            model_specification,
            attribute_states,
            configuration,
            check_states!,
            check_transition_matrix!,
            do_walk,
            find_absorption,
            find_next_state,
            get_absorbing_states,
            get_all_state_transitions,
            get_canonical,
            get_fundamental_matrix,
            get_probability_of_absorption_state,
            get_state_attribute,
            get_state_entropies,
            Q_default,
            simulate_absorptions,
            simulate_walks,
            swap_inds!,
            swap_states!,
            transition_matrix_from_data_frame
        )

    end
end


end
