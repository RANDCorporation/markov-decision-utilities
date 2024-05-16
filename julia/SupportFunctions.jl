"""
SupportFunctions includes a number of basic functions that facilitiate other
packages.

Includes:

    build_dict
    check_fields!
    check_keys!
    check_path
    check_trailing_num_field
    export_sub_file_csv
    get_table
    normalize_sparse_matrix
    parse_config
    print_header
    printnames
    print_valid_values
    read_csv
    repl_str_from_dict
    track_timer_db
    try_parse_float
    try_parse_int
    zip_cols
    get_rand_blocks
    rand_sample
    sort_row_ids
    write_and_compress_output

Use ?function to explore the function.
"""
module SupportFunctions

using DataFrames
using CSV
using Random
using SparseArrays

export build_dict
export check_fields!
export check_keys!
export check_path
export check_trailing_num_field
export export_sub_file_csv
export get_table
export normalize_sparse_matrix
export parse_config
export print_header
export printnames
export print_valid_values
export read_csv
export repl_str_from_dict
export track_timer_db
export try_parse_float
export try_parse_int
export zip_cols
export get_rand_blocks
export rand_sample
export sort_row_ids
export write_and_compress_output




##  FUNCTION build_dict IS USED TO QUICKLY BUILD DICTIONARIES FROM DATAFRAMES
"""
Build a dictionary from a data frame `df_in`. Maps the first column to an
    individual element (if there is only one additional column) or a tuple
    of row-wise elements (if there is more than one column.)

```
build_dict(
    df_in::DataFrame
)::Dict
```

"""
function build_dict(df_in::DataFrame)::Dict

    sz = size(df_in)
    n = sz[2]

    if n == 2
        #in two columns, simply map a to b
        dict_out = Dict([x for x in zip(df_in[:, 1], df_in[:, 2])])
    else
        #if there are more than two columns, leading columns are tuples
        fields_key = names(df_in)[1:(n - 1)]
        tups = [tuple(x...) for x in eachrow(Matrix(df_in[:, fields_key]))]
        dict_out = Dict([x for x in zip(tups, df_in[:, names(df_in)[n]])])
    end

    return dict_out
end


##  functions for checking fields/keys
"""
Check to see if `fields_check` are in DataFrame `data_frame`. Throws an
    error if any of the fields are not present.

## Constructs
```
check_fields!(
    data_frame::DataFrame,
    fields_check::Array{Symbol, 1}
)::Nothing
```

"""
function check_fields!(
    data_frame::DataFrame,
    fields_check::Array{Symbol, 1}
)::Nothing
    # verify fields
    fields_req = Set(String.(fields_check))

    if !issubset(fields_req, Set(names(data_frame)))
        str_missing_fields = setdiff(fields_req, Set(names(data_frame)))
        str_missing_fields = join("'".*str_missing_fields.*"'", ", ")
        error("Fields $(str_missing_fields) not found. Check the data frame.")
    end
end

##  check to ensure keys are present in the dictionary
"""
Check to see if `keys_check` are present in dictionary `dict`. Throws an
    error if any of the keys are not present.

    Optionally, you can specify `string_notify`, which repkaces the substring
    "##MISSINGHERE##" with missing keys, to specify a message for the error.

## Constructs

```
check_keys!(
    dict::Dict,
    keys_check::Array{String, 1},
    string_notify::String = "Keys ##MISSINGHERE## not found. Check the dictionary."
)::Nothing
```

```
check_keys!(
    dict::Dict,
    keys_check::Array{Symbol, 1},
    string_notify::String = "Keys ##MISSINGHERE## not found. Check the dictionary."
)::Nothing
```
"""
function check_keys!(
    dict::Dict,
    keys_check::Array{String, 1},
    string_notify::String = "Keys ##MISSINGHERE## not found. Check the dictionary."
)::Nothing
    # verify fields
    keys_req = Set(String.(keys_check))
    keys_avail = Set(String.(collect(keys(dict))))

    if !issubset(keys_req, keys_avail)
        str_missing_keys = setdiff(keys_req, keys_avail)
        str_missing_keys = join("'".*str_missing_keys.*"'", ", ")
        error(replace(string_notify, "##MISSINGHERE##" => str_missing_keys))
    end
end

function check_keys!(
    dict::Dict,
    keys_check::Array{Symbol, 1},
    string_notify::String = "Keys ##MISSINGHERE## not found. Check the dictionary."
)::Nothing
    # verify fields
    keys_req = Set(String.(keys_check))
    keys_avail = Set(String.(collect(keys(dict))))

    if !issubset(keys_req, keys_avail)
        str_missing_keys = setdiff(keys_req, keys_avail)
        str_missing_keys = join("'".*str_missing_keys.*"'", ", ")
        error(replace(string_notify, "##MISSINGHERE##" => str_missing_keys))
    end
end


"""
Check to see if `path` exists. If it does not, can specify `create_directory_q`
    to create the path if the path is a directory. Throws an error if the path
    does not exist and is not to be created.

## Constructs

```
check_path(
    path::String,
    create_directory_q::Bool
)::Nothing
```
"""
##  use check_path to return a path if it exists and throw an error otherwise; optional 'create_directory_q' can be used to create a directory
function check_path(
    path::String,
    create_directory_q::Bool
)::String
    if ispath(path)
        return path
    else
        if create_directory_q & !occursin(".", basename(path))
            mkdir(path)
            print("Created directory '$(path)'\n")
            return path
        else
            error("Path '$(path)' not found.")
        end
    end
end


##  FUNCTION check_trailing_num_field checks if the end of a field (split by "_") is numeric; allows for the return of that number
"""
Check a field for a trailing number (should be revised to use regular
    expressions). Returns type:

    * `Bool` if `return_value == "query"`, where the bool determines whether
        or not the field has a trailing number
    * `Int` or `Nothing` if `return_value == "num"` (`Int` if the number is
        found, `Nothing` is )

## Constructs

```
check_trailing_num_field(
    field::String,
    return_val::String = "query",
    field_delim::String = "_"
)::Union{Nothing, Bool, Int}
```
"""
function check_trailing_num_field(
    field::String,
    return_val::String = "query",
    field_delim::String = "_"
)
    nv = split(field, field_delim)
        tp = nothing
    if length(nv) > 1
        tp = tryparse(Int64, nv[length(nv)])
        rv_bool = (tp !== nothing)
    else
        rv_bool = false
    end

    if return_val == "query"
        return rv_bool
    elseif return_val == "num"
        return tp
    end
end



##  FUNCTION FOR EXPORTING A SUBSET OF A DATAFRAME AND CONVERTING IT TO A TAR GZ
function export_sub_file_csv(df_in::DataFrame, fields_save::Array{String, 1}, fp_csv_exp::String, conv_gz_q::Bool = true, rm_csv_q::Bool = false)

    print("\n\n###    Exporting drop fields to '$(fp_csv_exp)'... \n")

    #export unused data to CSV
    CSV.write(fp_csv_exp, df_in[:, fields_save])
    print("\n###    Export complete. ")

    if conv_gz_q
        print("Creating compressed file...\n")
        fp_targz_exp = replace(fp_csv_exp, ".csv" => ".tar.gz")
        fn_targz_exp = basename(fp_targz_exp)

        #remove existing tar.gz
        if ispath(fp_targz_exp)
            rm(fp_targz_exp)
        end

        #build tar.gz file
        dir_cur_tmp = pwd()
        comm_tgz = `tar -cvzf "$(fn_targz_exp)" "$(basename(fp_csv_exp))"`
        cd(dirname(fp_csv_exp))
        run(comm_tgz)
        cd(dir_cur_tmp)

        print("\n###    Compressed file complete.")

        if rm_csv_q
            print(" Removing CSV...\n")
            #remove the csv
            rm(fp_csv_exp)
        end

        print("\n###    Done.\n\n")
    end
end



##  define some functions for retrieving tables
function get_table(mode::String, fp_csv::String, query::String)

    # switch between csv/sql quickly with the "mode"
    if mode == "csv"
        # assume that the header is available
        df_in = read_csv(fp_csv, true)
    elseif mode == "sql"
        # placeholder for integration with sql
        df_in = DataFrame()
    else
        error("Mode $(mode) not found for get_table. Check to ensure it is specified as 'csv' or 'sql'.")
    end

    return df_in
end


## normalize a sparse matrix by row or column sums
function normalize_sparse_matrix(matrix::SparseMatrixCSC{Float64, Int64}, dim::Int64)
    # dim = 1 -> columns are normalized
    # dim = 2 -> rows are normalized
    if dim in [1, 2]
        sum_vec = sum(matrix, dims = dim)
        m, n = size(matrix)
        i, j, v = findnz(matrix)

        if dim == 1
            for v_ind in 1:length(v)
                v[v_ind] = sum_vec[j[v_ind]] == 0 ? 0.0 : v[v_ind]/sum_vec[j[v_ind]]
            end
        elseif dim == 2
            for v_ind in 1:length(v)
                v[v_ind] = sum_vec[i[v_ind]] == 0 ? 0.0 : v[v_ind]/sum_vec[i[v_ind]]
            end
        end

        return dropzeros(sparse(i, j, v, m, n))
    else
        error("Invalid matrix dimension $(dim) for normalization. Use '1' (columns) or '2' (rows).")
    end
end


## use parse_config to read a configuration file and return a dictionary
function parse_config(fp_config::String)

    if ispath(fp_config)
        # read lines from the file
        file_read = open(fp_config, "r")
        lines = readlines(file_read)
        close(file_read)

        # initialize the output dictionary
        dict_out = Dict{String, Any}()
        for i in 1:length(lines)
            line = lines[i]
            if length(line) > 0
                # check conditions
                if string(line[1]) != "#"
                    # drop comments
                    line = split(line, "#")[1]
                    # next, split on colon
                    line = split(line, ":")
                    key = strip(string(line[1]))
                    value = strip(string(line[2]))

                    # split on commas and convert to vector
                    value = split(value, ",")
                    # check first value
                    if tryparse(Float64, value[1]) != nothing
                        # convert to numeric
                        value = tryparse.(Float64, value)
                        value = filter(x -> (x != nothing), value)
                        # check for integer conversion
                        if round.(value) == value
                            value = Int64.(value)
                        end
                    else
                        # convert to character and return
                        value = string.(strip.(value))
                    end

                    # remove arrays and convert to boolean if applicable
                    if length(value) == 1
                        value = value[1]
                        if lowercase(string(value)) in ["false", "true"]
                            value = (lowercase(string(value)) == "true")
                        end
                    end

                    dict_out[key] = value
                end

            end
        end

        return dict_out
    else
        print("\nFile $(fp_config) not found. Please check for the configuration file.")
        return(-1)
    end
end



##  FUNCTION FOR PRINTING NICE HEADERS
function print_header(str_in::String)
    nchar = length(str_in)
    char_div = "#"

    n_space = 4
    n_buffer = 3

    buffer = char_div^n_buffer
    spacer = " "^n_space

    bridge = char_div^(nchar + 2*(n_space + n_buffer))
    river = buffer*(" "^(nchar + 2*n_space))*buffer
    txt = buffer*spacer*str_in*spacer*buffer

    div = "\n"*join([bridge, river, txt, river, bridge], "\n")*"\n"

    print(div)
end


##  FUNCTION printnames USED TO PRINT LONG LISTS QUICKLY
function printnames(list::Array{String})
    for nm in list
        print("$(nm)\n")
    end
end


##  print_valid_values can be used to print a list of values into a string
function print_valid_values(array::Array{String, 1})
    str_valid = join("'".*(string.(array)).*"'", ", ", ", and ")
    return str_valid
end
function print_valid_values(array::Array{Int64, 1})
    str_valid = join("'".*(string.(array)).*"'", ", ", ", and ")
    return str_valid
end
function print_valid_values(array::Array{Symbol, 1})
    str_valid = join("'".*(string.(array)).*"'", ", ", ", and ")
    return str_valid
end



##  FUNCTION read_csv IS USED TO READ IN CSVs TO DATAFRAMES
function read_csv(path_in::String, header_q::Bool, lim::Int64=-1, names_vec::Array{Symbol, 1} = Array{Symbol, 1}(), delimiter::String = ",")

    if lim == -1
        #get some data
        df_in = CSV.File(
            path_in,
            delim = delimiter,
            ignorerepeated = false,
            header = header_q
        ) |> DataFrame
    else
        #get some data
        df_in = CSV.File(
            path_in,
            delim = delimiter,
            ignorerepeated = false,
            header = header_q;
            limit = lim
        ) |> DataFrame
    end

    if (!header_q & (length(names_vec) == size(df_in)[2]))
        dict_rnm = Dict([x for x in zip(Symbol.(names(df_in)), names_vec)])
        rename!(df_in, dict_rnm)
    end

    return df_in
end


##  replace elements of a string using a dictionary
function repl_str_from_dict(str::String, dict_repl::Dict{String, String})
    str_out = str
    for k in keys(dict_repl)
        str_out = replace(str_out, k => dict_repl[k])
    end
    return str_out
end

function repl_str_from_dict(str::String, dict_repl::Dict{String, Float64})
    str_out = str
    for k in keys(dict_repl)
        str_out = replace(str_out, k => string(dict_repl[k]))
    end
    return str_out
end

function repl_str_from_dict(str::String, dict_repl::Dict{String, Int64})
    str_out = str
    for k in keys(dict_repl)
        str_out = replace(str_out, k => string(dict_repl[k]))
    end
    return str_out
end

function repl_str_from_dict(str::String, dict_repl::Dict{String, Symbol})
    str_out = str
    for k in keys(dict_repl)
        str_out = replace(str_out, k => string(dict_repl[k]))
    end
    return str_out
end

function repl_str_from_dict(str::Array{String, 1}, dict_repl::Dict)
    strs_out = copy(str)
    for i in 1:length(strs_out)
        strs_out[i] = repl_str_from_dict(strs_out[i], dict_repl)
    end
    return strs_out
end

##  FUNCTION FOR TRACKING TIME TO BUILD DB
function track_timer_db(checkpoint::Int64, time_0::Float64)
    print("\n Checkpoint (" * string(checkpoint)* ") complete at about " * string(round(time() - time_0)) * " seconds.\n\n")
end

# use try_parse_float to attempt to parse objects into a float; use multiple dispatch to avoid conditionals
function try_parse_float(x::String)
    tp = tryparse(Float64, x)
    if tp == nothing
        return missing
    else
        return tp
    end
end
function try_parse_float(x::Real)
    return Float64(x)
end
function try_parse_float(x::Missing)
    return missing
end

# use try_parse_int to attempt to parse an integer; if it's an integer, nothing happens; if it's float, it gets rounded; it it's a string, use tryparse

function try_parse_int(x::Int64)
    return x
end
function try_parse_int(x::Real)
    return Int64(round(x))
end
function try_parse_int(x::Missing)
    return missing
end
# string depends on other dispatches
function try_parse_int(x::String)
    tp = try_parse_float(x)
    tp = try_parse_int(tp)
    if (tp == nothing) | (tp === missing)
        return missing
    else
        return tp
    end
end



# use zip_cols to zip columns of a data frame or matrix together for iteration
function zip_cols(mat::DataFrame)
    return zip([mat[:, i] for i in 1:size(mat)[2]]...)
end
function zip_cols(mat::Matrix)
    return zip([mat[:, i] for i in 1:size(mat)[2]]...)
end



#############################################################
#    FUNCTIONS FOR SAMPLING AND OPTIMIZATION EXPERIMENTS    #
#############################################################

## FUNCTION get_rand_blocks SELECTS RANDOM BLOCKS (STRATIFIED BY POPULATION) TO USE TO ENFORCE CONSTRAINTS
function get_rand_blocks(vec::Array{Float64, 1}, q_vec::Array{Float64, 1} = [0, 0.25, 0.5, 0.75, 1], prob::Float64 = 0.1, min_num::Int64 = 5)

    # get the quantiles for "stratification"
    q = reverse(quantile(vec, sort(q_vec)))
    #add upper limit to ensure sampling is clean
    q[1] = q[1] + 1
    n_q = length(q) - 1

    # set the number of samples
    n = Int64(minimum([maximum([min_num, prob*length(vec)]), length(vec)]))

    n_per = Int64(floor(n/n_q))
    n_extra = n%n_q

    samp_out = Array{Int64, 1}(zeros(n))

    ind = 1

    for i in range(1, n_q, step = 1)
        if i <= n_extra
            n_samp = n_per + 1
        else
            n_samp = n_per
        end


        samp_space = findall(x -> (x < q[i]) & (x >= q[i + 1]), vec)
        samp_sub = rand_sample(samp_space, n_samp)

        samp_out[ind:(ind + n_samp - 1)] = samp_sub

        ind = ind + n_samp
    end

    return samp_out
end



##  FUNCTION rand_sample GENERATES A RANDOM SUBSAMPLE OF vec OF LENGTH n
function rand_sample(vec::AbstractArray, n::Int64)
    i = 1

    m = minimum([n, length(vec)])
    vec_out = Array{Any, 1}(zeros(m))

    while i <= m
        vec_out[i] = Random.rand(vec, 1)[1]
        vec = [x for x in vec if x != vec_out[i]]
        i += 1
    end

    return vec_out
end

##  sort row ids
function sort_row_ids(df_to_sort::DataFrame, fields_id::Array{Symbol, 1})

    # check id fields
    check_fields!(df_to_sort, fields_id)
    fields_in = Symbol.(names(df_to_sort))

    # cut data frame in two and sort fields, then re-concatenate
    df_sort_ids = copy(df_to_sort[:, fields_id])
    df_sort_ids = map(x -> sort(x), eachrow(Matrix(df_sort_ids)))
    df_sort_ids = DataFrame(permutedims(hcat(df_sort_ids...)), fields_id)

    df_to_sort = hcat(select(df_to_sort, Not(fields_id)), df_sort_ids)

    return df_to_sort[:, fields_in]

end


# write output to CSV, then compress to targz with same base name
function write_and_compress_output(fp_out::String, df_out::DataFrame)

    # get some file name components
    dir_write = dirname(fp_out)
    fn_write = basename(fp_out)
    ext = occursin(".", fn_write) ? split(fn_write, ".")[end] : nothing
    fn_base = (ext == nothing) ? fn_write : replace(fn_write, ("."*ext) => "")

    # switch directories
    dir_cur = pwd()
    cd(dir_write)

    # write the output
    print("\nWriting output to $(fp_out)...")
    CSV.write(fn_write, df_out)
    print("Done.\n")

    # get the tarball name and compress
    fn_tar = "$(fn_base).tar.gz"
    print("\nCompressing to $(joinpath(dir_write, fn_tar))...")
    ENV["COPYFILE_DISABLE"] = 1
    exit_status = run(`tar -cvzf $(fn_tar) $(fn_write)`)
    print("Done.\n")

    # return back
    cd(dir_cur)

    return exit_status
end


end
