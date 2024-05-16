# load support functions and data structures

module MDUSetups

using DataStructures
using SupportFunctions


export config
export dir_datasets
export dir_julia
export dir_out
export dir_proj
export dir_out
export fp_config_complexity
export fp_default_entropies
export fp_default_expected_outcome_by_edge
export fp_default_outcomes_out
export fp_default_outcomes_summary_out
export fp_default_simulation_walks_out
export fpt_default_analytical_convergence

#########################
#    SET DIRECTORIES    #
#########################

##  high level directory structure
dir_proj = dirname(dirname(@__FILE__))
@info("Setting project directory to '$(dir_proj)'")
dir_julia = check_path(joinpath(dir_proj, "julia"), false)
dir_out = check_path(joinpath(dir_proj, "out"), true)
dir_ref = check_path(joinpath(dir_proj, "ref"), false)


##  dependent directories
dir_datasets = check_path(joinpath(dir_ref, "datasets"), true)


##  get configuration and do some data cleaning on the dictionary
fp_config_complexity = check_path(joinpath(dir_proj, "complexity_modeling.config"), false)
fp_default_entropies = joinpath(dir_out, "state_entropies.csv")
fp_default_expected_outcome_by_edge = joinpath(dir_out, "expected_outcome_by_edge.csv")
fp_default_outcomes_out = joinpath(dir_out, "outcome_values.csv")
fp_default_outcomes_summary_out = joinpath(dir_out, "outcome_values_summary.csv")
fp_default_simulation_walks_out = joinpath(dir_out, "simulated_walks.csv")
fpt_default_analytical_convergence(x::Int64, field_scenario_id::Symbol) = joinpath(dir_out, "analytical_solutions_$(field_scenario_id)-$(x).csv")


##  initialize the configuration file
config = Configuration(fp_config_complexity)

end
