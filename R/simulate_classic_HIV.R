#' Classic HIV simulator for population level simulations.
#'
#' A simple agent-based simulation for HIV. Returns pairwise distance matrix samples
#' obtained from a simulated outbreak. This simulation models population level diversity,
#' but does not explicitly consider within-host diversity. A mutation clock rate, R0, and
#' several epidemiological parameters are required to describe the duration of infections.
#'
#' @return list with key "matrix"
#' @export
#'
#' @examples
#' parameters <- list("rate_function_parameters" <- list("R0"=5),
#'                     "minimum_population"=50, "maximum_population_target"=100)
#'
simulate_classic_HIV <- function(params) {  # nolint: object_name_linter
    # This gives an example "workflow" function for this package.
    # Simple integration tests should use similar workflows to test changes.

    # First, determine the rate function for offspring generation
    # This obtains the function used in [Graw et al. 2012], [Kupperman et al. 2022]
    biphasic_rate_function <- SEEPS::get_biphasic_HIV_rate(
        params = params[["rate_function_parameters"]])

    # Forward pass to generate transmission history
    simulator_result <- SEEPS::gen_transmission_history_exponential_constant(
        minimum_population=params[["minimum_population"]],
        offspring_rate_fn = biphasic_rate_function,
        total_steps = params[["total_steps_after_exp_phase"]],
        maximum_population_target = params[["maximum_population_target"]],
        spike_root = FALSE)
        # spike root option allows you to record distances to the initial infection

    # Next, determine the sample of ID's to subsample.
    # This function takes independent random samples, but you could do
    # something much more fancy here and change this call signature
    target_sample <- SEEPS::random_fixed_size_ids(
        active = simulator_result[["active"]],
        minimum_size = params[["minimum_population"]],
        spike_root = FALSE)

    # Obtain the transmission history for the subset of sampled individuals
    transmission_history <- SEEPS::reduce_transmission_history(
        samples = target_sample[["samples"]],
        parents = simulator_result[["parents"]],
        current_step = simulator_result[["t_end"]],
        spike_root = FALSE)

    # Convert time signals to # of mutations using a rate
    geneology <- SEEPS::stochastify_transmission_history(
        transmission_history=transmission_history$geneology,
        rate=params[["mutation_rate"]] / 12)  # Provide in rate per sequence per year

    # convert the geneology into a distance matrix
    distance_matrix <- SEEPS::geneology_to_distance_matrix(
        geneology=geneology$geneology,
        spike_root=FALSE)

    # Now return the result. Additional metadata used above can be collected
    # and added to the return list.
    return(list("matrix" = distance_matrix))
}
