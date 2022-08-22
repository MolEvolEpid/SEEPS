#' Modern HIV simulator
#'
#' A modern agent-based simulation for HIV that explicitly models within-host simulations using biophybreak extensions.
#'
#' @return list with keys "matrix" and "input_params"
#' @export
#'
#' @examples
#' parameters <- list("rate_function_parameters" = list("R0"=5), "a"=5, "b"=5,
#'                     "minimum_population"=15, "maximum_population_target"=500)
#'
simulate_modern_HIV <- function(params) {  # nolint: object_name_linter
    # This gives an example "workflow" function for this package.
    # Simple integration tests should use similar workflows to insert changes.

    # First, determine the rate function for offspring generation
    # This obtains the function used in [Graw et al. 2012], [Kupperman et al. 2022]
    biphasic_rate_function <- get_biphasic_HIV_rate(
        params = params[["rate_function_parameters"]])

    # Forward pass to generate transmission history
    simulator_result <- gen_transmission_history_exponential_constant(
        minimum_population=params[["minimum_population"]],
        offspring_rate_fn = biphasic_rate_function,
        total_steps = params[["total_steps_after_exp_phase"]],
        maximum_population_target = params[["maximum_population_target"]],
        spike_root = FALSE)

    target_sample <- random_prop_ids(
        active = simulator_result[["active"]],
        minimum_size = params[["minimum_population"]],
        proportion = 0.1,  # sample a 10% and report closest infections
        spike_root = FALSE)

    # Get a geneology compatable with biophybreak (bpp) layout.
    res <- reduce_transmission_history_bpb(
        samples = target_sample[["samples"]],
        parents = simulator_result$parents,
        current_step = simulator_result$t_end)

    # You do not need to provide a clock rate. Instead, provide values of
    # a and b for the within-host diversity models.
    phylogeny <- geneology_to_phylogeny_bpb(
        transmission_history = res$parents,
        infection_times = res$transmission_times,
        sample_times = res$sample_times,
        a = params[["a"]], b = params[["b"]])

    # Convert the time signals into a # of mutations using a rate
    phylogeny <- stochastify_transmission_history(
            transmission_history=phylogeny$phylogeny,
            rate=params[["mutation_rate"]] / 12)  # Provide in rate per sequence per year

    print(phylogeny)
    # Only reconstruct the nodes corresponding to the sampled tips
    distance_matrix <- geneology_to_distance_matrix_bpb(
        geneology = phylogeny$geneology,
        spike_root = FALSE)

    # Since we used a density sampler, we may want to extract a grouping of closely related
    # individuals of a fixed size.
    distance_matrix <- reduce_large_matrix(oversampled_matrix = distance_matrix,
        subsample_size = params[["minimum_population"]],
        spike_root = FALSE)
    return(list("matrix" = distance_matrix))
}
