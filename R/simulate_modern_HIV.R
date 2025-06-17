#' Modern HIV simulator
#'
#' A modern agent-based simulation for HIV that explicitly models within-host simulations using biophybreak extensions.
#'
#' @param params A list of parameters used in the simulation. This includes:
#' - rate_function_parameters: A list of parameters needed to evaluate the rate
#'  function. Expected as a `list(`R0`=...)``.
#' - a: The initial amount of diversity of sequences per host.
#' - b: The amount of diversity of sequences to increase per host per unit time.
#' - mutation_rate: The mutation rate in per nt per month per sequence.
#' - minimum_population: The minimum population size.
#' - maximum_population_target: The maximum population size.
#' - total_steps_after_exp_phase: The total number of steps after the exponential phase.
#' @return list with keys "matrix" and "input_params".
#'   The "matrix" key contains the resulting pairwise distance matrix.
#' @export
#'
#' @examples
#' parameters <- list(
#'   "rate_function_parameters" = list("R0" = 5), "a" = 5, "b" = 5,
#'   "minimum_population" = 15, "maximum_population_target" = 500
#' )
#'
simulate_modern_HIV <- function(params) { # nolint: object_name_linter
  # This gives an example "workflow" function for this package.
  # Simple integration tests should use similar workflows to insert changes.

  # First, determine the rate function for offspring generation
  # This obtains the function used in [Graw et al. 2012], [Kupperman et al. 2022]
  biphasic_rate_function <- SEEPS::get_biphasic_HIV_rate(
    params = params[["rate_function_parameters"]]
  )

  # Forward pass to generate transmission history
  simulator_result <- SEEPS::gen_transmission_history_exponential_constant(
    minimum_population = params[["minimum_population"]],
    offspring_rate_fn = biphasic_rate_function,
    total_steps = params[["total_steps_after_exp_phase"]],
    maximum_population_target = params[["maximum_population_target"]],
    spike_root = FALSE
  )

  # Simulate contact tracing on the contact network to identify a sample
  target_sample <- SEEPS::contact_traced_uniform_ids(
    active = simulator_result[["active"]],
    parents = simulator_result[["parents"]],
    minimum_sample_size = params[["minimum_population"]],
    p = 0.9
  ) # 90% probability of uncovering each contact

  # Get a geneology compatable with biophybreak (bpb) layout.
  res <- SEEPS::reduce_transmission_history_bpb(
    samples = target_sample[["samples"]],
    parents = simulator_result$parents,
    current_step = simulator_result$t_end
  )

  # Provide values of a and b for the within-host diversity models.
  phylogeny <- SEEPS::geneology_to_phylogeny_bpb(
    transmission_history = res$parents,
    infection_times = res$transmission_times,
    sample_times = res$sample_times,
    a = params[["a"]], b = params[["b"]],
    leaf_sample_ids = res$transformed_sample_indices
  )

  # Convert the time signals into a # of mutations using a rate
  phylogeny <- SEEPS::stochastify_transmission_history(
    transmission_history = phylogeny$phylogeny,
    rate = params[["mutation_rate"]] / 12
  ) # Provide in rate per sequence per year

  # Only reconstruct the nodes corresponding to the sampled tips
  distance_matrix <- SEEPS::geneology_to_distance_matrix(
    geneology = phylogeny$geneology,
    spike_root = FALSE
  )

  # Since we used a density sampler, we may want to extract a grouping of closely related
  # individuals of a fixed size.
  distance_matrix <- SEEPS::reduce_large_matrix(
    oversampled_matrix = distance_matrix,
    subsample_size = params[["minimum_population"]],
    spike_root = FALSE
  )

  return(list("matrix" = distance_matrix$matrix))
}
