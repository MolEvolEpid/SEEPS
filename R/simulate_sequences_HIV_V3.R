#' Modern HIV epidemic sequences simulator
#'
#' A modern agent-based simulation for HIV that explicitly models within-host simulations using biophybreak extensions.
#' This pre-build simulation returns a collection of sequences, rather than a matrix.
#' @return list with keys "matrix" and "input_params"
#' @export
#'
#' @examples
#' parameters <- list("rate_function_parameters" = list("R0"=5), "a"=5, "b"=5,
#'                     "mutation_rate" = 0.0067,  # units are per nt
#'                     "contact_tracing_discovery_probability" = 0.9,
#'                     "minimum_population"=15, "maximum_population_target"=500)
#'
simulate_sequences_HIV_V3 <- function(params) {  # nolint: object_name_linter
    # This gives an example "workflow" function for this package.
    # Simple integration tests should use similar workflows to insert changes.

    # First, determine the rate function for offspring generation
    # This obtains the function used in [Graw et al. 2012], [Kupperman et al. 2022]
    biphasic_rate_function <- SEEPS::get_biphasic_HIV_rate(
        params = params[["rate_function_parameters"]])

    # Forward pass to generate transmission history
    simulator_result <- SEEPS::gen_transmission_history_exponential_constant(
        minimum_population = params[["minimum_population"]],
        offspring_rate_fn = biphasic_rate_function,
        total_steps = params[["total_steps_after_exp_phase"]],
        maximum_population_target = params[["maximum_population_target"]],
        spike_root = FALSE)

    # Simulate contact tracing on the contact network to identify a sample
    target_sample <- SEEPS::contact_traced_uniform_ids(
        active = simulator_result[["active"]],
        parents = simulator_result[["parents"]],
        minimum_sample_size = params[["minimum_population"]],
        # probability of uncovering each contact
        p = params[["contact_tracing_discovery_probability"]])

    # Get a geneology compatable with biophybreak (bpp) layout.
    res <- SEEPS::reduce_transmission_history_bpb(
        samples = target_sample[["samples"]],
        parents = simulator_result$parents,
        current_step = simulator_result$t_end)

    # Provide values of a and b for the within-host diversity models.
    phylogeny <- SEEPS::geneology_to_phylogeny_bpb(
        transmission_history = res$parents,
        infection_times = res$transmission_times,
        sample_times = res$sample_times,
        a = params[["a"]], b = params[["b"]],
        leaf_sample_ids = res$transformed_sample_indices)


    # Convert the time signals into a # of mutations per site
    phylogeny <- SEEPS::stochastify_transmission_history(
            transmission_history = phylogeny$phylogeny,
            # Expect a rate in per nt per year
            # input distances are in months
            rate = params[["mutation_rate"]] / 12)

    # Get a the rate model for V3
    rate_model <- SEEPS::get_V3_rate_model()

    # Get the HIV1 reference sequence for V3
    V3_sequence <- SEEPS::lookup_sequence_by_name(organism_name = "HIV1",  # nolint: object_name_linter
                                                  region_name = "V3")

    # Call Seq-Gen to generate sequences from the provided reference sequence.
    sequences <- SEEPS::generate_sequences(
        phylogeny = phylogeny$geneology,
        branch_rate = params[["mutation_rate"]] / 12,
        root_sequence = V3_sequence,
        rate_model = rate_model,
        rate_per_nt = FALSE)  # Distances in phylogeny are already in per nt
    # The return is a string that we can write to a fasta file.

    # Instead, we'll convert it to a dataframe
    sequences <- SEEPS::fasta_string_to_dataframe(sequences)

    # Now build a pairwise distance matrix from the dataframe
    matrix <- SEEPS::build_distance_matrix_from_df(sequences)

    reduced_mat_data <- SEEPS::reduce_large_matrix(
        oversampled_matrix = matrix,
        subsample_size = params[["minimum_population"]],
        spike_root = FALSE)
    # Drop the sequences we didn't keep
    sequences <- sequences[reduced_mat_data$keep_indices, ]
    return(list("sequences" = sequences, "params" = params,
                "matrix" = reduced_mat_data$matrix))
}
