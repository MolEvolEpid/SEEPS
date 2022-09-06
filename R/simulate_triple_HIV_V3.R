#' Sample distance matrices for all three major simulation paradigms:
#' the transmission history, the phylogeny, and derived from simulated sequences (TN93)
#'
#' This simulation returns three distance matrices for the same set of individuals.
#' The first matrix (matrix_trans) is derived from the transmission history.
#' The second matrix (matrix_phylo) is derived from the phylogeny by simulating a
#' coalescent process on the transmission history.
#' The third matrix (matrix_seq) is derived from the sequences by computing the
#' pairwise distances using the TN93 model.
#'
#' @return list with keys "matrix" and "input_params"
#' @export
#'
#' @examples
#' parameters <- list("rate_function_parameters" = list("R0"=5), "a"=5, "b"=5,
#'                     "mutation_rate" = 0.0067,  # units are per nt
#'                     "contact_tracing_discovery_probability" = 0.9,
#'                     "minimum_population"=15, "maximum_population_target"=500)
#'
simulate_all_paradigms_HIV_V3 <- function(params) {  # nolint: object_name_linter
    # This gives an example "workflow" function for this package.
    # Simple integration tests should use similar workflows to insert changes.

    ############ Simulate transmission history and contact tracing ############
    # This obtains the function used in [Graw et al. 2012], [Kupperman et al. 2022]

    # Get the HIV1 reference sequence for V3
    V3_sequence <- SEEPS::lookup_sequence_by_name(organism_name = "HIV1",  # nolint: object_name_linter
                                                  region_name = "V3")
    seq_len <- length(V3_sequence)

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
    target_sample <- SEEPS::contact_traced_uniform_restarts_ids(
        active = simulator_result[["active"]],
        parents = simulator_result[["parents"]],
        minimum_sample_size = params[["minimum_population"]],
        # probability of uncovering each contact
        p = params[["contact_tracing_discovery_probability"]])


    ###################### Transmission history to matrix ######################


    # Get a geneology from the transmission history for the first matrix
    geneology_transmission <- SEEPS::reduce_transmission_history(
        samples = target_sample[["samples"]],
        parents = simulator_result$parents,
        current_step = simulator_result$t_end)

    # Convert time signals to # of mutations using a rate
    geneology <- SEEPS::stochastify_transmission_history(
        transmission_history = geneology_transmission$geneology,
        # Rate here is in per nt per month per sequence
        rate = params[["mutation_rate"]] / 12 * seq_len)

    # convert the geneology into a distance matrix
    matrix_trans <- SEEPS::geneology_to_distance_matrix(
        geneology = geneology$geneology,
        spike_root = FALSE)

    # Take a subsample of closest neighbors
    reduced_matrix_trans <- SEEPS::reduce_large_matrix(
        oversampled_matrix = matrix_trans,
        subsample_size = params[["minimum_population"]],
        spike_root = FALSE)

    ###### Simulate phylogeny from transmission history for second matrix ######

    # Get a geneology compatable with biophybreak (bpb) layout.
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
            rate = params[["mutation_rate"]] / 12 * seq_len)

    # Only reconstruct the nodes corresponding to the sampled tips
    matrix_phylo <- SEEPS::geneology_to_distance_matrix(
        geneology = phylogeny$geneology,
        spike_root = FALSE)

    # Take a subsample of closest neighbors
    reduced_matrix_phylo <- SEEPS::reduce_large_matrix(
        oversampled_matrix = matrix_phylo,
        subsample_size = params[["minimum_population"]],
        spike_root = FALSE,
        index_id = reduced_matrix_trans$sampled_index,
        sort_order = reduced_matrix_trans$sort_indices)

    ############ Simulate sequences from phylogeny for third matrix ############

    # Get a the rate model for V3
    rate_model <- SEEPS::get_V3_rate_model()


    # Call Seq-Gen to generate sequences from the provided reference sequence.
    sequences <- SEEPS::generate_sequences(
        phylogeny = phylogeny$geneology,
        branch_rate = params[["mutation_rate"]] / 12,  # Convert time distances to per subs per site
        root_sequence = V3_sequence,
        rate_model = rate_model,
        rate_per_nt = FALSE)  # Distances in phylogeny are already in per nt
    # The return is a string that we can write to a fasta file.

    # Instead, we'll convert it to a dataframe
    sequences <- SEEPS::fasta_string_to_dataframe(sequences)

    # Now build a pairwise distance matrix from the dataframe
    matrix_seq <- SEEPS::build_distance_matrix_from_df(sequences, model = "TN93")

    ################# Reduce matrices through nearest neighbor #################

    # The sequence based distance matrix is probably the most accurate.
    reduced_mat_data <- SEEPS::reduce_large_matrix(
        oversampled_matrix = matrix_seq,
        subsample_size = params[["minimum_population"]],
        spike_root = FALSE,
        index_id = reduced_matrix_trans$sampled_index,
        sort_order = reduced_matrix_trans$sort_indices)
    # Drop the sequences we didn't keep
    sequences <- sequences[reduced_mat_data$keep_indices, ]
    # Drop the rows/cols of the trans and phylo matrices we didn't keep



    return(list("sequences" = sequences, "params" = params,
                "matrix_seqs" = reduced_mat_data$matrix,
                "matrix_trans" = reduced_matrix_trans$matrix,
                "matrix_phylo" = reduced_matrix_phylo$matrix))
}
