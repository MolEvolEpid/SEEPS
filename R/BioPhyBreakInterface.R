#' Generate a phylogeny from a transmission history using BioPhyBreak's coalescent simulator
#'
#' Run a backwards simulation using the transmission history to generate a phylogeny.
#' See `geneology_to_phylogeny_bpb` for details on how the input should be structured.
#' This function assumes that all coalescent events occur at or after the initial infection,
#' that there is a single introduction and ancestral sequence.
#'
#'
#' @param transmission_history A transmission history matrix
#' @param infection_times A vector of infection times
#' @seealso reduce_transmission_history_bpb
#' @seealso generate_sequences
#' @importFrom ape read.tree
#' @export
geneology_to_phylogeny_bpb <- function(transmission_history,
        infection_times, leaf_sample_ids,
        sample_times, a = 5, b = 5, make_plot = FALSE) {

    n <- length(transmission_history)
    labels <- 1:n
    tt <- data.frame(donor = as.vector(transmission_history),
            names = labels,
            # time_infection
            t_inf = as.vector(infection_times) / 12,  # Convert from months to years
            #time_sampling
            t_sam = as.vector(sample_times) / 12)  # Convert from months to years
    # This calls the bpb code
    tree_and_tips <- sim.coal.tree.m(tt = tt, a = rep(a, n), b = rep(b, n),
                                     tree_name = 1,
                                     leaf_raw_ids = leaf_sample_ids)

    # Need to pass in a seed argument here
    # Need to check here that ape is installed first. Else, warn and continue.
    # todo: check that ape is available before converting
    tree <- ape::read.tree(text = tree_and_tips[["tree_newick"]])
    if (make_plot) plot(tree)

    return(list("names_newick" = tree_and_tips[["names_newick"]],
                "phylogeny" = tree_and_tips[["phylogeny"]],
                "newick_string" = tree_and_tips[["names_newick"]], "tree" = tree))
}
