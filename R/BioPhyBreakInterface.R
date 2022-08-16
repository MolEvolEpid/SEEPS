# Public facing utility functions for interfacing with included BioPhyBreak functionss
# call BioPhyBreak
#' @importFrom ape read.tree
#' @export
geneology_to_phylogeny_bpb <- function(transmission_history,
        infection_times,
        sample_times, a=5, b=5, make_plot=FALSE) {

    n <- length(transmission_history)
    labels <- 1:n
    tt <- data.frame(donor = as.vector(transmission_history),
            names = labels,
            # time_infection
            t_inf = as.vector(infection_times) / 12,  # Convert from months to years
            #time_sampling
            t_sam = as.vector(sample_times) / 12)  # Convert from months to years
    # This calls the bpb code
    tree_and_tips <- sim.coal.tree.m(tt=tt, a = rep(a, n), b = rep(b, n), tree_name = 1)

    # Need to pass in a seed argument here
    # Need to check here that ape is installed first. Else, warn and continue.
    # todo: check that ape is available before converting
    tree <- ape::read.tree(text = tree_and_tips[["tree_newick"]])
    if(make_plot) plot(tree)
    # print(tree_and_tips[["phylogeny"]])
    return(list("names_newick" = tree_and_tips[["names_newick"]],
                "phylogeny" = tree_and_tips[["phylogeny"]],
                "newick_string" = tree_and_tips[["names_newick"]], "tree" = tree))
}
