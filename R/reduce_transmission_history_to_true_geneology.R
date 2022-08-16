#' Reduce simulation output to transmission history for a subset.
#'
#' @return
#' @export
#'
#' @importFrom stats rpois
#' @examples
reduce_transmission_history <- function(samples, parents,
                                        current_step, spike_root = FALSE) {
    observation_size <- length(samples)
    geneology <- matrix(0, 2 * observation_size, 5)
    geneology[1:(2 * observation_size - 1), 1] <- c(1:(2 * observation_size - 1))

    # Taking into account of the sampling time
    geneology[1:observation_size, 3] <- current_step + rpois(observation_size, 6)
    if (spike_root) {
        # Sample date for the initial infection is 0, no tip adjustment
        geneology[observation_size, 3] <- 0 # starts at the origin
    }
    nNodes <- observation_size
    IDs <- 1:observation_size
    pp <- observation_size + 1


    while (nNodes > 1) {
        # Get the newest infection
        ind <- which.max(parents[samples, 2])
        # Get the parent of the newest infection
        parent <- parents[samples[ind], 1]
        # Is the parent is in the sample?
        coall <- which(samples == parent)
        if (length(coall) == 0) { # no
            # Update the parent
            samples[ind] <- parent
        } else { # We have a coalescent event
            # Store the sample time for the parent
            geneology[pp, 3] <- parents[samples[ind], 2]
            # Write the merge/coal/join information into the array
            geneology[IDs[ind], 2] <- pp
            # Record it for both parents
            geneology[IDs[coall], 2] <- pp
            if (!is.na(IDs[coall])) { # If we have sequences left, write the join time
                # Only happens when we try to backprop from the first infection
                if (geneology[IDs[coall], 3] < parents[samples[ind], 2]) {
                    geneology[IDs[coall], 3] <- parents[samples[ind], 2]
                }
            }
            IDs[coall] <- pp
            samples <- samples[-ind]
            IDs <- IDs[-ind]
            pp <- pp + 1
            nNodes <- nNodes - 1
        }
    }
    inds <- 1:(2 * observation_size - 1)
    geneology[inds, 3] <- geneology[inds, 3] - geneology[(2 * observation_size - 1), 3]
    # Subtract off from the root
    geneology[inds, 3] <- max(geneology[inds, 3]) - geneology[inds, 3]
    # Happens to place the root back at 0 if we spiked in the origin


    inds <- 1:(2 * observation_size - 2)
    # Convert absolute times to local (relative) differences.
    geneology[inds, 4] <- geneology[geneology[inds, 2], 3] - geneology[inds, 3]
    geneology[inds, 5] <- geneology[inds, 4] / (sum(geneology[inds, 4]))
    return(list("geneology" = geneology))
}

###########################################################
#
# The above algorithm recovers the tree only for a
# subset of active individuals. We need the full network
# for the subset for within-host dynamics.
# This function is designed for biophylobreak compatability.
#
###########################################################


#'
#' @export
reduce_transmission_history_bpb <- function( # nolint:object_length_linter
                                            samples, parents, current_step) {
    leaves <- samples

    # Infections increase monotonically. Adding new sequences gives increased index.
    # Largest index in sample is largest number of individuals we need to consider.
    # pre-allocate for speed, these could be large lists and we want to avoid a copy
    parents_tree <- rep(-1, length = max(samples))
    offspring_times <- rep(-1, length = max(samples))
    # We'll use -1 as a status code for "not in the history. We'll need this for cleanup later.


    for (sample in samples) {
        exit_flag <- TRUE
        while (exit_flag) {
            # Get which infection we are going to back-prop  to
            parent_index <- parents[sample, 1] #
            # Record the parent for the newly found infection as the cause of the infection
            parents_tree[sample] <- parent_index
            # Record the birth time
            if (offspring_times[sample] == -1) { # We don't know the parent, detect coalesence
                offspring_times[sample] <- parents[sample, 2] # col 2 is for sample time
                sample <- parent_index
            } else { # We hit a coalecent event, we can stop now
                exit_flag <- FALSE
            }
            if (sample == 0) { # We hit the root of the outbreak.
                exit_flag <- FALSE
            }
        }
    }

    # TODO split this function here

    # Now "trim" the lists and compute observation times.
    trimmed_parents_tree <- parents_tree
    trimmed_offspring_times <- offspring_times
    leaf_locations_mask <- parents_tree == -2 # Create vector of FALSE
    leaf_locations_mask[leaves] <- TRUE # Record locations
    # print(leaf_locations_mask)

    exit_flag <- TRUE
    index <- 1
    while (exit_flag) {
        value <- trimmed_parents_tree[index]
        if (value == -1) {
            # decrement the tail
            mask <- trimmed_parents_tree >= index
            # We cannot have values larger than `position` before the position
            # no need to filter by position, but that might give a small speedup for large vectors
            trimmed_parents_tree[mask] <- trimmed_parents_tree[mask] - 1
            # Now shift the list and drop the -1
            trimmed_parents_tree <- trimmed_parents_tree[-index]
            trimmed_offspring_times <- trimmed_offspring_times[-index]
            # Also update the vector of leaf information
            leaf_locations_mask <- leaf_locations_mask[-index]
        }
        if (index == length(trimmed_parents_tree)) {
            # At the end. Exit now
            exit_flag <- FALSE
        } else {
            if (value != -1) {
                # If we modified the list, the position has a new value
                # Only increase the index if we didn't change the list
                index <- index + 1
            }
        }
    }

    # Function 3

    # Now set the sample times
    # For leaves, set the sample time to current_step
    # for internal nodes, we need thier last offspring

    internal_nodes_mask <- as.logical(1 - leaf_locations_mask)
    end_times <- rep(-1, length(trimmed_parents_tree))
    end_times[leaf_locations_mask] <- current_step
    internal_nodes <- which(internal_nodes_mask)

    for (node in internal_nodes) {
        # Get the offspring
        offspring_indices <- which(trimmed_parents_tree == node)
        # Get their generation times
        times <- trimmed_offspring_times[offspring_indices]
        # Record the largest generation time
        end <- min(max(times) + 1e-2, current_step) # When to stop simulating
        # Add in an epsilon for numerical stability,
        end_times[node] <- end
    }

    return(list(
        "parents" = trimmed_parents_tree,
        "transmission_times" = trimmed_offspring_times,
        # Record when samples were taken
        "sample_times" = end_times,
        # Record how many samples  there are. Either 0 (internal) or 1(leaf)
        "samples_available" = 1 * leaf_locations_mask
    ))
}
