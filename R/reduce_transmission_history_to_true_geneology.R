#' Reduce simulation output to transmission history for a subset.
#'
#' Reduce a large simulation output to a smaller transmission history for a subset
#' by tracing back the ancestory of each individual in the sample. If `spike_root`
#' is `TRUE`, then the root of the tree is included in the geneology.
#' We call this a geneology, rather than a phylogeny as it assumes that the
#' transmission history is the phylogeny and no within-host diversity occurs.
#'
#' To include within host diversity, use `reduce_transmission_history_bpb` to
#' extract a subst of the transmission history that we need for the phylogeny,
#' and see  `geneology_to_phylogeny_bpb` to simulate within-host diversity
#' and recover a true phylogeny.
#'
#' @param samples A vector of individuals (integers) to include in the sample.
#' @param parents A matrix of parental individuals that encodes the transmission
#'             history and sample times.
#' @param current_step The current (absolute) time step in the simulation.
#' @param spike_root A boolean indicating whether the geneology should
#'            include the root of the outbreak or not. Default is `FALSE`.
#'            This should be specified even if the founding infection is sampled,
#'           as the root of the outbreak will have evolved since the founding event.
#'
#' @seealso reduce_transmission_history_bpb
#' @return A list with 1 element: "geneology" a matrix of transmission history
#'   that encodes an evolutionary tree.
#' @export
#'
#' @importFrom stats rpois
reduce_transmission_history_mt <- function(samples, parents,
                                        current_step, spike_root = FALSE) {
    samples_all <- unlist(samples)
    sample_times <- unlist(
        sapply(
            1:length(samples),
            function(i) {
                rep(current_step[[i]], length(samples[[i]]))
            }
        )
    )
    samples <- samples_all  # Over-write the samples, now that we know when
    # they were sampled

    observation_size <- length(samples)
    geneology <- matrix(0, 2 * observation_size, 7)
    # Columns:
    # 1. Local Individual ID (1:n, 0)
    # 2. Parent ID (1:n, 0)
    # 3. Time of infection (non-negative integers)
    # 4. Branch length (non-negative integers)
    # 5. Distances (normalized here, but will be filled later via sampling)
    # 6. Leaf status (bool)
    # 7. Absolute index - leaves only. Number from samples

    geneology[1:(2 * observation_size - 1), 1] <- c(1:(2 * observation_size - 1))

    # Taking into account of the sampling time
    geneology[1:observation_size, 3] <- sample_times  # + rpois(observation_size, 6)
    geneology[1:observation_size, 6] <- 1  # These are leaves
    if (spike_root) {
        # Sample date for the initial infection is 0, no tip adjustment
        geneology[observation_size, 3] <- 0 # starts at the origin
    }
    # Set the initial inds for each node
    original_samples <- samples
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

            # Store the sample ID in column 7
            # We need to not have written a value, and be a leaf
            if (geneology[IDs[ind], 7] == 0 && geneology[IDs[ind], 6] == 1) {
                geneology[IDs[ind], 7] <- original_samples[ind]
            }
            if (geneology[IDs[coall], 7] == 0 && geneology[IDs[coall], 6] == 1) {
                geneology[IDs[coall], 7] <- original_samples[coall]
            }


            IDs[coall] <- pp
            samples <- samples[-ind]
            original_samples <- original_samples[-ind]  # remove from copy
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

#' Reduce simulation output to transmission history for a subset.
#'
#' Reduce a large simulation output to a smaller transmission history for a subset
#' by tracing back the ancestory of each individual in the sample. If `spike_root`
#' is `TRUE`, then the root of the tree is included in the geneology.
#' We call this a geneology, rather than a phylogeny as it assumes that the
#' transmission history is the phylogeny and no within-host diversity occurs.
#'
#' To include within host diversity, use `reduce_transmission_history_bpb` to
#' extract a subst of the transmission history that we need for the phylogeny,
#' and see  `geneology_to_phylogeny_bpb` to simulate within-host diversity
#' and recover a true phylogeny.
#'
#' @param samples A vector of individuals (integers) to include in the sample.
#' @param parents A matrix of parental individuals that encodes the transmission
#'             history and sample times.
#' @param current_step The current (absolute) time step in the simulation.
#' @param spike_root A boolean indicating whether the geneology should
#'            include the root of the outbreak or not. Default is `FALSE`.
#'            This should be specified even if the founding infection is sampled,
#'           as the root of the outbreak will have evolved since the founding event.
#'
#' @seealso reduce_transmission_history_bpb
#' @return A list with 1 element: "geneology" a matrix of transmission history
#'   that encodes an evolutionary tree.
#' @export
#'
#' @importFrom stats rpois
reduce_transmission_history <- function(samples, parents,
                                        current_step, spike_root = FALSE) {
    observation_size <- length(samples)
    geneology <- matrix(0, 2 * observation_size, 7)
    # Columns:
    # 1. Local Individual ID (1:n, 0)
    # 2. Parent ID (1:n, 0)
    # 3. Time of infection (non-negative integers)
    # 4. Branch length (non-negative integers)
    # 5. Distances (normalized here, but will be filled later via sampling)
    # 6. Leaf status (bool)
    # 7. Absolute index - leaves only. Number from samples

    geneology[1:(2 * observation_size - 1), 1] <- c(1:(2 * observation_size - 1))

    # Taking into account of the sampling time
    geneology[1:observation_size, 3] <- current_step  # + rpois(observation_size, 6)
    geneology[1:observation_size, 6] <- 1  # These are leaves
    if (spike_root) {
        # Sample date for the initial infection is 0, no tip adjustment
        geneology[observation_size, 3] <- 0 # starts at the origin
    }
    # Set the initial inds for each node
    original_samples <- samples
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

            # Store the sample ID in column 7
            # We need to not have written a value, and be a leaf
            if (geneology[IDs[ind], 7] == 0 && geneology[IDs[ind], 6] == 1) {
                geneology[IDs[ind], 7] <- original_samples[ind]
            }
            if (geneology[IDs[coall], 7] == 0 && geneology[IDs[coall], 6] == 1) {
                geneology[IDs[coall], 7] <- original_samples[coall]
            }


            IDs[coall] <- pp
            samples <- samples[-ind]
            original_samples <- original_samples[-ind]  # remove from copy
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


#' Reduce simulation output to transmission history for a subset to include
#'   within host diversity.
#'
#' For a detailed explenation of inputs, see `reduce_transmission_history`, which
#' is intended to reconstruct back only until the most recent common ancestor
#'  of the sample, and return a tree.
#'
#' To include within host diversity, use `reduce_transmission_history_bpb` to
#' extract a subset of the transmission history that we need for the phylogeny,
#' and see `geneology_to_phylogeny_bpb` to simulate within-host diversity
#' and recover a true phylogeny.
#'
#' @param samples A vector of individuals (integers) to include in the sample.
#' @param parents A matrix of parental individuals that encodes the transmission
#'             history and sample times.
#' @param current_step The current (absolute) time step in the simulation.
#' @param spike_root A boolean indicating whether the geneology should
#'            include the root of the outbreak or not. Default is `FALSE`.
#'            This should be specified even if the founding infection is sampled,
#'           as the root of the outbreak will have evolved since the founding event.
#'
#' @seealso reduce_transmission_history_bpb
#' @return A list with 4 elements:
#' `parents` A vector of parents of each infection in the sample until the root
#' `times` A vector of times of sampling times in the tree. Sample times for
#' internal nodes are after the last offspring generation time needed to
#' reconstruct the sample.
#' `transmission_times` A vector of transmission times of each infection in the tree.
#' `samples_available` A boolean vector (mask) of which samples are leaves in the tree.
#'    Used by the coalescent simulation to know which individuals should be assigned
#'    detected sequences.
#'
#' @export
#'
#' @importFrom stats rpois
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
            } else { # We hit a coallecent event, we can stop now
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
    leaf_indices <- as.integer(parents_tree == -2)  # Create vector of 0
    leaf_indices[leaves] <- samples

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
            # Also update the vectors of leaf information
            leaf_locations_mask <- leaf_locations_mask[-index]
            leaf_indices <- leaf_indices[-index]
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
        end <- min(max(times) + 1e-3, current_step) # When to stop simulating
        # Add in an epsilon for numerical stability,
        end_times[node] <- end
    }

    return(list(
        "parents" = trimmed_parents_tree,
        "transmission_times" = trimmed_offspring_times,
        # Record when samples were taken
        "sample_times" = end_times,
        # Record how many samples  there are. Either 0 (internal) or 1(leaf)
        "samples_available" = 1 * leaf_locations_mask,
        # The index of the sample in the original list
        "transformed_sample_indices" = leaf_indices
    ))
}

############################

# A second version of the above function that allows for samples to be taken
# at different times. This is useful for simulating coalescent trees with
# multiple samples taken.

############################

#' @export
reduce_transmission_history_bpb2 <- function( # nolint:object_length_linter
                                            samples, parents, current_step) {
    # samples is a list of samples taken at different times
    # sample_times is a vector of times at which samples were taken
    # samples[[1]] is a vector of samples taken at current_step[[1]]
    leaves <- unlist(samples)
    samples_list <- samples  # Store the list of samples
    samples <- unlist(samples)  # Flatten the list of samples


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
            } else { # We hit a coallecent event, we can stop now
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
    # Make a mask of the sample times, so we can account for the trimming
    sample_time_vector <- rep(-1, length = length(parents_tree))  # -1 as status code
    for (i in seq_along(samples_list)) {
        # Store the index of the sample in the original list
        sample_time_vector[samples_list[[i]]] <- i
        # Don't convert to sample times yet, we need to trim the list first
    }
    leaf_indices <- as.integer(parents_tree == -2)  # Create vector of 0
    leaf_indices[leaves] <- samples


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
            # Also update the vectors of leaf information
            leaf_locations_mask <- leaf_locations_mask[-index]
            leaf_indices <- leaf_indices[-index]
            sample_time_vector <- sample_time_vector[-index]
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
    sample_batch <- rep(-1, length(trimmed_parents_tree))
    # end_times[leaf_locations_mask] <- current_step
    # loop over the pairs of samples and set the end_times
    for (i in seq_along(samples_list)) {
        mask <- sample_time_vector == i
        end_times[mask] <- current_step[[i]]
        sample_batch[mask] <- i
    }
    internal_nodes <- which(internal_nodes_mask)

    for (node in internal_nodes) {
        # Get the offspring
        offspring_indices <- which(trimmed_parents_tree == node)
        # Get their generation times
        times <- trimmed_offspring_times[offspring_indices]
        # Record the largest generation time
        end <- min(max(times) + 1e-3, Inf) # When to stop simulating
        # Add in an epsilon for numerical stability,
        end_times[node] <- end
    }


    return(list(
        "parents" = trimmed_parents_tree,
        "transmission_times" = trimmed_offspring_times,
        # Record when samples were taken
        "sample_times" = end_times,
        # Record how many samples  there are. Either 0 (internal) or 1(leaf)
        "samples_available" = 1 * leaf_locations_mask,
        # Denote which time point the samples came from.
        "sample_group_vectors" = sample_batch,
        # The index of the sample in the original list
        "transformed_sample_indices" = leaf_indices
    ))
}
