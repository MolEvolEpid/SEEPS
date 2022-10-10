gen_transmission_history_balanced_tree <- function(population_size, #nolint: object_length_linter
                              spike_root = FALSE) {
    # Simulate a balanced tree wih `population_size` leaves.

    depth <- floor(log2(population_size + 1))
    num_nodes <- 2^(depth)  # internal nodes, only leaf nodes if population_size is a power of 2
    if (num_nodes == population_size) {
        # If population_size is a power of 2, we don't need to add any extra nodes
        K <- 0
    } else {
        K <- population_size  # Number of extra leaves we need
    }

    num_nodes <- 2 * num_nodes + K
    parents <- matrix(0, num_nodes, 2)  # DAG for offspring

    # Declare counters and pointers
    flag <- TRUE
    parents[2,] <- c(1,1)  # Setup the index case (1) at time 0 from case 1
     # Setup the index case
    birth_step <- c(1)
    # calculate the depth of a binary tree needed for population_size leaves
    # We're going to do a 2-index approch.
    parent_vector <- c(1)
    offspring_vector <- c(2,3)
    pointer <- 2
    current_depth <- 1
    flag <- TRUE
    while (flag) {
        # Shuffle the pointers for the index vectors
        parent_vector <- offspring_vector
        # Build out a new offspring vector
        offspring_vector <- (2 ** current_depth):((2 ** (current_depth + 1)) - 1)
        # Enforce the cap
        vector_length <- length(offspring_vector)
        if (vector_length >= population_size) {  # We have extra population
            offspring_vector <- offspring_vector[1:(population_size)]
            # This pass through will fill us up, so don't go again
            flag <- FALSE
        }

        addition_length <- length(offspring_vector)
        # For each parent, add the appropriate number of offspring
        if (addition_length > 1) {  # If we have only one offspring to add,
        # we don't need to do this skip ahead to the edge case handler

            for (i in 1:(addition_length %/% 2)) {  # Loop over all parents

                index <- pointer + 1
                parents[index, 1] <- parent_vector[i]  # Record the parent
                parents[index, 2] <- current_depth  # Record the birth time
                # Try to add in another
                index <- pointer + 2
                parents[index, 1] <- parent_vector[i]  # Record the parent
                parents[index, 2] <- current_depth  # Record the birth time
                # Update the global counter
                pointer <- pointer + 2
            }
        }

        # If we have an odd number of offspring, add one more.
        # Only an issue with the final layer of leaves
        if (addition_length %% 2 == 1) {
            index <- pointer + 1
            parents[index, 1] <- parent_vector[(addition_length %/% 2) + 1]  # Record the parent
            parents[index, 2] <- current_depth  # Record the birth time
            # Update the global counter
            pointer <- pointer + 1
        }
        current_depth <- current_depth + 1
    }
    active <- offspring_vector  # The newly leaf nodes are the active nodes
    # Add on any parent nodes we didn't attach offspring to
    if (2 * length(parent_vector) < length(offspring_vector)) {
        # We have extra parents, get the tail we didn't assign offspring to
        extra_parents <- parent_vector[(length(offspring_vector) + 1):length(parent_vector)]
        active <- c(active, extra_parents)  # Not worth optimizing this concatenate out
    }
    return(parents)
}

# Should have
# 0 0
# 1 1
# 2 2
# 3 2
# 4 3
# 4 3
# 5 3
# 6 3
# 7 4
# 7 4
# 8 4
# 8 4
# 9 4
# 9 4
# 10 4
# 10 4
# 11 5
# 11 5
# 12 5
# 12 5
# 13 5
# 13 5
# 14 5
# 14 5
# 15 5
# 15 5
# 16 5
# 16 5
# 17 5
# 17 5
# 18 5
# 18 5
