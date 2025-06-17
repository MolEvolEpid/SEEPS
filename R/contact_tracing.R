################################################################################
#
# Contact tracing functionality for generating samples correlation
#
# Users provide a discovery rate function that describes the information obtained
# from a contact tracing interview. The function is used to accept or reject
# new detections, new detections are iteratively contact traced until a termination
# condition is reached, or we run out of new detections to trace.
#
################################################################################

#' Obtain a sample using a iterative contact tracing with a uniform discovery rate
#'
#' Perform iterative contact tracing on a simulated contact network.
#' Randomness within contact tracing comes from the probability of discovering
#' a new infection. This function assumes the infection rate is uniform across
#' all possible contacts and individuals.
#'
#' Contact tracing is terminated when 2 * min_sample_size active nodes are discovered,
#' or when there are no more detections to trace. At least `min_sample_size`
#' individuals will always be returned.
#'
#' To determine the initial detection, we loop over a the list of active nodes
#' until we complete a successful contact tracing.
#'
#' All contact tracing algorithms store both the ids of all discovered individuals
#' and the ids of discovered active individuals.
#' @param active A vector of active individuals
#' @param parents A matrix encoding the transmission history
#' @param minimum_sample_size The minimum number of individuals to form a sample
#' @param p The probability of discovering each contact during tracing
#' @param max_attempts Maximum number of attempts to perform to obtain a sample.
#'   If no sample is found after this number of attempts, `FALSE` is returned
#' @return A list with three fields: "status", "samples", and "found"
#' @export
contact_traced_uniform_ids <- function(active, parents, minimum_sample_size, p,
                                       max_attempts = 1) {
    discovery_function <- uniform_discovery_factory(p = p)
    termination_function <- sufficient_data_data_factory(
        minimum_size = 2 * minimum_sample_size
    )
    # Call the contact tracing engine
    initial_detections <- sample(active, length(active), replace = FALSE)
    result <- contact_tracing_engine(
        detected_id = initial_detections,
        active = active,
        parents = parents,
        minimum_size = minimum_sample_size,
        discovery_function = discovery_function,
        termination_function = termination_function,
        max_attempts = max_attempts
    )
    # Check the output before we return
    if (result[["success"]]) {
        return(result)
    } else {
        warning("Failed to find a sample after", max_attempts, "attempts")
        # print(result)
        return(FALSE)
    }
}

# TODO: Refactor this with the above function to a single function with dispatcher
#' Obtain a sample using a iterative contact tracing with a uniform discovery rate
#' If a minimum number is not found, the algorithm is restarted with a new initial detection
#'
#' Perform iterative contact tracing on a simulated contact network.
#' Randomness within contact tracing comes from the probability of discovering
#' a new infection. This function assumes the infection rate is uniform across
#' all possible contacts and individuals.
#'
#' Contact tracing is terminated when 2 * min_sample_size active nodes are discovered,
#' or when there are no more detections to trace. At least `min_sample_size`
#' individuals will always be returned.
#'
#' To determine the initial detection, we loop over a the list of active nodes
#' until we complete a successful contact tracing.
#'
#' All contact tracing algorithms store both the ids of all discovered individuals
#' and the ids of discovered active individuals.
#' @param active A vector of active individuals
#' @param parents A matrix encoding the transmission history
#' @param minimum_sample_size The minimum number of individuals to form a sample
#' @param p The probability of discovering each contact during tracing
#' @param initial_detections An optional list of initial detections to start the algorithm with.
#'   If not specified, the algorithm will randomly select a set of initial detections.
#' @return A list with four fields: "status", "samples", "success", and "found"
#' if the algorithm fails to find a sample, FALSE is returned instead
#' @export
contact_traced_uniform_restarts_ids <- function(active, parents, # nolint: object_name_linter
                                                minimum_sample_size, p,
                                                initial_detections = NULL) {
    discovery_function <- uniform_discovery_factory(p = p)
    termination_function <- sufficient_data_data_factory(
        minimum_size = minimum_sample_size
    )
    # Call the contact tracing engine
    if (is.null(initial_detections)) {
        initial_detections <- sample(active, length(active), replace = FALSE)
    } else {  # check the input
        # check we got a numeric vector of ids
        if (!is.numeric(initial_detections)) {
            stop("Initial detections must be a numeric vector of ids")
        }

        # check that the ids are in the active list
        if (any(!initial_detections %in% active)) {
            stop("Initial detections must be a subset of the active nodes")
        }
    }
    results <- list()
    group_id_counter <- 1 # Track which group each individual was sampled from
    group_ids <- c()
    attempt_counter <- 1
    samples <- c()
    target_size <- minimum_sample_size # We modify this every time we restart
    repeat { # do-while pattern
        # Remove any discovered nodes from initial_detections
        initial_detections <- initial_detections[!(initial_detections %in% samples)]
        result <- contact_tracing_engine(
            detected_id = initial_detections,
            active = active,
            parents = parents,
            minimum_size = target_size,
            discovery_function = discovery_function,
            termination_function = termination_function,
            max_attempts = 1
        ) # Do not try and repeat nodes
        # Check the output before we return
        # A bit of a hack to get the group ids
        result[["group_ids"]] <- group_id_counter + 0 * result[["samples"]]
        group_id_counter <- group_id_counter + 1
        results[[attempt_counter]] <- result
        attempt_counter <- attempt_counter + 1
        # concatenate the results
        sampleStruct <- clean_sample_structure(results)
        # Unpack
        samples <- sampleStruct$samples
        group_ids <- sampleStruct$group_ids

        # Now update the target size by subtracting off the number of samples we have
        target_size <- target_size - length(samples)
        # while
        if (length(samples) >= minimum_sample_size) {
            # We have enough data, take the first set we discovered
            samples <- samples[1:minimum_sample_size]
            group_ids <- group_ids[1:minimum_sample_size]
            found <- unique(unlist(sapply(results, function(x) x[["found"]])))
            break
        }
    }

    # Now update the list we're going to return
    result[["success"]] <- TRUE
    result[["samples"]] <- samples
    result[["group_ids"]] <- group_ids
    result[["found"]] <- found
    result[["status"]] <- paste(
        "Found enough data with",
        attempt_counter - 1, "attempt(s)"
    )
    if (result[["success"]]) {
        return(result)
    } else {
        return(FALSE)
    }
}
#' Clean the sample structure
#'
#' Given a list of results from contact tracing, return a list of unique samples
#' along with the group id for each of them in two lists.
#'
#' @param results A list of results from contact tracing. See
#' `SEEPS::contact_tracing_engine` for more details.
#'
#' @export
clean_sample_structure <- function(results) {
    # Return a list of unqiue samples along with the group id for each of them in two lists
    # Create
    all_samples <- unlist(lapply(results, function(x) x[["samples"]]))
    all_group_ids <- unlist(lapply(results, function(x) x[["group_ids"]]))
    df_merged <- data.frame("samples" = I(all_samples), "group_ids" = I(all_group_ids))
    df_merged <- df_merged[!duplicated(df_merged$samples), ]

    # Return the two numeric vectors from the columns of df_merged by key
    return(list(
        "samples" = as.numeric(unlist(df_merged$samples)),
        "group_ids" = as.numeric(unlist(df_merged$group_ids))
    ))
}

################################################################################
#
# Core algorithms
#
################################################################################

#' Perform contact tracing on a simulated contact network
#'
#' The engine for simulating contact tracing. Provide network information
#' (`active`, `parents`), a list or vector of individuals to begin tracing with
#' (`detected_id`), functions to control the discovery and termination
#' of tracing (`discovery_function`, `termination_function`), and parameters about
#' iteration (`max_attempts`) and fallback termination conditions
#' (`minimum_size`).
#'
#' @param detected_id A vector of individuals to begin tracing with
#' @param active A vector of active individuals
#' @param parents A matrix encoding the transmission history
#' @param minimum_size The minimum number of individuals to form a sample
#' @param discovery_function A function to determine which individuals to trace
#' @param max_attempts The maximum number of attempts to perform to obtain a sample, before failing.
#' @param termination_function A function to determine when to terminate tracing a contact tracing.
#' @return A list with three fields: "status", "samples", and "found"
#' @export
contact_tracing_engine <- function(detected_id,
                                   active, parents,
                                   discovery_function,
                                   max_attempts,
                                   termination_function,
                                   minimum_size = 3) {
    # Try to perform iterative contact tracing. If the length of the active
    # nodes is more than the minimum_size, collect into a list and return the output.
    # Else, try again.

    # Type safety
    # Convert to a vector. Handles both lists and single values.
    detected_id <- c(unlist(detected_id))
    attempts <- 1
    while (TRUE) {
        # try
        start_index <- (attempts %% length(detected_id)) + 1
        starting_id <- detected_id[[start_index]]
        result <- contact_tracing_core(
            detected_id = starting_id,
            active = active, parents = parents,
            discovery_function = discovery_function,
            termination_function = termination_function
        )
        if (length(result$samples) >= minimum_size) {
            result[["success"]] <- TRUE
            break
        }
        if (attempts >= max_attempts) {
            result[["success"]] <- FALSE
            break
        }
        # Incremenet the attempt counter and try again
        attempts <- attempts + 1
    }
    return(result)
}

#' Core contact tracing algorithm
#'
#' Given a set of active nodes and the transmision history (parents),
#' sample a set of nodes that are contact traced using the `discovery_function`
#' parameter to determine the transmission rate.
#'
#' @param detected_id The id of the node used to start the tracing
#' @param active A vector of active individuals
#' @param parents A matrix of transmission history
#' @param discovery_function A function that takes a list of nodes and relative
#'   transmission times and determines which will be included
#' @param termination_function A function that takes the list of discovered nodes
#'  and determines if the tracing should stop
#' @return A list with three fields: "found", "samples", and "status"
#'
#' @export
contact_tracing_core <- function(detected_id, active, parents,
                                 discovery_function, termination_function) {
    # Core algorithm for contact tracing
    # We'll build a list of all discovered nodes, and then sample from this list.
    found <- c(detected_id)
    seeds <- c(-1) # Store the seed we used to find each new infection
    # Allows us to be efficient (DAG property) and not have to search the entire
    # list of active nodes each time to preserve uniqueness. Use -1 as a flag
    # to indicate the initial detection event.

    index <- 1
    while (TRUE) {
        current_individual <- found[[index]]
        source_individual <- seeds[[index]] # Preserve uniqueness
        # get all all candidates from the tree

        all_can <- get_connections(seed = current_individual, parents = parents)
        # Determine which we will discover
        discovered <- discovery_function(
            parent_node = all_can[["parent"]],
            secondary_nodes = all_can[["secondary_infections"]],
            infection_time = all_can[["infection_time"]],
            secondary_time = all_can[["secondary_times"]]
        )
        # It's fine if we discover the same node twice, but we don't want to
        # contact trace from that node twice
        discovered <- discovered[discovered != source_individual]
        # Don't discover the root (0) of the inital infection (1)
        discovered <- discovered[discovered != 0]

        # Add the discovered nodes to the list of found nodes
        found <- c(found, discovered)
        # Record which node we used to find each new infection to prevent
        # duplications in our detection list
        seeds <- c(seeds, rep(current_individual, length(discovered)))
        # Check provided termination condition
        if (termination_function(found = found, active = active, index = index)) {
            status <- "Termination condition satisfied"
            break
        }

        # increment the index
        index <- index + 1
        # Check if we have any more discovered nodes to trace
        if (index > length(found)) {
            status <- "Out of nodes to explore"
            break
        }
        # If you made it this far, we have another node to explore
    }
    # Determine the nodes that are still active
    active_nodes <- active[active %in% found]
    return(list("found" = found, "samples" = active_nodes, status = status))
}

# Get the nodes connected to a (seed) node.
# Expect a parent and a list of offspring.
# Time data may not be known or reliable in pratice, but for simulations,
# This is useful information for the discovery function.
#' @export
get_connections <- function(seed, parents) {
    # get the parent of the seed
    parent <- parents[seed, 1]
    infection_time <- parents[seed, 2]
    # Get all offspring whose parent is seed
    if (seed %in% parents[, 1]) { # Check that we have data
        # This returns their absolute ID (row number)
        secondary_infections <- which(parents[, 1] == seed)
        # Get the infection times of secondary infections
        secondary_time <- parents[secondary_infections, 2]
    } else {
        # Report no secondary infections if the seed has no children
        secondary_infections <- NULL
        secondary_time <- NULL
    }
    return(list(
        parent = parent,
        infection_time = infection_time,
        secondary_infections = secondary_infections,
        secondary_times = secondary_time
    ))
}

#' A factory function to discover connections with uniform probability
#'
#' Obtain a discovery function that determines which adjacent nodes in the
#' contact network will be revealed with uniform probability.
#'
#' @param p The discovery probability. A float between 0 and 1.
#'
#' @return A function that takes a parent node, a list of secondary nodes,
#'  and their infection times, and returns a list of discovered nodes. Uses the
#'  uniform probability distribution to determine who is sampled.
#' @importFrom stats rbinom
#' @seealso  contact_traced_uniform_ids
#' @export
uniform_discovery_factory <- function(p) {
    # Uniform discovery function

    uniform_discovery <- function(parent_node, secondary_nodes,
                                  infection_time, secondary_times) {
        # Treat all nodes/contacts the same
        nodes <- c(secondary_nodes, parent_node)
        # Place secondary nodes first so we'll search through them first.
        # Use a uniform probability distribution for node discovery.
        mask <- as.logical(rbinom(n = length(nodes), prob = p, size = 1))
        # return the nodes that were discovered
        return(nodes[mask])
    }

    return(uniform_discovery)
}

#' Factory function for a detection to check termination condition based on
#' A sufficient amount of individuals
#' Terminate when we have found a `minimum_size` number of active individuals
#' @param minimum_size The minimum number of discovered active individuals to
#' terminate
#' @return A function that takes a list of found nodes, a list of active nodes,
#' and an index, and returns a boolean indicating if the termination condition
#' has been met. Returns true if the number of active nodes found is greater
#' than or equal to `minimum_size`.
#' @export
sufficient_data_data_factory <- function(minimum_size) {
    terminate_when_enough_data <- function(found, active, index) {
        active_found <- active[active %in% found]
        if (length(active_found) >= minimum_size) {
            return(TRUE)
        } else {
            return(FALSE)
        }
    }
    return(terminate_when_enough_data)
}

#' Factory function to never terminate contact tracing until all known nodes
#'   have been traced
#'
#' When contact tracing, we don't want to stop until we have traced all known.
#' This is computationally more expensive, but better describes reality.
#'
#' @return A function that takes a list of found nodes, a list of active nodes,
#' and an index, and returns a boolean indicating if the termination condition
#' has been met. Always returns False.
#' @export
never_terminate_early_factory <- function() {
    never_terminate_early <- function(found, active, index) {
        return(FALSE)
    }
    return(never_terminate_early)
}
