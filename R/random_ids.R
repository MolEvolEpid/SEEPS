# Random ids module for sampling a list of (possible) detections
# from a transmission history.

random_ids <- function(
    active, minimum_size = FALSE,
    proportional = FALSE, spike_root = FALSE) {
    # Either minimum_size or proportional must be specified.
    if (!(minimum_size || proportional)) {
        stop("Either a minimum size or proportional sampling must be specified")
    }
    if (proportional) {
        observation_size <- length(active) * proportional
        observation_size <- floor(observation_size) # enforce type safety
    }
    if (minimum_size > 0) observation_size <- as.integer(minimum_size)

    samples <- sample(active, observation_size, replace = FALSE)
    if (spike_root) {
        samples <- c(samples, 0) # add the root in here
        observation_size <- observation_size + 1 # increase the number of samples
    }
    return(list("observation_size" = observation_size, "samples" = samples))
}

#' Sample a proportion of the active individuals
#'
#' A helper function for sampling a proportion of the active individuals.
#' Select a proportion of individuals from the active set. Does not use
#' any information about the transmission history to select the sample.
#'
#' @param active A vector of active individuals
#' @param proportion A float between 0 and 1.
#'  The proportion of active individuals to sample.
#' @param minimum_size The minimum number of individuals to sample.
#'  If the proportional sample is smaller than the minimum size,
#'  the proprotional size is used.
#' @param spike_root Whether to include the root in the sample.
#' @seealso random_fixed_size_ids
#' @export
random_prop_ids <- function(active, proportion, minimum_size,
                            spike_root = FALSE) {
    return(random_ids(active,
        proportional = proportion,
        minimum_size = minimum_size
    ))
}

#' Sample a fixed number of of individuals randomly from the population.
#'
#' A helper function for sampling a fixed number of individuals from the
#' population. Select a fixed number of individuals from the active set.
#' Does not use any information about the transmission history or the
#' active population size to select the sample.
#'
#' @param active A vector of active individuals
#' @param minimum_size The minimum number of individuals to sample.
#' @param spike_root Whether to include the root in the sample.
#' @seealso random_prop_ids
#' @export
random_fixed_size_ids <- function(active, minimum_size, spike_root = FALSE) {
    return(random_ids(active,
        minimum_size = minimum_size,
        spike_root = spike_root
    ))
}
