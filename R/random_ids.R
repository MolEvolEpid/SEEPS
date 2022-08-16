# Helper function for sampling. If proportional, we specify a proportion of the active.


# TODO Candidate to split function into two
#' Generate random ids through sampling
#'
#' @export
random_ids <- function(active,minimum_size = FALSE,
proportional = FALSE, spike_root = FALSE) {
    # Either minimum_size or proportional must be specified.
    if (!(minimum_size || proportional)) {
        stop("Either a minimum size or proportional sampling must be specified")
    }
    if (proportional) {
        print("Trying to set the number of samples dynamically")
        observation_size <- length(active) * proportional
        observation_size <- floor(observation_size)  # enforce type safety

    }
    if (minimum_size > 0) observation_size <- as.integer(minimum_size)

    samples <- sample(active, observation_size, replace = FALSE)
    if (spike_root) {
        samples <- c(samples, 0)  # add the root in here
        observation_size <- observation_size + 1  # increase the number of samples
    }
    return(list("observation_size" = observation_size, "samples" = samples))
}