#' Generate transmission history for full population
#'
#' Forward stochastic simulation to generate transmission history for a
#' complete population with an exponential growth then constant population structure.
#'
#' This simulation is intended for HIV, but may be broadly applicable if the lifespan distribution
#' is adjusted to other distributions/parameters.
#'
#' @export
gen_transmission_history_exponential_constant <- function(minimum_population, #nolint: object_length_linter
        offspring_rate_fn, maximum_population_target, total_steps,
                              spike_root = FALSE) {

    # Store results in here. This is usually large enough,
    parents <- matrix(0, max(total_steps, 12) * maximum_population_target, 2)

    # Declare counters and pointers
    flag <- TRUE
    start_step <- 1000
    active_index <- 2
    curr_step <- 1
     # Setup the index case
    active <- c(1)
    birth_step <- c(1)
    end_step <- sample(13:36, 1)
    while (curr_step < total_steps + start_step) {

        inds_rem <- NULL  # Indices to remove

        length_active <- length(active)  # pre-compute for efficiency
        if (flag && length_active > 0.9 * maximum_population_target && length_active > minimum_population) {
            # Third check needed for edge case of min_pop >= 0.9 maximum_population_target
            flag <- FALSE  # Do this chunk of code exactly once
            start_step <- curr_step  # Update the starting step of the constant phase
        }

        if (maximum_population_target - length_active > 0) {  # If we can/need to grow
            # Sample offspring distribution
            offsprings <- offspring_rate_fn(current_step = curr_step, birth_step = birth_step)
            # Cap offspring at the target
            tot_offsprings <- ceiling(min(maximum_population_target - length_active, sum(offsprings)))
            # Assign parents using the sampled rates/weights
            prts <- sample(1:length_active, tot_offsprings, replace = TRUE, prob = offsprings)
            # Record parents
            if ((active_index + tot_offsprings - 1) > dim(parents)[1]) {
                # Double the table size. This should not need to occur, but will break the simulation
                # if we overflow the table. Default values are chosen with HIV dynamics in mind.
                parents <- rbind(parents, matrix(0, max(total_steps, 12) * maximum_population_target, 2))
            }
            parents[active_index:(active_index + tot_offsprings - 1), 1] <- active[prts]
            parents[active_index:(active_index + tot_offsprings - 1), 2] <- curr_step
            # Update state vectors needed for stepping the model
            active <- c(active, active_index:(active_index + tot_offsprings - 1))
            birth_step <- c(birth_step, rep(curr_step + 1, tot_offsprings))
            end_step <- c(end_step, curr_step + sample(13:36, tot_offsprings, replace = TRUE)

            )
        } else {
            # If we can't grow, set the update increments to 0
            tot_offsprings <- 0
        }

        # Remove old individuals from active population
        inds_rem <- which(end_step <= curr_step)  # remove anything we should
        len_rem <- length(inds_rem)
        # If the removal would remove the entire population, don't remove any individuals and do another step.
        # Safeguard for when R0 is close to or less than 1 and the population is small to avoid accidental extinction.
        if (len_rem > 0 && len_rem < length_active) {
            active <- active[-inds_rem]
            end_step <- end_step[-inds_rem]
            birth_step <- birth_step[-inds_rem]
        }
        # Update global indexes/pointers
        curr_step <- curr_step + 1
        active_index <- active_index + tot_offsprings
    }

    returned_data <- list("parents" = parents, "active" = active,
                          "t_end" = curr_step, "total_offspring" = tot_offsprings)

    return(returned_data)
}
