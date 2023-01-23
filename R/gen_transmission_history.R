# The previous simulation. Here for QA & unit testing.
gen_transmission_history_exponential_constant_classic <- function(minimum_population, #nolint: object_length_linter
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
    # Top level/user function. Users should start with this function.

    parameters <- wrap_parameters(minimum_population, offspring_rate_fn,
                                  maximum_population_target, total_steps,
                                  spike_root)

    state <- initialize(parameters)
    validate_state(state)  # Check that the state is valid
    state <- gen_exp_phase(state, parameters)
    state <- gen_const_phase(state, parameters, num_steps=total_steps)
    return(clean_up(state))
}
#' Utility function for wrapping parameters
#' @export
#' @noRd
wrap_parameters <- function(minimum_population, offspring_rate_fn,
                            maximum_population_target, total_steps,
                            spike_root = FALSE) {
    parameters <- list("minimum_population" = minimum_population,
                       "offspring_rate_fn" = offspring_rate_fn,
                       "maximum_population_target" = maximum_population_target,
                       "total_steps" = total_steps,
                       "spike_root" = spike_root)
    return(parameters)
}

#' Utility function for initializing the simulation
#' @export
initialize <- function(parameters) {
    # Initialize the simulation parameters (a "state" in SEEPS parlance).
    # This is called once at the beginning of the simulation.

    parents <- matrix(0, max(parameters$total_steps, 12) * parameters$maximum_population_target, 2)

    # Declare counters and pointers
    # start_step <- 0
    active_index <- 2
    curr_step <- 1
    # Setup the index case
    active <- c(1)
    birth_step <- c(1)
    end_step <- sample(13:36, 1)
    # add length_active
    length_active <- length(active)
    # Pack it together and return it
    return(list("parents" = parents, "active" = active,
                "active_index" = active_index, "end_step"= end_step,
                "birth_step" = birth_step, "length_active" = length_active,
                "curr_step" = curr_step, "total_offspring" = 1))
}

#' Simulate an exponential growth.
#' @export
gen_exp_phase <- function(state, parameters) {
    # Perform exponential growth until the population reaches the target size.
    # Then return the new state of the simulation.

    while(TRUE) {
        state <- step(state, parameters)
        # Check if we can stop
        if (state$length_active > 0.9 * parameters$maximum_population_target
            && state$length_active > parameters$minimum_population) {

            # If both, we have sufficient population to draw a sample,
            # and we have reached the target population size, we can stop
            # the exponential growth phase.
            break
        }
    }
    return(state)
}
#' Simulate a fixed number of steps
#'
#' Usually, we want to simulate a fixed number of steps after the exponential growth phase.
#' This function performs that step.
#'
#' @export
gen_const_phase <- function(state, parameters, num_steps) {
    # Perform constant population growth until the simulation is complete.
    # Specify the number of steps to perform.
    if (num_steps > 0) {
        # Do the steps
        for (i in 1:num_steps) {
            print("Performing a step")
            state <- step(state, parameters)
        }
    }

    return(state)
}

#' Clean up the simulation state and return the relevant data.
#'
#' @export
clean_up <- function(state){
    # Extract the relevant data from the simulation and return it. We do not need
    # all of the internal state to be exposed to the user.
    # This is called once at the end of the simulation.

    returned_data <- list("parents" = state$parents, "active" = state$active,
                          "t_end" = state$curr_step,
                          "total_offspring" = state$tot_offsprings)
    return(returned_data)
}

#' Perform a single step of the simulation
#'
#' @export
step <- function(state, parameters) {
    # Perform a single step of the simulation. This is the core algorithm. It is used
    # by both the exponential and constant phases which provide hyperparameter inputs
    # This is called repeatedly until the simulation is complete. This represents a single step.

    state$length_active <- length(state$active)  # pre-compute for efficiency
    if (parameters$maximum_population_target - state$length_active > 0) {  # If we can/need to grow
        # Sample offspring distribution
        offsprings <- parameters$offspring_rate_fn(current_step = state$curr_step, birth_step = state$birth_step)
        # Cap offspring at the target
        state$tot_offsprings <- ceiling(min(parameters$maximum_population_target - state$length_active, sum(offsprings)))
        # Assign parents using the sampled rates/weights
        prts <- sample(1:state$length_active, state$tot_offsprings, replace = TRUE, prob = offsprings)
        # Record parents
        if ((state$active_index + state$tot_offsprings - 1) > dim(state$parents)[1]) {
            # Double the table size. This should not need to occur, but will break the simulation
            # if we overflow the table. Default values are chosen with HIV dynamics in mind.
            state$parents <- rbind(state$parents, matrix(0, max(parameters$total_steps, 12) * parameters$maximum_population_target, 2))
        }
        state$parents[state$active_index:(state$active_index + state$tot_offsprings - 1), 1] <- state$active[prts]
        state$parents[state$active_index:(state$active_index + state$tot_offsprings - 1), 2] <- state$curr_step

        # Update state vectors needed for stepping the model
        state$active <- c(state$active, state$active_index:(state$active_index + state$tot_offsprings - 1))
        state$birth_step <- c(state$birth_step, rep(state$curr_step + 1, state$tot_offsprings))
        state$end_step <- c(state$end_step, state$curr_step + sample(13:36, state$tot_offsprings, replace = TRUE))

    } else {
        # If we can't grow, set the update increments to 0
        state$tot_offsprings <- 0
    }
    # Remove old individuals from active population
    inds_rem <- which(state$end_step <= state$curr_step)  # remove anything we should
    len_rem <- length(inds_rem)
    # If the removal would remove the entire population, don't remove any individuals and do another step.
    # Safeguard for when R0 is close to or less than 1 and the population is small to avoid accidental extinction.
    if ((len_rem > 0) && (len_rem < state$length_active)) {
        # Only remove if we have at least one, and we are not removing the entire population...
        state$active <- state$active[-inds_rem]
        state$end_step <- state$end_step[-inds_rem]
        state$birth_step <- state$birth_step[-inds_rem]
    }
    # Update global indexes/pointers
    state$curr_step <- state$curr_step + 1
    state$active_index <- state$active_index + state$tot_offsprings
    # return the new state
    return(state)
}

#' Remove samples from the simulation
#'
#' Remove a set of samples from the simulation. Intended to be used with contact tracing, we expect
#' individuals identified through contact tracing to be removed from the population of active (uncontrolled or undiagnosed)
#' cases.
#'
#' Provide a state object and a vector of samples to remove. A new state object is returned with the samples removed.
#'
#' @export
remove_samples <- function(state, samples) {
    # 'Remove' samples from the simulation. This is done by setting their end_step to the current time.
    # This emulates the removal of individuals due to sampling.
    ids_rem <- which(state$active %in% samples)
    state$active <- state$active[-ids_rem]
    state$end_step <- state$end_step[-ids_rem]
    state$birth_step <- state$birth_step[-ids_rem]
    return(state)
}

#' Utility function to check that a state is (computationally) valid
#'
#' Check that the state object is a valid and has the correct list elements.
#'
#' @export
validate_state <- function(state) {
    # Check that the state list has the correct types.
    # We're trying to avoid the OO system and keep data in a list between transforms.

    # Check that the state is a list
    if (!is.list(state)) {
        stop("State must be a list.")
    }
    # Expect that state has all the attributes declared in initialize
    if (!all(c("parents", "active", "active_index", "end_step", "birth_step", "length_active", "curr_step", "total_offspring") %in% names(state))) {
        stop("State is missing one or more of the required attributes.")
    }
    # Check that the attributes are the correct type
    if (!is.matrix(state$parents)) {
        stop("Parents must be a matrix.")
    }
    if (!is.integer(state$active)) {
        stop("Active must be an integer vector.")
    }
    if (!is.integer(state$active_index)) {
        stop("Active index must be an integer.")
    }
    if (!is.integer(state$end_step)) {
        stop("End step must be an integer vector.")
    }
    if (!is.integer(state$birth_step)) {
        stop("Birth step must be an integer vector.")
    }
    if (!is.integer(state$curr_step)) {
        stop("The current step must be an integer.")
    }
    if (!is.integer(state$total_offspring)) {
        stop("Total offspring must be an integer.")
    }
    if (!is.integer(state$length_active)) {
        stop("length_active must be an integer.")
    }
    # Check that only these attributes are present
    if (length(names(state)) != 8) {
        stop("State has the wrong number of values.")
    }
    # length_active

}
