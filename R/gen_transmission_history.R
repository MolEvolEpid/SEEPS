#' A fallback simulation routine
#' Depricated code. May be removed at any point. New code should use the refactored API, or
#' `SEEPS::gen_transmission_history_exponential_constant` for a supported version.
#'
#' @param minimum_population The minimum population size. The simulation will not terminate until this size is reached.
#' @param offspring_rate_fn A function that returns the offspring rate for each individual at each time step.
#' @param maximum_population_target The maximum population size.
#' @param total_steps The total number of time steps to simulate after the exponential growth phase.
#' @param spike_root Ignored. Default FALSE.
#' @export
gen_transmission_history_exponential_constant_classic <- function(minimum_population, # nolint: object_length_linter
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
        inds_rem <- NULL # Indices to remove

        length_active <- length(active) # pre-compute for efficiency
        if (flag && length_active > 0.9 * maximum_population_target && length_active > minimum_population) {
            # Third check needed for edge case of min_pop >= 0.9 maximum_population_target
            flag <- FALSE # Do this chunk of code exactly once
            start_step <- curr_step # Update the starting step of the constant phase
        }

        if (maximum_population_target - length_active > 0) { # If we can/need to grow
            # Sample offspring distribution
            offsprings <- offspring_rate_fn(current_step = curr_step, birth_step = birth_step)
            # Cap offspring at the target
            total_offsprings <- ceiling(min(maximum_population_target - length_active, sum(offsprings)))
            # Assign parents using the sampled rates/weights
            prts <- sample(1:length_active, total_offsprings, replace = TRUE, prob = offsprings)
            # Record parents
            if ((active_index + total_offsprings - 1) > dim(parents)[1]) {
                # Double the table size. This should not need to occur, but will break the simulation
                # if we overflow the table. Default values are chosen with HIV dynamics in mind.
                parents <- rbind(parents, matrix(0, max(total_steps, 12) * maximum_population_target, 2))
            }
            parents[active_index:(active_index + total_offsprings - 1), 1] <- active[prts]
            parents[active_index:(active_index + total_offsprings - 1), 2] <- curr_step
            # Update state vectors needed for stepping the model
            active <- c(active, active_index:(active_index + total_offsprings - 1))
            birth_step <- c(birth_step, rep(curr_step + 1, total_offsprings))
            end_step <- c(end_step, curr_step + sample(13:36, total_offsprings, replace = TRUE))
        } else {
            # If we can't grow, set the update increments to 0
            total_offsprings <- 0
        }

        # Remove old individuals from active population
        inds_rem <- which(end_step <= curr_step) # remove anything we should
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
        active_index <- active_index + total_offsprings
    }

    returned_data <- list(
        "parents" = parents, "active" = active,
        "t_end" = curr_step, "total_offspring" = total_offsprings
    )

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
#' @param minimum_population The minimum population size. The simulation will not terminate until this size is reached.
#' @param offspring_rate_fn A function that returns the offspring rate for each individual at each time step.
#' @param maximum_population_target The maximum population size.
#' @param total_steps The total number of time steps to simulate after the exponential growth phase.
#' @param spike_root Passed to `wrap_parameters`. Default FALSE.
#' @export
gen_transmission_history_exponential_constant <- function(minimum_population, # nolint: object_length_linter
                                                          offspring_rate_fn, maximum_population_target, total_steps,
                                                          spike_root = FALSE) {
    # Top level/user function. Users should start with this function.

    parameters <- wrap_parameters(
        minimum_population, offspring_rate_fn,
        maximum_population_target, total_steps,
        spike_root
    )

    state <- initialize(parameters)
    validate_state(state) # Check that the state is valid
    state <- gen_exp_phase(state, parameters)
    state <- gen_const_phase(state, parameters, num_steps = total_steps)
    return(clean_up(state))
}

#' Utility function for wrapping parameters
#'
#' New users should review the tutorial vignette [Simulation_API] for a more
#' detailed discussion.
#'
#' @param minimum_population The minimum population size. The simulation will not terminate until this size is reached.
#' @param offspring_rate_fn A function that returns the offspring rate for each individual at each time step.
#' @param maximum_population_target The maximum population size.
#' @param total_steps The total number of time steps to simulate after the exponential growth phase.
#' @param spike_root Whether or not to keep the root infection at time zero in the phylogeny. Default FALSE.
#' May lead to unstable behavior, and is not fully implemented.
#'
#' @export
wrap_parameters <- function(minimum_population, offspring_rate_fn,
                            maximum_population_target, total_steps,
                            spike_root = FALSE) {
    parameters <- list(
        "minimum_population" = minimum_population,
        "offspring_rate_fn" = offspring_rate_fn,
        "maximum_population_target" = maximum_population_target,
        "total_steps" = total_steps,
        "spike_root" = spike_root
    )
    return(parameters)
}

#' Utility function for initializing the simulation
#'
#' @param parameters A list of parameters to run the simulation.
#' See [SEEPS::wrap_parameters] for the expected structure.
#'
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
    return(list(
        "parents" = parents, "active" = active,
        "active_index" = active_index, "end_step" = end_step,
        "birth_step" = birth_step, "length_active" = length_active,
        "curr_step" = curr_step, "total_offspring" = 1
    ))
}

#' Simulate an exponential growth.
#'
#' A utility function for simulating exponential growth. The simulation is
#' allowed to grow until the population reaches the target size.
#' Specifically, the simulation terminates once the population exceeds 90% of
#' the target size, and `parameters$minimum_population`individuals are present.
#'
#' @param state The state object (list) from the simulation.
#' See `initialize` for the expected structure.
#' @param parameters A list of parameters to run the simulation.
#' See `wrap_parameters` for the expected structure.
#' @param all_states A boolean flag to return all states during the simulation.
#' Default `FALSE`. If `TRUE`, returns a list of all states during the simulation.
#' If `FALSE`, returns only the final state.
#' @return Either a state or list of states, depending on `all_state` flag.
#' @export
gen_exp_phase <- function(state, parameters, all_states = FALSE) {
    # Perform exponential growth until the population reaches the target size.
    # Then return the new state of the simulation.
    if (all_states) {
        state_list <- list()
    }
    while (TRUE) {
        state <- step(state, parameters)
        if (all_states) {
            state_list[[length(state_list) + 1]] <- state
        }
        # Check if we can stop
        if (
            state$length_active > 0.9 * parameters$maximum_population_target &&
                state$length_active > parameters$minimum_population) {
            # If both, we have sufficient population to draw a sample,
            # and we have reached the target population size, we can stop
            # the exponential growth phase.
            break
        }
    }
    if (all_states) {
        return(state_list)
    } else {
        return(state)
    }
}

#' Simulate a fixed number of steps
#'
#' A utility function for simulating a fixed number of time steps. Each step is
#' performed regardless of the number of active individuals.
#'
#' @param state The state object (list) from the simulation.
#' See `initialize` for the expected structure.
#' @param parameters A list of parameters to run the simulation.
#' See `wrap_parameters` for the expected structure.
#' @param all_states A boolean flag to return all states during the simulation.
#' Default `FALSE`. If `TRUE`, returns a list of all states during the simulation.
#' If `FALSE`, returns only the final state.
#' @param num_steps The number of steps to perform. Must be a positive integer.
#' @param verbose A boolean flag to print a message at each step. Default `FALSE`.
#' @return Either a state or list of states, depending on `all_state` flag.
#'
#' @export
gen_const_phase <- function(state, parameters, num_steps, verbose = FALSE, all_states = FALSE) {
    # Perform constant population growth until the simulation is complete.
    # Specify the number of steps to perform.
    if (all_states) {
        state_list <- list()
    }
    if (num_steps > 0) {
        # Do the steps
        for (i in 1:num_steps) {
            if (verbose) {
                print("Performing a step")
            }

            state <- step(state, parameters)
            if (all_states) {
                state_list[[length(state_list) + 1]] <- state
            }
        }
    }
    if (all_states) {
        return(state_list)
    } else {
        return(state)
    }
}

#' Clean up the simulation state and return the relevant data.
#'
#' @param state The state object (list) from the simulation. See `initialize` for the expected structure.
#' @export
clean_up <- function(state) {
    # Extract the relevant data from the simulation and return it. We do not need
    # all of the internal state to be exposed to the user.
    # This is called once at the end of the simulation.

    returned_data <- list(
        "parents" = state$parents, "active" = state$active,
        "t_end" = state$curr_step,
        "total_offspring" = state$total_offsprings
    )
    return(returned_data)
}

#' Perform a single step of the simulation
#'
#' The workhorse function of SEEPS. Step forward a state one time step.
#'
#' @param state The state object (list) from the simulation. See `initialize` for the expected structure.
#' @param parameters A list of parameters to run the simulation. See `wrap_parameters` for the expected structure.
#' @export
step <- function(state, parameters) {
    # Perform a single step of the simulation. This is the core algorithm. It is used
    # by both the exponential and constant phases which provide hyperparameter inputs
    # This is called repeatedly until the simulation is complete. This represents a single step.

    state$length_active <- length(state$active) # pre-compute for efficiency
    if (parameters$maximum_population_target - state$length_active > 0) { # If we can/need to grow
        # Sample offspring distribution
        offsprings <- parameters$offspring_rate_fn(current_step = state$curr_step, birth_step = state$birth_step)
        # Cap offspring at the target
        state$total_offsprings <- ceiling(min(
            parameters$maximum_population_target - state$length_active,
            sum(offsprings)
        ))
        # Assign parents using the sampled rates/weights
        prts <- sample(1:state$length_active, state$total_offsprings, replace = TRUE, prob = offsprings)
        # Record parents
        if ((state$active_index + state$total_offsprings - 1) > dim(state$parents)[1]) {
            # Double the table size. This should not need to occur, but will break the simulation
            # if we overflow the table. Default values are chosen with HIV dynamics in mind.
            state$parents <- rbind(
                state$parents,
                matrix(0, max(parameters$total_steps, 12) * parameters$maximum_population_target, 2)
            )
        }
        state$parents[state$active_index:(state$active_index + state$total_offsprings - 1), 1] <- state$active[prts]
        state$parents[state$active_index:(state$active_index + state$total_offsprings - 1), 2] <- state$curr_step

        # Update state vectors needed for stepping the model
        state$active <- c(state$active, state$active_index:(state$active_index + state$total_offsprings - 1))
        state$birth_step <- c(state$birth_step, rep(state$curr_step + 1, state$total_offsprings))
        state$end_step <- c(state$end_step, state$curr_step + sample(13:36, state$total_offsprings, replace = TRUE))
    } else {
        # If we can't grow, set the update increments to 0
        state$total_offsprings <- 0
    }
    # Remove old individuals from active population
    inds_rem <- which(state$end_step <= state$curr_step) # remove anything we should
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
    state$active_index <- state$active_index + state$total_offsprings
    # return the new state
    return(state)
}

#' Remove samples from the simulation
#'
#' Remove a set of samples from the simulation. Intended to be used with contact
#' tracing, we expect individuals identified through contact tracing to be
#' removed from the population of active (uncontrolled or undiagnosed) cases.
#'
#' Provide a state object and a vector of samples to remove. A new state object
#' is returned with the samples removed.
#'
#' @param state The state object (list) from the simulation. See `initialize`
#'  for the expected structure.
#' @param samples A vector of individual id's to remove from the simulation.
#'
#' @export
remove_samples <- function(state, samples) {
    # 'Remove' samples from the simulation. This is done by setting their
    # end_step to the current time.
    # This emulates the removal of individuals due to sampling.
    ids_rem <- which(state$active %in% samples)
    state$active <- state$active[-ids_rem]
    state$end_step <- state$end_step[-ids_rem]
    state$birth_step <- state$birth_step[-ids_rem]
    return(state)
}

#' Keep a set of samples in the simulation
#'
#' Remove all other samples from the simulation. Intended to be used to represent a masking event,
#' Where only a subset of the population is propagated forward in time.
#'
#' Provide a state object and a vector of samples to keep. A new state object is returned with only the samples kept.
#'
#' @param state The state object (list) from the simulation. See `initialize` for the expected structure.
#' @param samples A vector of individual id's to keep in the simulation. All other id's are removed.
#'
#' @return A new state object with only the samples kept.
#' @export
keep_samples <- function(state, samples) {
    # 'Keep' samples in the simulation. This is done by setting the end_step of all other samples to the current time.
    # This emulates the removal of individuals due to sampling.
    ids_rem <- which(!(state$active %in% samples))
    state$active <- state$active[-ids_rem]
    state$end_step <- state$end_step[-ids_rem]
    state$birth_step <- state$birth_step[-ids_rem]
    return(state)
}

#' Utility function to check that a state is (computationally) valid
#'
#' Check that the state object is a valid and has the correct list elements.
#'
#' @param state The state object (list) from the simulation. See `initialize` for the expected structure.
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
    if (!all(c(
        "parents", "active", "active_index", "end_step", "birth_step",
        "length_active", "curr_step", "total_offspring"
    ) %in% names(state))) {
        stop("State is missing one or more of the required attributes.")
    }
    # Check that the attributes are the correct type
    if (!is.matrix(state$parents)) {
        stop("Parents must be a matrix.")
    }
    if (!is_integer(state$active)) {
        warning("Active must only contain integers. This is not enforced by the type system.")
    }
    if (!is_integer(state$active_index)) {
        warning("Active index must only contain integers. This is not enforced by the type system.")
    }
    if (!is_integer(state$end_step)) {
        warning("End step must only contain integers. This is not enforced by the type system.")
    }
    if (!is_integer(state$birth_step)) {
        warning("Birth step must only contain integers. This is not enforced by the type system.")
    }
    if (!is_integer(state$curr_step)) {
        warning("The current step must be an integer. This is not enforced by the type system.")
    }
    if (!is_integer(state$total_offspring)) {
        warning("Total offspring must be an integer. This is not enforced by the type system.")
    }
    if (!is_integer(state$length_active)) {
        warning("length_active must be an integer. This is not enforced by the type system.")
    }
    # Check that only these attributes are present
    if (length(names(state)) != 8) {
        warning("State has the wrong number of values. This is not enforced by the type system.")
    }
}

is_integer <- function(x) {
    all(x == floor(x))
}
