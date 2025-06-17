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
                                                          offspring_rate_fn,
                                                          maximum_population_target, total_steps,
                                                          spike_root = FALSE, ...) {
    # Top level/user function. Users should start with this function.
    parameters <- wrap_parameters(
        minimum_population = minimum_population,
        offspring_rate_fn = offspring_rate_fn,
        maximum_population_target = maximum_population_target,
        total_steps = total_steps, spike_root = spike_root, ...
    )

    state <- SEEPS::initialize(parameters)
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
#' Take optional arguments for rate_selector and lifespan_sampler_selector
#'
#' @param randomize_parental_assignments If `TRUE`, offspring will be assigned to parents randomly. Default is `TRUE`.
#' Also round the number of secondary infections generated up to the nearest integer. This allows for the simulation to
#' have individuals contributing fractional infections to the total number of infections. This approximation is
#' asymtpotically correct but should not be used for small populations. Set to `FALSE`, and use a dispersion sampler to
#' randomize the number of infections that each individual generates per step. See also
#' [SEEPS::poisson_dispersion_sampler()].
#' @export
wrap_parameters <- function(minimum_population, offspring_rate_fn,
                            maximum_population_target,
                            total_steps, randomize_parental_assignments = TRUE,
                            spike_root = FALSE, allow_extinction = FALSE, ...) {
    parameters <- list(
        "minimum_population" = minimum_population,
        "offspring_rate_fn" = offspring_rate_fn,
        "maximum_population_target" = maximum_population_target,
        "total_steps" = total_steps,
        "spike_root" = spike_root,
        "randomize_parental_assignments" = randomize_parental_assignments,  # If individuals do not contribute integer
        "allow_extinction" = allow_extinction
    )

    # if we have a rate_selector, add it to the parameters
    params <- list(...)

    # now, we're going to check if there are any extra arguments. For most of these,
    # we provide default behavior, so the user does not need to configure them.
    # these presets should have minimum complexity (constants, uniform, poisson, etc)

    # First, check to see if we have rates and lifespan samplers
    # if we have a lifespan_sampler, add it to the parameters
    if ("lifespan_sampler" %in% names(params)) {
        if (!is.list(params$lifespan_sampler)) {
            parameters$lifespan_sampler <- list(params$lifespan_sampler)
        } else {  # is a list
            parameters$lifespan_sampler <- params$lifespan_sampler
        }
    } else {  # equip a default lifespan_sampler
        message("No lifespan_sampler function provided. Assigning uniform lifespans on 13 to 36 months.")
        parameters$lifespan_sampler <- list(SEEPS::uniform_lifespan_sampler(13, 36))
    }

    if (!is.list(parameters$offspring_rate_fn)) {
        parameters$offspring_rate_fn <- list(parameters$offspring_rate_fn)
    } else {
        parameters$offspring_rate_fn <- parameters$offspring_rate_fn
    }

    if ("rate_selector" %in% names(params)) {
        parameters$rate_selector_function <- params$rate_selector
    } else {
        message("No rate_selector function provided. Assigning offspring their parent's type")
        # Use a default that will only be used if the user does not provide one
        parameters$rate_selector_function <- SEEPS::assign_parental_factory(parameters$offspring_rate_fn)
    }

    if ("lifespan_selector" %in% names(params)) {
        parameters$lifespan_selector <- params$lifespan_selector
    } else {
        message("No lifespan_selector function provided. Assigning offspring their parent's type.")
        parameters$lifespan_selector <- SEEPS::assign_parental_factory(parameters$lifespan_sampler)
    }


    if ("dispersion_sampler" %in% names(params)) {
        # check if we got a list. If not, we need to wrap it into a list
        if (!is.list(params$dispersion_sampler)) {
            params$dispersion_sampler <- list(params$dispersion_sampler)
        } else {
            parameters$dispersion_sampler <- params$dispersion_sampler
        }
    } else {
        message("No dispersion_sampler function provided. Assigning zero variance dispersion.")
        parameters$dispersion_sampler <- list(function(means, current_step, birth_step, params = list()) {
            return(means)
        })
    }

    if ("dispersion_selector" %in% names(params)) {
        parameters$dispersion_selector <- params$dispersion_selector
    } else {
        message("No dispersion_selector function provided. Assigning offspring their parent's type.")
        parameters$dispersion_selector <- SEEPS::assign_parental_factory(parameters$dispersion_sampler)
    }

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

    # Store parental linkages
    parents <- matrix(0, max(parameters$total_steps, 12) * parameters$maximum_population_target, 2)
    # Store the types of lifespans and spread rates
    lifespan_types <- matrix(0, max(parameters$total_steps, 12) * parameters$maximum_population_target, 1)
    # Control which offspring rate function is used
    spread_rate_types <- matrix(0, max(parameters$total_steps, 12) * parameters$maximum_population_target, 1)
    # similarly, we'll need to set the dispersion rate types
    dispersion_rate_types <- matrix(0, max(parameters$total_steps, 12) * parameters$maximum_population_target, 1)

    # set the type for the first individual to 1
    lifespan_types[1] <- 1
    spread_rate_types[1] <- 1
    dispersion_rate_types[1] <- 1

    # Declare counters and pointers
    active_index <- 2
    curr_step <- 1
    # Setup the index case
    active <- c(1)
    birth_step <- c(1)
    end_step <- parameters$lifespan_sampler[[1]](0, 0, c(1))
    # add length_active
    length_active <- length(active)
    # Pack it together and return it
    state <- list(
        "parents" = parents, "active" = active,
        "active_index" = active_index, "end_step" = end_step,
        "birth_step" = birth_step, "length_active" = length_active,
        "curr_step" = curr_step, "total_offspring" = 1,
        "spread_rate_types" = spread_rate_types,
        "dispersion_rate_types" = dispersion_rate_types,
        "lifespan_types" = lifespan_types,
        "extinct" = FALSE
    )
    return(state)
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
    num_steps <- 0
    if (state$extinct) {
        warning("The simulation is already extinct. No further steps can be taken. Exiting.")
        # if we are extinct, we can't do any more steps. Break the loop early.
        if (all_states) {
            state_list[[length(state_list) + 1]] <- state
            return(state_list)
        } else {
            return(state)
        }
    }
    while (TRUE) {
        state <- step(state, parameters)
        num_steps <- num_steps + 1
        if (all_states) {
            state_list[[length(state_list) + 1]] <- state
        }
        # check if we went extinct
        if (state$extinct) {
            # if we are extinct, we can't do any more steps. Break the loop early.
            break
        }
        # Check if we can stop
        if (state$length_active > 0.9 * parameters$maximum_population_target &&
                state$length_active >= parameters$minimum_population) {
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
        # first, check if we are extinct
        if (state$extinct) {
            warning("The simulation is already extinct. No further steps can be taken. Exiting.")
            # if we are extinct, we can't do any more steps. Break the loop early.
            if (all_states) {
                return(state_list)
            } else {
                return(state)
            }
            # if we are not, do the loop
        }

        for (i in 1:num_steps) {
            if (verbose) {
                print("Performing a step")
            }

            state <- step(state, parameters)
            if (all_states) {
                state_list[[length(state_list) + 1]] <- state
            }
            # check if we went extinct
            if (state$extinct) {
                # if we are extinct, we can't do any more steps. Break the loop early.
                break
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
#' @export
clean_up <- function(state, trim = FALSE, ...) {
    # Extract the relevant data from the simulation and return it. We do not need
    # all of the internal state to be exposed to the user.
    # This is called once at the end of the simulation.

    # Drop the last rows of the parents matrix that are all zeros
    if (isTRUE(trim)) {
        state$parents <- state$parents[1:state$trim, ]
        # "lifespan_types" = state$lifespan_types,
        lifespan_types <- state$lifespan_types[1:state$trim]
        # "spread_rate_types" = state$spread_rate_types,
        spread_rate_types <- state$spread_rate_types[1:state$trim]
        # "dispersion_rate_types" = state$dispersion_rate_types
        dispersion_rate_types <- state$dispersion_rate_types[1:state$trim]
    } else {
        parents <- state$parents
        lifespan_types <- state$lifespan_types
        spread_rate_types <- state$spread_rate_types
        dispersion_rate_types <- state$dispersion_rate_types
    }
    returned_data <- list(
        "parents" = parents, "active" = state$active,
        "t_end" = state$curr_step,
        "total_offspring" = state$tot_offsprings,
        "lifespan_types" = lifespan_types,
        "spread_rate_types" = spread_rate_types,
        "dispersion_rate_types" = dispersion_rate_types
    )
    return(returned_data)
}

#' Step function for the simulation
#' @export
step <- function(state, parameters, ...) {
    state$length_active <- length(state$active)  # pre-compute for efficiency
    pop_to_cap <- parameters$maximum_population_target - state$length_active  # can add up to this many individuals
    # TODO: Extract this into a separate function
    if (pop_to_cap > 0) {
        # Calculate the expected number of offspring
        offsprings <- compute_offsprings(state, parameters, ...)
        offsprings_dispersed <- compute_dispersion(offsprings, state, parameters, ...)  # RNG HERE, possibly
        # determine how many offspring can be generated
        state$tot_offsprings <- ceiling(min(pop_to_cap,  sum(offsprings_dispersed)))
        # Assign parents using the sampled rates/weights
        # parents stores the index for which individual in state$active is the parent
        if (state$tot_offsprings > 0) {
            if (isTRUE(parameters$randomize_parental_assignments)) {
                # assign the parents based on who is available
                parents <- sample(1:state$length_active, state$tot_offsprings,
                                  replace = TRUE, prob = offsprings_dispersed)  # RNG HERE
            } else {  # assume offsprings_dispersed stores how many offspring each individual generates

                parents <- rep(seq_along(offsprings_dispersed), offsprings_dispersed)
                # shuffle them

                # don't allow sampling of more offspring than we can have
                pop_to_resample <- min(length(parents), state$tot_offsprings)
                # Population cap should be enforced here.
                parents <- sample(parents, pop_to_resample)  # RNG HERE

            }
        } else {
            parents <- NULL
        }
        #END TODO: Extract this into a separate function
        state <- check_safe_sizes(state, parameters, ...)  # NO RNG
        state <- assign_inherit(parents, state, parameters, ...)  # RNG HERE, possibly
        state <- record_parents(parents, state, parameters, ...)  # NO RNG
        state <- update_vital_statistics(parents, state, parameters, ...)  # RNG HERE

    } else {
        state$tot_offsprings <- 0
    }

    state <- remove_expiring_individuals(state, parameters$allow_extinction)
    state <- step_sim_indexes(state)
    return(state)
}


#' @export
step_sim_indexes <- function(state) {
    state$curr_step <- state$curr_step + 1
    state$active_index <- state$active_index + state$tot_offsprings
    return(state)
}


#' @param state The current state of the simulation
#' @param allow_extinction If `TRUE`, the simulation will allow for the population to go extinct. Default is `FALSE`.
#' If FALSE, will not remove older individuals from the simulation until the population produces an offspring.``
#' @export
remove_expiring_individuals <- function(state, allow_extinction = FALSE) {
    #TODO: Check me
    # Remove individuals who have reached the end of their lifespan
    inds_rem <- which(state$end_step <= state$curr_step)  # remove anything we should
    len_rem <- length(inds_rem)
    # If the removal would
    # remove
    # the entire population, don't remove any individuals and do another step.
    # Safeguard for when R0 is close to or less than 1 and the population is small to avoid accidental extinction.
    # (len_rem > 0) && (len_rem < state$length_active)

    # print out some diagnostic info about the current state
    if (len_rem == state$length_active && allow_extinction) {
        state$extinct <- TRUE
    }
    if ((len_rem < state$length_active) && (allow_extinction == FALSE && len_rem > 0)) {
        # Only remove if we have at least one, and we are not removing the entire population...
        state$active <- state$active[-inds_rem]
        state$end_step <- state$end_step[-inds_rem]
        state$birth_step <- state$birth_step[-inds_rem]

    } else if(allow_extinction && len_rem > 0) {  # remove old individuals
        # If we would
        state$active <- state$active[-inds_rem]
        state$end_step <- state$end_step[-inds_rem]
        state$birth_step <- state$birth_step[-inds_rem]
    }

    return(state)
}


#' @export
check_safe_sizes <- function(state, parameters) {
    # make sure that we have enough space to insert the new individuals
    # if we don't, add on at enough space to do the smaller of 12 generations or total_steps generations
    if ((state$active_index + state$tot_offsprings - 1) > dim(state$parents)[1]) {
        # Double the table size. This should not need to occur, but will break the simulation
        # if we overflow the table. Default values are chosen with HIV dynamics in mind.
        state$parents <- rbind(state$parents,
            matrix(0, max(parameters$total_steps, 12) * parameters$maximum_population_target, 2))

        state$lifespan_types <- c(state$lifespan_types,
                                  rep(0, max(parameters$total_steps, 12) * parameters$maximum_population_target))

        state$spread_rate_types <- c(state$spread_rate_types,
                                     rep(0, max(parameters$total_steps, 12) * parameters$maximum_population_target))

        state$dispersion_rate_types <- c(state$dispersion_rate_types,
                                         rep(0, max(parameters$total_steps, 12) * parameters$maximum_population_target))
    }
    return(state)
}

#' @export
update_vital_statistics <- function(parents, state, parameters, ...) {
    # 1. Add new births to active
    # 2. Update birth step for new individuals
    # 3. Set the end step for new individuals. Uses the lifespan samplers

    # Update state vectors needed for stepping the model
    if(state$tot_offsprings > 0) {  # only attempt sampling if we have new individuals to get dates for
    state$active <- c(state$active, state$active_index:(state$active_index + state$tot_offsprings - 1))
    state$birth_step <- c(state$birth_step, rep(state$curr_step + 1, state$tot_offsprings))
        new_lifespans <- SEEPS::multi_sampler(
            birth_times = rep(state$curr_step + 1, state$tot_offsprings),
            current_time = state$curr_step,
            sampler_indexing = state$lifespan_types[state$active_index:(state$active_index + state$tot_offsprings - 1)],
            samplers = parameters$lifespan_sampler,
            params = list()
        )
        state$end_step <- c(state$end_step, state$curr_step + new_lifespans)
    }

    return(state)
}

#' @export
record_parents <- function(parents, state, parameters, ...) {
    # First, check if we have any parents to record
    if (state$tot_offsprings == 0 || length(parents) == 0) {
        # if not, just return the state
        return(state)
    }

    # We have some parents to assign
    if ((state$active_index + state$tot_offsprings - 1) > dim(state$parents)[1]) {
        # Double the table size. This should not need to occur, but will break the simulation
        # if we overflow the table. Default values are chosen with HIV dynamics in mind.
        state$parents <- rbind(
            state$parents,
            matrix(0, max(parameters$total_steps, 12) * parameters$maximum_population_target, 2)
        )
    }

    # Record parents

    state$parents[state$active_index:(state$active_index + state$tot_offsprings - 1), 1] <- state$active[parents]
    state$parents[state$active_index:(state$active_index + state$tot_offsprings - 1), 2] <- state$curr_step


    return(state)
}

#' @export
assign_inherit <- function(parents, state, parameters, ...) {
    # First, check if we have any offspring to assign properties to

    if (state$tot_offsprings == 0) {
        return(state)
    } # ELSE: we have offspring to assign properties to, let's proceed

    # Need to assign rates for:
    # 1. `lifespan_types`
    new_lifespan_indexes <- parameters$lifespan_selector(
        state$birth_step,
        current_time = state$curr_step,
        rate_indexing = state$lifespan_types[state$active, ],
        rate_functions = parameters$lifespan_sampler,
        parents = parents
    )

    state$lifespan_types[state$active_index:(state$active_index + state$tot_offsprings - 1)] <-
        new_lifespan_indexes

    # 2. `spread_rate_types`

    new_rate_indexes <- parameters$rate_selector_function(
        state$birth_step,
        current_time = state$curr_step,
        rate_indexing = state$spread_rate_types[state$active, ],
        rate_functions = parameters$offspring_rate_fn,
        parents = parents)

    state$spread_rate_types[state$active_index:(state$active_index + state$tot_offsprings - 1)] <-
        new_rate_indexes

    # 3. `dispersion_rate_types`

    new_dispersion_indexes <- parameters$dispersion_selector(
        state$birth_step,
        current_time = state$curr_step,
        rate_indexing = state$dispersion_rate_types[state$active, ],
        rate_functions = parameters$dispersion_sampler,
        parents = parents
    )


    state$dispersion_rate_types[state$active_index:(state$active_index + state$tot_offsprings - 1)] <-
        new_dispersion_indexes

    return(state)
}


#' @export
compute_dispersion <- function(offsprings, state, parameters, ...) {
    if ("dispersion_selector" %in% names(parameters)) {
        offsprings <- SEEPS::dispersion_calculator(
            mean = offsprings,
            current_step = state$curr_step,
            birth_step = state$birth_step,
            dispersion_type = state$dispersion_rate_types[state$active],
            samplers = parameters$dispersion_sampler,
            params = list()
        )
    }
    return(offsprings)
}


#' @export
compute_offsprings <- function(state, parameters, ...) {
    # offsprings <- parameters$offspring_rate_fn(current_step = state$curr_step, birth_step = state$birth_step)
    offsprings <- SEEPS::rate_calculator(
        birth_times = state$birth_step,
        current_time = state$curr_step,
        rate_indexing = state$spread_rate_types[state$active],
        rate_functions = parameters$offspring_rate_fn,
        params = list()
    )
    return(offsprings)
}


# Candidate version of the code
#' Perform a single step of the simulation
#'
#' The workhorse function of SEEPS. Step forward a state one time step.
#'
#' @param state The state object (list) from the simulation. See `initialize` for the expected structure.
#' @param parameters A list of parameters to run the simulation. See `wrap_parameters` for the expected structure.
#' @export
step2 <- function(state, parameters) {
    # Perform a single step of the simulation. This is the core algorithm. It is used
    # by both the exponential and constant phases which provide hyperparameter inputs
    # This is called repeatedly until the simulation is complete. This represents a single step.

    state$length_active <- length(state$active) # pre-compute for efficiency
    if (parameters$maximum_population_target - state$length_active > 0) { # If we can/need to grow
        # Sample offspring distribution
        rate_fn <- parameters$offspring_rate_fn[[1]]
        offsprings <- rate_fn(current_step = state$curr_step, birth_step = state$birth_step)

        # offsprings <- parameters$offspring_rate_fn(current_step = state$curr_step, birth_step = state$birth_step)
        # Cap offspring at the target
        state$total_offsprings <- ceiling(min(
            parameters$maximum_population_target - state$length_active,
            sum(offsprings)
        ))
        # Assign parents using the sampled rates/weights
        prts <- sample(1:state$length_active, state$tot_offsprings, replace = TRUE, prob = offsprings)
        # Record parents
        if ((state$active_index + state$total_offsprings - 1) > dim(state$parents)[1]) {
            # Double the table size. This should not need to occur, but will break the simulation
            # if we overflow the table. Default values are chosen with HIV dynamics in mind.
            state$parents <- rbind(state$parents,
                                   matrix(0, max(parameters$total_steps, 12) * parameters$maximum_population_target, 2))
        }
        state$parents[state$active_index:(state$active_index + state$tot_offsprings - 1), 1] <- state$active[prts]
        state$parents[state$active_index:(state$active_index + state$tot_offsprings - 1), 2] <- state$curr_step

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
#' Remove a set of samples from the simulation. Intended to be used with contact tracing, we expect
#' individuals identified through contact tracing to be removed from the population of active
#' (uncontrolled or undiagnosed) cases.
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
    if (!all(c("parents", "active", "active_index", "end_step", "birth_step", "length_active",
               "curr_step", "total_offspring") %in% names(state))) {
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
        warning(paste0("End step must only contain integers. This is not enforced by the type system. End step is ",
                       state$end_step))
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
    if (length(names(state)) < 8) {
        warning("State has the wrong number of values. This is not enforced by the type system.")
    }
}

is_integer <- function(x) {
    all(x == floor(x))
}
