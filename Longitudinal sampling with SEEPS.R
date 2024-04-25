## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup--------------------------------------------------------------------
library(SEEPS)
set.seed(1945)


## -----------------------------------------------------------------------------
states_sampled <- list()
samples_all <- list()

samples_per_year <- 5
# Discover each contact with 80% probability
contact_tracing_parameter <- 0.8

sim_parameters <- SEEPS::wrap_parameters(
    minimum_population = 6,
    offspring_rate_fn = SEEPS::get_biphasic_HIV_rate(params = list("R0" = 3)),
    maximum_population_target = 20,
    total_steps = 0,  # We are going to manually control this so the value will not be used
    spike_root = FALSE)

state <- SEEPS::initialize(sim_parameters)  # Setup
state <- SEEPS::gen_exp_phase(state, sim_parameters)  # Run until we have at least 90% of the maximum population

# Take a first sample at the end of the exponential phase
sample_taken <- SEEPS::contact_traced_uniform_restarts_ids(
    state$active, state$parents, samples_per_year, contact_tracing_parameter)
# Record the sampling
samples_all[[1]] <- list("curr_step" = state$curr_step, "samples" = sample_taken$samples)
# Remove from the population
state <- SEEPS::remove_samples(state, sample_taken$samples)

# Now simulate forward
for (time_pointer in 1:3) {
    state <- SEEPS::gen_const_phase(state, sim_parameters, num_steps = 12)

    sample_taken <- SEEPS::contact_traced_uniform_restarts_ids(
        state$active, state$parents, samples_per_year, contact_tracing_parameter)
    # Record the sampling
    samples_all[[time_pointer + 1]] <- list("curr_step" = state$curr_step, "samples" = sample_taken$samples)
    # Remove from the population
    state <- SEEPS::remove_samples(state, sample_taken$samples)
}

# Now that we have collected all of the sampling data, we can flatten it into a single list

sample_times_ <- list()
sample_ids <- list()
ptr <- 1
for (i in seq_along(samples_all)) {
    sample_times_[[ ptr]] <- samples_all[[i]]$curr_step
    sample_ids[[ptr]] <- samples_all[[i]]$samples
    ptr <- ptr + 1
}

# Convert to

genealogy <- SEEPS::reduce_transmission_history_bpb2(
    samples = sample_ids,
    parents = state$parents,
    current_step = sample_times_)

# Convert to a phylogeny

phylogeny <- SEEPS::geneology_to_phylogeny_bpb(
    transmission_history = genealogy$parents,
    infection_times = genealogy$transmission_times,
    sample_times = genealogy$sample_times,
    a = 5, b = 5,
    leaf_sample_ids = genealogy$transformed_sample_indices)

# Convert to a newick string for plotting
newick_tree_string <- SEEPS::phylogeny_to_newick(phylogeny$phylogeny, mode = "mean", label_mode = "abs")

# Plot the tree
library(ape)
tree <- drop.tip(read.tree(text = newick_tree_string), "")
plot(tree, cex = 0.5)
axis(1, las=1, col.axis="black", col.ticks="black")


## -----------------------------------------------------------------------------
gen_transmission_history_exponential_constant_local <- function(minimum_population, #nolint: object_length_linter
        offspring_rate_fn, maximum_population_target, total_steps,
                              spike_root = FALSE) {
    # Top level/user function. Users should start with this function.

    parameters <- SEEPS::wrap_parameters(minimum_population, offspring_rate_fn,
                                  maximum_population_target, total_steps,
                                  spike_root)

    state <- SEEPS::initialize(parameters)
    SEEPS::validate_state(state)  # Check that the state is valid
    all_states <- list()
    all_states[[1]] <- state
    state <- gen_exp_phase(state, parameters, all_states = TRUE)
    all_states <- c(all_states, state)

    state <- gen_const_phase(all_states[[length(all_states)]], parameters, num_steps = total_steps, all_states = TRUE)
    all_states <- c(all_states, state)
    return(lapply(all_states, SEEPS::clean_up))
}


## -----------------------------------------------------------------------------
offspring_rate_fn <- get_biphasic_HIV_rate_function(50, 3, 24, params = list("R0" = 1.1))


## -----------------------------------------------------------------------------
maximum_population_target <- 30
simulator_result <- gen_transmission_history_exponential_constant_local(
        minimum_population = 6,
        offspring_rate_fn = offspring_rate_fn,
        # offspring_rate_fn = SEEPS::get_biphasic_HIV_rate(params = list("R0" = 1.3)),
        maximum_population_target = maximum_population_target,
        total_steps = 48)

removed_list <- list()
sample_times <- list()
birth_list <- list()


## -----------------------------------------------------------------------------
ptr <- 1
removed_list[[ptr]] <- c()
sample_times[[ptr]] <- c()
birth_list[[ptr]] <- c(1)
ptr <- 2


## -----------------------------------------------------------------------------

for (i in 2:length(simulator_result)) {
    # Get the active nodes at i-1 and i
    active_i <- simulator_result[[i]]$active
    active_i_1 <- simulator_result[[i-1]]$active

    # Compare the active nodes at i-1 and i
    # If the active node is in i-1 but not i, add it to a "removed" list at position i
    removed_active <- active_i_1[!(active_i_1 %in% active_i)]
    new_active <- active_i[!(active_i %in% active_i_1)]
    if (length(removed_active)  == 0 && length(new_active) == 0) {
        next  # Nothing to add to the samples list
    }

    # Add to removed list
    removed_list[[ptr]] <- removed_active
    birth_list[[ptr]] <- new_active
    # Record the sample time
    sample_times[[ptr]] <- simulator_result[[i]]$t_end
    ptr <- ptr + 1
}

# Final step - get the active nodes at the last step
active_final <- simulator_result[[length(simulator_result)]]$active
# Store it
removed_list[[ptr]] <- active_final
# Record the sample time
sample_times[[ptr]] <- simulator_result[[length(simulator_result)]]$t_end


## -----------------------------------------------------------------------------
# Now record the transmission history
# print(removed_list)
# print(sample_times)
full_transmission_history <- reduce_transmission_history_mt(parents=simulator_result[[length(simulator_result)]]$parents,
    samples=removed_list,
    current_step=sample_times)
# ```
#
# Now we can visualize the entire transmission history.
#
# ```{r, fig.height = 12}
# string <- SEEPS::phylogeny_to_newick(full_transmission_history$geneology, mode="mean", label_mode="abs")
# tree <- read.tree(text = string)
#
# plot(ladderize(tree))
# axis(1, las = 1, col.axis = "black", col.ticks = "black")
# ```
#
# We can also plot out the number of new infections per time step.
#
# ```{r}
# births_ <- lapply(birth_list, length)
# deaths <- lapply(removed_list, length)
# # births_ <- births[1:length(births) - 1]  # Don't report the final step
# deaths_ <- deaths[1:length(deaths) - 1]  # Don't report the final step
# sample_times_ <- sample_times[1:length(sample_times) - 1]  # Don't report the final step
# deaths_ <- unlist(deaths_)
# births_ <- unlist(births_)
# print(length(births_))
# print(length(deaths_))
# plot(seq_along(births_), births_, label = "Births")
# lines(seq_along(deaths_), deaths_, col="red", label="Deaths")
# legend(39,y=NULL, legend=c("Births", "Deaths"), col=c("black", "red"), lty=1:2, cex=0.8)
# ```
#
# ```{r}
# print(length(births_))
# print(length(deaths_))
# ```
#
# ```{r}
# deltas <- births_ - deaths_
# print(length(deltas))
# plot(1:length(deltas), cumsum(deltas), ylab="Number of active infections", xlab="Time index [steps]", ylim=c(-1, maximum_population_target + 2))
# # Draw a horizontal line at y=20 with a dashed line
# abline(h = maximum_population_target, col = "red", lty = 2)
