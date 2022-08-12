
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bluefern

<!-- badges: start -->
<!-- badges: end -->

A modern and modular simulator for phylogenetics and phylodynamics.

<!-- One short paragraph about simulation -->

## Installation

You can install the released version of bluefern from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("bluefern")
```

## Example

Obtain pairwise distance matrices with as little as 3 lines of code:

``` r
install.packages("bluefern")
parameters <- list("rate_function_parameters" = list("R0" = 5),
    "minimum_population_size" = 15, "maximum_population_target" = 100)
simultion_results <- bluefern::simulate_classic_HIV(parameters)
```

## Modularity

Bluefern offers complete modularity in designing simulations without an
additional class structure. Simulation steps are interchangable
functions. You can insert, modify, or remove workflow steps at ease. The
`simulate_classic_HIV` function we used above expands to:

``` r
simulate_classic_HIV <- function(params) {
    # This gives an example "workflow" function for this package.
    # Simple integration tests should use similar workflows to insert changes.

    # Step 1: Determine the rate function for offspring generation
    # This obtains the function used in [Graw et al. 2012], [Kupperman et al. 2022]
    biphasic_rate_function <- get_biphasic_HIV_rate(
        params = params[["rate_function_parameters"]])

    # Step 2: Forward pass to generate transmission history
    simulator_result <- gen_transmission_history_exponential_constant(
        minimum_population=params[["minimum_population"]],
        offspring_generation_rate_function = biphasic_rate_function,
        NPop = params[["maximum_population_target"]],
        spike_root = FALSE)
        # spike root option allows you to record distances to the initial infection

    # Step 3: determine the sample of ID's to subsample.
    # This function takes independent random samples, but you could do
    # something much more fancy here and change this call signature
    target_sample <- random_ids(
        active = simulator_result[["active"]],
        minimum_size = params[["minimum_population"]],
        proportional = FALSE,  # sample a proportion of active infections
        spike_root = FALSE)

    # Step 4: Obtain the transmission history for the subset of sampled individuals
    transmission_history <- reduce_transmission_history(
        active = simulator_result[["active"]],
        parents = simulator_result[["parents"]],
        observation_size = target_sample[["observation_size"]],
        current_step = simulator_result[["t_end"]],
        tot_offsprings = simulator_result[["total_offspring"]],
        spike_root = FALSE)

    #  Step 5: Use the geneology simulation method in [Authors et al. Y###]
    # This makes synthetic data look more "realistic"
    geneology <- transmissions_to_geneology(
        join_tree = transmission_history[["ancestry tree"]],
        internode_distances = transmission_history[["branch_data"]],
        spike_root=FALSE)

    # Step 6: Reduce the geneology to a pairwise distance matrix
    distance_matrix <- geneology_to_distance_matrix(
        geneology=geneology[["geneology"]],
        spike_root=FALSE)

    # Now return the result to the user. Additional metadata used above
    # can be collected and added to the returned list.
    return(list("matrix" = distance_matrix, "input_params" = params))
}
```

Showing 6 steps in this simulation:

1.  Obtain a rate function for offspring generation. We use a simple 2
    phase rate function `get_biphasic_HIV_rate`, to model an increased
    offspring generation rate following initial infection, then a low
    rate after.
2.  Simulate the transmission network. This simulation requires an
    offspring generation rate function (step 1), and min/max population
    sizes.
3.  Obtain a random sampling of active infections that are sampled. It
    is computationally expensive to reconstruct large trees.
4.  Reduce the simulation output to a transmission history tree.
5.  Convert the transmission history tree to a phylogeny using \[Author
    et al, Y###\] to obtain more realistic results
6.  Reduce the phylogeny to a pairwise distance matrix.

## Additional options and features

-   Tree subsampling. Between steps 5 and 6 above, further drop
    individuals using knowledge of the reconstructed tree.
-   We provide support for a broad class of discretized rate functions,
    the most general being `get_Kphasic_hiv_rate_function`. This
    function accepts a list of rates for each simulation step and the
    length of each phase. <!-- This needs an example -->
-   Cluster-driven sampling. An expensive way to obtain a cluster is to
    recover the history for a large proportion of the population, then
    take only the closest `k-1` indiivduals around a randomly sampled
    individual. The result is a sample of `k` closely related
    individuals, but with the cost of reconstructing `n > k` individuals
    in a larger tree.
-   Root-to-tip distances. Changing `spike_root` to`TRUE` above will add
    and track the root infection through the framework. The last
    row/column in the distance matrix will correspond to the distance to
    the root infection. For large simulations, this may increase
    computational time as the tree will be rooted at the index case,
    rather than the most recent common ancestor (MRCA).

## Todo list

In no particular order.

-   Fast subsampling options for obtaining clusters without needing to
    sample a large proportion and down-sample the tree. This approach is
    very very expensive for over 500 sequences.
-   Calculus utilities to convert almost any function into a rate
    function (numerical integration).
-   Support joining clusters together from individual simulations to
    form large distance matrices.
