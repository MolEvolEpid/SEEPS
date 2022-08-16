
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
    # Simple integration tests should use similar workflows to test changes.

    # First, determine the rate function for offspring generation
    # This obtains the function used in [Graw et al. 2012], [Kupperman et al. 2022]
    biphasic_rate_function <- get_biphasic_HIV_rate(
        params = params[["rate_function_parameters"]])

    # Forward pass to generate transmission history
    simulator_result <- gen_transmission_history_exponential_constant(
        minimum_population=params[["minimum_population"]],
        offspring_rate_fn = biphasic_rate_function,
        total_steps = params[["total_steps_after_exp_phase"]],
        maximum_population_target = params[["maximum_population_target"]],
        spike_root = FALSE)
        # spike root option allows you to record distances to the initial infection

    # Next, determine the sample of ID's to subsample.
    # This function takes independent random samples, but you could do
    # something much more fancy here and change this call signature
    target_sample <- random_ids(
        active = simulator_result[["active"]],
        minimum_size = params[["minimum_population"]],
        proportional = FALSE,  # sample a proportion of active infections
        spike_root = FALSE)

    # Obtain the transmission history for the subset of sampled individuals
    transmission_history <- reduce_transmission_history(
        samples = target_sample[["samples"]],
        parents = simulator_result[["parents"]],
        current_step = simulator_result[["t_end"]],
        spike_root = FALSE)

    # Convert time signals to # of mutations using a rate
    geneology <- stochastify_transmission_history(
        transmission_history=transmission_history$geneology,
        rate=params[["mutation_rate"]] / 12)  # Provide in rate per sequence per year

    # convert the geneology into a distance matrix
    distance_matrix <- geneology_to_distance_matrix(
        geneology=geneology$geneology,
        spike_root=FALSE)

    # Now return the result. Additional metadata used above can be collected
    # and added to the return list.
    return(list("matrix" = distance_matrix))
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
