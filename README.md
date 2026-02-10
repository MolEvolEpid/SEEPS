
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SEEPS: Sequence Evolution and Epidemiological Process Simulator

<!-- badges: start -->

<!-- badges: end -->

A modern and modular simulator for phylogenetics and phylodynamics
written in pure R, SEEPS provides advanced simulation capabilities for
population dynamics, sampling techniques, genealogy and phylogeny
extraction, can generate trees, sequences, or pairwise distance
matrices.

Check out the
[tutorial](https://mol-evo-epid.github.io/SEEPS/articles/SEEPS.html) for
more details.

If you use SEEPS in your research, please cite the following paper:

> Kupperman MD, Ke R, Leitner T. 2025. Identifying Impacts of Contact Tracing on HIV Epidemiological Inference from Phylogenetic Data. Virus Evolution 11(1): veaf068.
>  https://doi.org/10.1093/ve/veaf068

## Installation

You can install the released version of SEEPS from
[R-universe](https://r-universe.dev/) with:

``` r
install.packages("SEEPS", repos = "https://r-universe.dev")
```

Or install the development version from
[GitHub](https://github.com/MolEvolEpid) with:

``` r
devtools::install_github("MolEvolEpid/SEEPS")
```

## Example

Obtain pairwise distance matrices with as little as 3 lines of code:

``` r
install.packages("SEEPS", repos = "https://r-universe.dev")
parameters <- list(
    "rate_function_parameters" = list("R0" = 5),
    "minimum_population" = 15, "maximum_population_target" = 1000,
    "total_steps_after_exp_phase" = 0, "mutation_rate" = 0.0067 * 300
)
simultion_results <- SEEPS::simulate_classic_HIV(parameters)
```

## Modularity

SEEPS offers complete modularity in designing simulations without an
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
        params = params[["rate_function_parameters"]]
    )

    # Forward pass to generate transmission history
    simulator_result <- gen_transmission_history_exponential_constant(
        minimum_population = params[["minimum_population"]],
        offspring_rate_fn = biphasic_rate_function,
        total_steps = params[["total_steps_after_exp_phase"]],
        maximum_population_target = params[["maximum_population_target"]],
        spike_root = FALSE
    )
    # spike root option allows you to record distances to the initial infection

    # Next, determine the sample of ID's to subsample.
    # This function takes independent random samples, but you could do
    # something much more fancy here and change this call signature
    target_sample <- random_fixed_size_ids(
        active = simulator_result[["active"]],
        minimum_size = params[["minimum_population"]],
        spike_root = FALSE
    )

    # Obtain the transmission history for the subset of sampled individuals
    transmission_history <- reduce_transmission_history(
        samples = target_sample[["samples"]],
        parents = simulator_result[["parents"]],
        current_step = simulator_result[["t_end"]],
        spike_root = FALSE
    )

    # Convert time signals to # of mutations using a rate
    geneology <- stochastify_transmission_history(
        transmission_history = transmission_history$geneology,
        rate = params[["mutation_rate"]] / 12
    ) # Provide in rate per sequence per year

    # convert the geneology into a distance matrix
    distance_matrix <- geneology_to_distance_matrix(
        geneology = geneology$geneology,
        spike_root = FALSE
    )

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
5.  Convert the transmission history tree to a phylogeny to obtain more
    realistic results.
6.  Reduce the phylogeny to a pairwise distance matrix.

## Additional options and features

- Generate sequences with
  [`seq-gen`](http://tree.bio.ed.ac.uk/software/seqgen/) \[Rambaut and
  Grassly\]. Provide a phylogeny matrix and a root sequence to
  `generate_sequences`. Pre-built rate models and reference sequences
  are available. See the tutorial for more details.
  <!-- todo: write this tutorial -->
- Contact tracing. Use knowledge of the contact network to determine a
  sample of individuals, rather than take a random sample.
- Tree subsampling. Between steps 5 and 6 above, further drop
  individuals using knowledge of the reconstructed tree.
- We provide support for a broad class of discretized rate functions,
  the most general being `get_Kphasic_hiv_rate_function`. This function
  accepts a list of rates for each simulation step and the length of
  each phase. <!-- This needs an example -->
- Cluster-driven sampling with variable intensity. An approximation of a
  hybrid between contact tracing and random sampling is to randomly
  sample a proportion (a *sampling intensity* between 1% and 100%) of
  the active infections, then select a small subset of closely related
  sequences. For pathogens (SARS-CoV2) with little or no variation
  across most transmission events, a lower density will produce more
  diverse samples, while high sampling intensity will result in closely
  related or many identical sequences being obtained. This functionality
  is provided through the `proportional` option in `random_ids` and the
  `reduce_large_matrix` functions.
- Root-to-tip distances. Changing `spike_root` to`TRUE` above will add
  and track the root infection through the framework. The last
  row/column in the distance matrix will correspond to the distance to
  the root infection. For large simulations, this may increase
  computational time as the tree will be rooted at the index case,
  rather than the most recent common ancestor (MRCA).

## Todo list

In no particular order. Please reach out or open a pull request if you
would like to contribute to any of these features:

- Support for multiple types within the agent-based simulation
- Support for individual level contact tracing performance (give each
  individual a discovery rate)
- Temporal contact tracing, where contact tracing is a dynamic process
  through time. This would be a more realistic model of contact tracing,
  since it is not instantaneous.
- Support for joining simulations together. This would form the basis of
  a metapopulation model, and enable primative multithreading/multicore
  simulations.

Think we’re missing a core feature? Open an issue and let us know!

## License

SEEPS is licensed under a modified BSD-3 license.
