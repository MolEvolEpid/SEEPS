% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reduce_transmission_history_to_true_geneology.R
\name{reduce_transmission_history_bpb}
\alias{reduce_transmission_history_bpb}
\title{Reduce simulation output to transmission history for a subset to include
  within host diversity.}
\usage{
reduce_transmission_history_bpb(samples, parents, current_step)
}
\arguments{
\item{samples}{A vector of individuals (integers) to include in the sample.}

\item{parents}{A matrix of parental individuals that encodes the transmission
history and sample times.}

\item{current_step}{The current (absolute) time step in the simulation.}

\item{spike_root}{A boolean indicating whether the geneology should
 include the root of the outbreak or not. Default is `FALSE`.
 This should be specified even if the founding infection is sampled,
as the root of the outbreak will have evolved since the founding event.}
}
\value{
A list with 4 elements:
`parents` A vector of parents of each infection in the sample until the root
`times` A vector of times of sampling times in the tree. Sample times for
internal nodes are after the last offspring generation time needed to
reconstruct the sample.
`transmission_times` A vector of transmission times of each infection in the tree.
`samples_available` A boolean vector (mask) of which samples are leaves in the tree.
   Used by the coalescent simulation to know which individuals should be assigned
   detected sequences.
}
\description{
For a detailed explenation of inputs, see `reduce_transmission_history`, which
is intended to reconstruct back only until the most recent common ancestor
 of the sample, and return a tree.
}
\details{
To include within host diversity, use `reduce_transmission_history_bpb` to
extract a subset of the transmission history that we need for the phylogeny,
and see `geneology_to_phylogeny_bpb` to simulate within-host diversity
and recover a true phylogeny.
}
\seealso{
reduce_transmission_history_bpb
}
