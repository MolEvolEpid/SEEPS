% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contact_tracing.R
\name{uniform_discovery_factory}
\alias{uniform_discovery_factory}
\title{A factory function to discover connections with uniform probability}
\usage{
uniform_discovery_factory(p)
}
\arguments{
\item{p}{The discovery probability. A float between 0 and 1.}
}
\description{
Obtain a discovery function that determines which adjacent nodes in the
contact network will be revealed with uniform probability.
}
\seealso{
contact_traced_uniform_ids
}
