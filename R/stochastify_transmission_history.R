#' Sample evolutionary distances to an integer number of mutations per branch
#'
#' Convert temporal distances for a geneology to an integer number of mutations.
#' This step is extracted so a single geneology can be re-sampled and create multiple possible matrices.
#'
#' @param transmission_history A matrix output by `reduce transmission history`.
#' @param rate Clock rate for evolutionary process with units of substitions per sequence per year.
#' @export
stochastify_transmission_history <- function(transmission_history, rate) {
    inds <- 1:(dim(transmission_history)[1] - 2)
    transmission_history[inds, 5] <- rpois(length(inds), transmission_history[inds, 4] * rate)
    # taking lambd =  here is not a problem - R describes it as a point mass at 0
    return(list("geneology" = transmission_history))
}
