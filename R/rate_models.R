# Rate models for evolution

###############################################################################
#
# Rate models from Swedish transmission chain. See
# [Leitner et al. (1997)](https://doi.org/10.1128/JVI.71.6.4761-4770.1997)
# for details. For parameters not reported in the paper, the values are
# obtained using the same calculations as described in the paper.
#
###############################################################################


#' Rate model for env V3 estimated from the Swedish transmission chain
#'
#' Estimates for GTR+I+G rate model parameter from [Leitner et al. 1997] on the
#' Swedish transmission chain on envelope V3 region.
#'
#' @return A list of rate model parameters.
#' @seealso generate_sequences generate_rate_model
#' @examples
#' rate_model <- get_V3_rate_model()
#'
#' @export
get_V3_rate_model <- function() { # nolint: object_name_linter
    rate_model <- list(
        a2c = 0.594, a2g = 1, # parameter f = 1 is the reference
        a2t = 0.186, c2g = 0.125,
        c2t = 0.822, g2t = 0.235,
        # Estimates in Leitner et al. 1997 (Table 1) do not
        # sum to 1 due to rounding, so we normalize the estimate.
        # Confidence intervals for the estimate are >5% of the
        # estimate we use here.
        fa = 0.4627 / 1.001, fc = 0.1474 / 1.001,
        fg = 0.1598 / 1.001, ft = 0.2302 / 1.001,
        # i is observed fraction. With a larger sample,
        # this will be too low.
        i = 0.68,
        # alpha is the shape parameter for the discrete gamma
        # ncat is the number of categories
        alpha = 0.384, ncat = 8
    )
    return(rate_model)
}


#' Rate model for gag p17 estimated from the Swedish transmission chain
#'
#' Estimates for GTR+I+G rate model parameter from [Leitner et al. 1997] on the
#' Swedish transmission chain on p17. Rate estimates for p17 follow from a
#' secondary analysis.
#'
#' @return A list of rate model parameters.
#' @seealso generate_sequences generate_rate_model
#' @examples
#' rate_model <- get_p17_rate_model()
#' @export
get_p17_rate_model <- function() {
    stop("This model is not yet implemented")
    rate_model <- list(
        a2c = -1, a2g = 1, # Need to get numbers here
        a2t = -1, c2g = -1,
        c2t = -1, g2t = -1,
        fa = 0.3995, fc = 0.1677,
        fg = 0.2487, ft = 0.1841,
        # i is observed fraction. With a larger sample,
        # this will be too low.
        i = 0.85,
        # alpha is the shape parameter for the discrete gamma
        # ncat is the number of categories
        alpha = 0.257, ncat = 8
    )
    return(rate_model)
}

###############################################################################
#
# Utility functions and helper functions for rate models
#
###############################################################################

#' Generate a rate model from provided paramters
#'
#' Accept a set of parameters for GTR+I+G model. See below for details of the
#' parameterization. To disable I, set the parameter i to 0. To disable G, set
#' alpha=1 and ncat=1 for an exponential distribution.
#'
#' For an overview of GTR+I+G and substitition models, see
#' [here](https://www.ccg.unam.mx/~vinuesa/Model_fitting_in_phylogenetics.html).
#' For a more detailed construction, see
#' [Yang 1994](https://doi.org/10.1007/BF00178256), and
#' [Yang 1996](https://doi.org/10.1093/oxfordjournals.molbev.a025625).
#'
#'
#' @param a2c The rate of adenosine to cytosine transversions.
#' @param a2g The rate of adenosine to cytosine transitions.
#' @param a2t The rate of adenosine to thymine transversions.
#' @param c2g The rate of cytosine to guanine transversions.
#' @param c2t The rate of cytosine to thymine transitions.
#' @param g2t The rate of guanine to thymine transversions.
#' @param fa The fraction of adenosine at equilibrum.
#' @param fc The fraction of cytosine at equilibrum.
#' @param fg The fraction of guanine at equilibrum.
#' @param ft The fraction of thymine at equilibrum.
#' @param i  The proportion of invariant sites. Traditionally estimated by the
#'   proportion of sites with no observed mutations.
#' @param alpha The shape parameter for the discrete gamma distribution.
#' @param ncat The number of categories in the discrete gamma distribution.
#' @return A list of rate model parameters.
#'
#' @seealso  generate_sequences get_V3_rate_model get_p17_rate_model check_rate_model
#' @export
#' @examples
#' rate_model <- generate_rate_model(
#'     a2c = 1, a2g = 1, a2t = 1,
#'     c2g = 1, c2t = 1, g2t = 2,
#'     fa = 0.25, fc = 0.25,
#'     fg = 0.25, ft = 0.25,
#'     i = 0.1, alpha = 0.25, ncat = 8
#' )
#'
generate_rate_model <- function(a2c, a2g, a2t, c2g, c2t, g2t, fa, fc, fg, ft, i, alpha, ncat) {
    rate_model <- list(
        a2c = a2c, a2g = a2g, a2t = a2t, c2g = c2g, c2t = c2t, g2t = g2t,
        fa = fa, fc = fc, fg = fg, ft = ft,
        i = i,
        alpha = alpha, ncat = ncat
    )
    return(rate_model)
}
#' Check that a rate model is mathematically valid
#'
#' A helper function to check that a rate model paramterization is
#' mathematically valid.
#'
#' Check that all rates are positive, that the fraction of bases add up to 1,
#' the fraction of invariant sites is in [0, 1), and the discrete gamma model
#' is defined for `alpha` and `ncat`.
#'
#' @param rate_model A list of rate model parameters.
#' @seealso generate_rate_model generate_sequences
check_rate_model <- function(rate_model) {
    if (rate_model$a2c < 0 || rate_model$a2g < 0 || rate_model$a2t < 0 ||
        rate_model$c2g < 0 || rate_model$c2t < 0 || rate_model$g2t < 0 ||
        rate_model$fa < 0 || rate_model$fc < 0 || rate_model$fg < 0 ||
        rate_model$ft < 0 || rate_model$i < 0 || rate_model$alpha < 0 ||
        rate_model$ncat < 0) {
        stop("All rate model parameters must be positive")
    }
    if (rate_model$fa + rate_model$fc + rate_model$fg + rate_model$ft != 1) {
        stop("Fractions of each base must sum to 1")
    }
    if (rate_model$i >= 1) {
        stop("i must be strictly less than 1")
    }
    if (rate_model$ncat %% 1 != 0) {
        stop("ncat must be an integer")
    }
    # If you made it this far, the rate model is valid and return TRUE
    return(TRUE)
}
