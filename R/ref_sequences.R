# Interact with builtin reference sequences.
# See also: /R/load_sequences.R for construction of the data
#' Obtain a sequence from the builtin reference sequences using coordinates
#'
#' Provide an interval and a reference sequence name.
#' Currently supported are: "HIV1" (HXB2 reference) and "toy" (a poly A/C/G/T
#' example sequence for testing).
#'
#' @param start The start coordinate of the sequence
#' @param end The end coordinate of the sequence
lookup_sequence_by_index <- function(organism_name, start, stop) {
    # Get the sequence
    sequence <- pkgenv$seq_store[[organism_name]]$sequence  # nolint: object_usage_linter
    # Get the sequence length
    sequence_length <- length(sequence)
    # Check the start and stop are valid
    if (start < 1 || start > sequence_length || stop < 1 || stop > sequence_length) {
        msg <- paste("The start and stop positions are invalid.", "start is ", start,
                     "stop is ", stop, ". The querried length is ", sequence_length)
        stop(msg)
    }  # else:
    # Return the sequence
    return(sequence[start:stop])
}
#' Obtain a sequence from the builtin reference sequences using an name and annotated region
#'
#' Return a portion of a reference genome using a standard annotation name.
#'
#' Currently supported are: "HIV1" (HXB2 reference) and "toy" (a poly A/C/G/T
#' example sequence for testing). Supported HIV1 regions include the short annotated
#' list, see (here)[https://www.hiv.lanl.gov/components/sequence/HIV/search/help.html#region]
#' for the full list. Clinical regions (p17, V3) are supported as "p17-clinical" and
#' "v3-clinical".
#'
#' For a more detailed lookup proceedure with the reference sequences using user-provided
#' coordinates, see `lookup_sequence_by_index`.
#'
#' @seealso lookup_sequence_by_index
#' @export
#'
lookup_sequence_by_name <- function(organism_name, region_name) {
    # Get the organism
    organism <- pkgenv$seq_store[[organism_name]]  # nolint: object_usage_linter
    # Get the region
    region_map <- organism$region
    # Find the correct region
    ids <- which(region_map$region == region_name)
    # Filter
    region <- region_map[ids, ]
    # Get the sequence
    sequence <- organism$sequence[region$start:region$end]
    # Return the sequence
    return(sequence)
}

#' Show the list of available reference sequences
#'
#' A list of all available organism/sequence names is printed to the screen.
#' @export
show_available_sequences <- function() {
    # Get the organism names
    organism_names <- names(pkgenv$seq_store)  # nolint: object_usage_linter
    # Print the organism names
    cat(organism_names, "\n")
    return(NULL)
}

#' Show the list of available regions for a provided reference sequence
#'
#' @param organism_name The name of the reference sequence
#' @param chart bool (default TRUE). If TRUE, a chart of the available regions
#' is printed to the screen.
#' @return A list of available regions for the named reference sequence.
#' If the organism_name is not found, an NULL is printed.
#' @export
show_available_regions <- function(organism_name, chart = TRUE) {
    # Get the organism
    organism <- pkgenv$seq_store[[organism_name]]  # nolint: object_usage_linter
    # Get the region map
    region_map <- organism$region
    # Get the region names
    if (chart) {
        output <- region_map
    } else {
        output <- region_map$region
    }
    print(output)
}
