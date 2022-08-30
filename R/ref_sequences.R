# Interact with builtin reference sequences.
# See also: /data-raw/load_sequences.R for construction of the data object

lookup_sequence_by_index <- function(organism_name, start, stop) {
    # Get the sequence
    sequence <- pkgenv$seq_store[[organism_name]]$sequence
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

lookup_sequence_by_name <- function(organism_name, region_name) {
    # Get the organism
    organism <- pkgenv$seq_store[[organism_name]]
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
