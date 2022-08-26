check_2line_fasta <- function(fasta_string) {
    # Check that the fasta file is a 2-ine fasta file
    # iterate through pairs of rows and check that the first line
    # is the name of the sequence and the second line is the sequence
    lines <- unlist(strsplit(fasta_string, split = "\n"))
    line_counter <- 1
    for (line in lines[1:length(lines)]) { # nolint
        if (line_counter %% 2 == 0) {
            # Sequence line
            # Check that all chatacters are ACTG
            if (!all(grepl("[ACGT]", line))) {
                stop("invalid character detected")
            }
        } else {
            # Check that the first character is a >
            if (substring(line, 1, 1) != ">") {
                message("Not a valid definition line")
            }
        }
        line_counter <- line_counter + 1
    }
    # If you made it this far, return TRUE
    return(TRUE)
}
