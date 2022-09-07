#' Generate sequences from a phylogeny using Seq-Gen through phyclust
#'
#' Given a phylogeny, generate possible sequences using Seq-Gen [Rambaut & Grassley, 1997].
#' A root sequence for the simulation is required. The root sequence is placed
#' at the MRCA of the phylogeny unless `spike_root = TRUE` is specified when the phylogeny
#' is constructed.
#'
#' If the root sequence is a character vector, it is flattened into a single string.

#' @param phylogeny A phylogeny object
#' @param root_sequence A root sequence
#' @param rate_per_nt Flag to indicate whether the mutation rate is
#'   per nucleotide (nt) or per sequence. If per-sequence (`rate_per_nt=FALSE`),
#' the mutation rate is normalized by the length of the root sequence.
#' @param rate_model A list of GTR+I+G model parameters. Expects a list
#'   of 13 parameters:  (6 rate parameters) `a2c`, `a2g`, `a2t`, `c2g`,
#'  `c2t`, `g2t`, (nucleotide frequencies:) `fa, fc, fg, ft`,
#'  (proportion of sites with no variation:) `i`, (site specific
#'  heterogeneity shape parameter) `alpha`, and number of categories for
#'  discretized gamma heterogeneity (`ncat`).
#' @seealso geneology_to_phylogeny_bpb
#' @importFrom rngtools RNGseed
#' @export
generate_sequences <- function(phylogeny, branch_rate, root_sequence,
                               rng_seed = -1, rate_model,
                               rate_per_nt = FALSE) {
    # We need to convert the phylogeny into a newick tree for seq-gen
    # use mutations per site as branch lengths
    # RNG seed is not supported in phyclust. Instead, the R rng is used.
    if (rng_seed != -1) {
        # Get R's current RNG seed
        rng_seed_store <- rngtools::RNGseed()  # Get the current RNG seed
        on.exit(rngtools::RNGseed(rng_seed_store))  # Restore the RNG seed
        set.seed(rng_seed)
        # Need to fix this. Make seed explicit and explicitly manage state
    }
    # Ensure the root sequence provided is a single character string, not a vector
    if (length(root_sequence) != 1) {
        # Convert to a single character string
        root_sequence <- paste0(root_sequence, collapse = "")
    }
    phylogeny_local <- phylogeny  # Don't modify the input argument
    if (rate_per_nt) {
        # If the given mutation rate is in mutations per sequence, we need to
        # re-normalize it to mutations per site.
        phylogeny_local[, 4] <- phylogeny_local[, 4] / length(root_sequence)
    }
    # Convert branch lengths in time to in E[subs per site]
    phylogeny_local[, 4] <- phylogeny_local[, 4] * branch_rate

    # Build seq-gen call components
    newick_string <- SEEPS::phylogeny_to_newick(phylogeny = phylogeny_local,
                                                mode = "mean",
                                                label_mode = "abs")

    newick_string <- SEEPS::add_root_to_newick(newick_string)  # Add in the root node
    model_string <- " -mGTR"
    rate_string <- paste0("-r", rate_model["a2c"], ",", rate_model["a2g"], ",",
                          rate_model["a2t"], ",", rate_model["c2g"], ",",
                          rate_model["c2t"], ",", rate_model["g2t"])
    freq_string <- paste0(" -f", rate_model["fa"], ",", rate_model["fc"], ",",
                          rate_model["fg"], ",", rate_model["ft"])
    i_string <- paste0(" -i", rate_model["i"])
    gamma_string <- paste0(" -a", rate_model["alpha"], " ", " -g", rate_model["ncat"])
    # relaxed phylip format is easier to parse,
    fmt_string <- paste(" -op", "-k1", "-q")  # -q for quiet

    input <- paste0("1 ", nchar(root_sequence), "\nroot ", root_sequence,
                    "\n1\n", newick_string)
    # Build the option string
    call_opts <- paste(model_string, rate_string, freq_string,
                      i_string, gamma_string, fmt_string)

    # Call seq-gen
    # This uses the phyclust C library under-the-hood to avoid an issue with
    # the input parser in phyclust's R interface.
    seqgen_result <- seqgen(input = input, opts = call_opts)

    # Convert to a fasta file
    fasta <- phylip_to_fasta(seqgen_result)
    return(fasta)
}

phylip_to_fasta <- function(phylip_string) {
    # convert a reduced phylip file to a fasta file
    lines <- unlist(strsplit(phylip_string, split = "\n"))
    fasta_vector <- character(length(lines))
    index <- 1
    for (line in lines[2:length(lines)]) {
        if (length(line) > 0) {
            tokens <- unlist(strsplit(line, split = " "))
            seq_name <- tokens[1]
            seq <- tokens[length(tokens)]  # get the last token
            fasta_string <- paste0(">", seq_name, "\n", seq, "\n")
            fasta_vector[index] <- fasta_string
            index <- index + 1
        }
    }
    return(fasta_vector)
}

seqgen <- function(input, opts) {
    # If phyclust is not loaded, we need to load it here

    # Setup a temp file for the input data
    tmp_fname_data <- paste0("seqgen.", Sys.getpid(), ".temp.data")
    tmp_fname_work <- paste0("seqgen.", Sys.getpid(), ".temp.work")
    # Tie name to PID so that it is unique over all instances of R
    tmp_file_data <- tempfile(tmp_fname_data)
    tmp_file_work <- tempfile(tmp_fname_work)
    write(input, file = tmp_file_data, sep = "\n")
    sopts <- unlist(strsplit(opts, " "))
    sopts <- sopts[sopts != ""]  # Drop empty strings
    cmd <- c("seq-gen", sopts, tmp_file_data)
    # Call seq-gen
    .Call("R_seq_gen_main", cmd, tmp_file_work, PACKAGE = "phyclust")
    data <- readLines(con = tmp_file_work, warn = FALSE)
    # Drop any empty lines from the output
    data <- data[data != ""]
    # We're done. Release the temp files.
    unlink(tmp_file_data)
    unlink(tmp_file_work)
    return(data)
}

#' Convert a fasta string output by seq-gen to a dataframe
#'
#' A utility function to build a dataframe from a fasta string.
#'
#'  @param fasta_string A string containing the fasta output from seq-gen
#'
#' @return A dataframe with the "seq" and "name" columns
#' @export
fasta_string_to_dataframe <- function(fasta_string, trim = TRUE, include_root = FALSE) {
    # Convert a fasta string to a dataframe
    lines <- unlist(strsplit(fasta_string, split = "\n"))
    data <- data.frame(seq = I(list()), name = I(list()))
    index <- 1
    skip_flag <- FALSE
    for (line in lines) {
        if (length(line) > 0) {
            if (substr(line, 1, 1) == ">") {
                # This is a header line
                # if trim , check that the name is not zero and has characters
                # if not, skip this line
                if (trim) {
                    # Check name before we add it
                    if (substring(line, 2, nchar(line)) == "") {
                        skip_flag <- TRUE
                        next
                    }
                    if (substring(line, 2, nchar(line)) == "0_") {
                        skip_flag <- TRUE
                        next
                    }
                    if (!include_root && substring(line, 2, nchar(line)) == "root") {
                        skip_flag <- TRUE
                        next
                    }
                }
                data[[index, "name"]] <- substring(line, 2, nchar(line) - 1)
                # Strip off the leading ">" and tailing "_"
            } else if (skip_flag) {
                # We've already skipped this line, so skip the sequence
                skip_flag <- FALSE
                # Reset the flag
                next
                # Go to next iteration
            } else {
                # This is a sequence line
                data[[index, "seq"]] <- line
                index <- index + 1
            }
        }
    }
    return(data)
}
