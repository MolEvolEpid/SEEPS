# Utility functions to load in data to the package.
# This directory contains any raw data. To add an organism to the
# reference set, provide three files with names:
# `seq_{NAME}.{csv, fasta, metadata}`.
# The fasta file should contain exactly one reference sequence.
# The csv file should contain a list of annotated regions for the sequence
# These will be exposed to the user. Each annotated region contains a name
# Each annotated region to the sequence and a start and end position.
# The retrieval will obtain both the start and end position.
# Avoid long names, quotes are strictly prohibited.
#
# The metadata file contains any other important/legal information such as
# citations or reference terms.
# Acession IDs should be provided within the fasta file.
#' @importFrom tools file_path_sans_ext
#' @importFrom ape read.dna
#' @importFrom utils read.csv
load_raw_data <- function() {
    # Bind the data and force the loading of reference sequence data
    pkgenv$seq_store <- get0("seq_store", envir = asNamespace("SEEPS"))
}