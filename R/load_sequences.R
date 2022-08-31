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
    # Read all files with prefix seq_ in from current directory
    files <- list.files("./data", pattern = "seq_.*\\.(csv|fasta|metadata)")
    # Get the unique organism names
    organism_names <- unique(strsplit(tools::file_path_sans_ext(files), "seq_"))
    # Loop over organism names and load in the data
    for (organism_name in organism_names) {
        # Get the organism name
        organism_name <- organism_name[2]
        # Get the organism metadata
        metadata <- paste(readLines(paste0("data/seq_", organism_name, ".metadata")),
                          collapse = "\n")
        # Get the organism fasta file
        fasta <- ape::read.dna(paste0("data/seq_", organism_name, ".fasta"),
                                 as.character = TRUE, format = "fasta")
        sequence <- as.character((unlist(as.list(fasta))))

        # Get the organism csv file
        csv <- utils::read.csv(paste0("data/seq_", organism_name, ".csv"),
                          header = TRUE, sep = ",")
        # Get the organism accession ID
        accession_id <- rownames(fasta)
        # Get the organism sequence
        sequence <- sequence
        # Get the organism regions
        # Add the organism to the package
        add_organism(organism_name, accession_id, sequence, csv, metadata)
    }
    return(TRUE)
}

#' Add an organism to the package
add_organism <- function(organism_name, accession_id, sequence, region_map, metadata) {
    # Add the organism to the package's internal environment
    pkgenv$seq_store[[organism_name]] <<- list(  # nolint: object_usage_linter
        accession_id = accession_id,
        sequence = toupper(sequence),  # ape likes to load in lowercase
        region = region_map,
        metadata = metadata
    )
    return(NULL)
}
