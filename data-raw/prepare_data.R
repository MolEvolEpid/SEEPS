load_local_data <- function() {
    # Read all files with prefix seq_ in from current directory
    files <- list.files("data-raw", pattern = "seq_.*\\.(csv|fasta|metadata)",
                        recursive = TRUE)
    # Get the unique organism names
    organism_names <- unique(strsplit(tools::file_path_sans_ext(files), "seq_"))
    # Loop over organism names and load in the data
    seq_store <- list()
    for (organism_name in organism_names) {
        # Get the organism name
        organism_name <- organism_name[2]
        # Get the organism metadata
        metadata <- paste(readLines(paste0("data-raw/seq_", organism_name, ".metadata")),
                          collapse = "\n")
        # Get the organism fasta file
        fasta <- ape::read.dna(paste0("data-raw/seq_", organism_name, ".fasta"),
                                 as.character = TRUE, format = "fasta")
        sequence <- as.character((unlist(as.list(fasta))))

        # Get the organism csv file
        csv <- utils::read.csv(paste0("data-raw/seq_", organism_name, ".csv"),
                          header = TRUE, sep = ",")
        # Get the organism accession ID
        accession_id <- rownames(fasta)
        # Get the organism sequence
        sequence <- sequence
        # Get the organism regions
        # Add the organism to the package
        seq_store <- add_organism(seq_store, organism_name, accession_id,
                                  sequence, csv, metadata)
    }
    # Save to internal data
    usethis::use_data(seq_store, overwrite = TRUE, internal = TRUE)
}

#' Add an organism to the package
add_organism <- function(seq_store, organism_name, accession_id,
                         sequence, region_map, metadata) {
    # Add the organism to the package's internal environment
    seq_store[[organism_name]] <- list(  # nolint: object_usage_linter
        accession_id = accession_id,
        sequence = toupper(sequence),  # ape likes to load in lowercase
        region = region_map,
        metadata = metadata
    )
    return(seq_store)
}

load_local_data()
