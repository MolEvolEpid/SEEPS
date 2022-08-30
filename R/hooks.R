# Hooks for package loading
# Any additional work that should be done on package load or attach go here
pkgenv <- new.env(parent = emptyenv())
pkgenv[["seq_store"]] <- list()

.onLoad <- function(libname, pkgname) {
    # We need to load in internal state (data)
    # Load in the fasta and map data so it can obtain sequences

    # This is a function defined in the R/load_sequences.R file
    load_raw_data()  # nolint: object_usage_linter
}

.onAttach <- function(libname, pkgname) {
    # Pretty print a load message
}
