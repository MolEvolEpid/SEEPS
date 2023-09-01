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
    print_load_message()
}

print_load_message <- function() {
    cat(
"                              __   \n",
"                             \\ \\  \n",
" _____________________________\\ \\ \n",
"/_____/_____/_____/_____/_____/ / \n",
"                             /_/  \n",
"   _____ ________________  _____  \n",
"  / ___// ____/ ____/ __ \\/ ___/  \n",
"  \\__ \\/ __/ / __/ / /_/ /\\__ \\   \n",
" ___/ / /___/ /___/ ____/___/ /   \n",
"/____/_____/_____/_/    /____/    \n",
"  __                              \n",
" / /                              \n",
"/ / ______________________________\n",
"\\ \\/_____/_____/_____/_____/_____/\n",
" \\_\\                              \n",
"                                  \n")
}

