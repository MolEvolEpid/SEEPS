# Test based on the "toy" model in the data

test_that("toy model is correctly retrieved with `lookup_sequence_by_name()`", {
    # Simple cases with 10 bases, repeated
    expect_equal(toupper(lookup_sequence_by_name(organism_name = "toy", region_name = "polyA")),
                 rep("A", 10))
    expect_equal(toupper(lookup_sequence_by_name(organism_name = "toy", region_name = "polyC")),
                 rep("C", 10))
    expect_equal(toupper(lookup_sequence_by_name(organism_name = "toy", region_name = "polyG")),
                 rep("G", 10))
    expect_equal(toupper(lookup_sequence_by_name(organism_name = "toy", region_name = "polyT")),
                 rep("T", 10))

    # case AC
    expect_equal(toupper(lookup_sequence_by_name(organism_name = "toy", region_name = "polyAC")),
                 c(rep("A", 10), rep("C", 10)))
    # case CG
    expect_equal(toupper(lookup_sequence_by_name(organism_name = "toy", region_name = "polyCG")),
                 c(rep("C", 10), rep("G", 10)))
    # case GT
    expect_equal(toupper(lookup_sequence_by_name(organism_name = "toy", region_name = "polyGT")),
                 c(rep("G", 10), rep("T", 10)))
    # case ACG
    expect_equal(toupper(lookup_sequence_by_name(organism_name = "toy", region_name = "polyACG")),
                 c(rep("A", 10), rep("C", 10), rep("G", 10)))
    # case CGT
    expect_equal(toupper(lookup_sequence_by_name(organism_name = "toy", region_name = "polyCGT")),
                 c(rep("C", 10), rep("G", 10), rep("T", 10)))

})

test_that("toy model is correctly retrieved with `lookup_sequence_by_index()", {
    # Simple cases with 10 bases, repeated
    expect_equal(toupper(lookup_sequence_by_index(organism_name = "toy", start = 1, stop = 10)),
                 rep("A", 10))
    expect_equal(toupper(lookup_sequence_by_index(organism_name = "toy", start = 11, stop = 20)),
                 rep("C", 10))
    expect_equal(toupper(lookup_sequence_by_index(organism_name = "toy", start = 21, stop = 30)),
                 rep("G", 10))
    expect_equal(toupper(lookup_sequence_by_index(organism_name = "toy", start = 31, stop = 40)),
                 rep("T", 10))

    # case AC
    expect_equal(toupper(lookup_sequence_by_index(organism_name = "toy", start = 1, stop = 20)),
                 c(rep("A", 10), rep("C", 10)))
    # case CG
    expect_equal(toupper(lookup_sequence_by_index(organism_name = "toy", start = 11, stop = 30)),
                 c(rep("C", 10), rep("G", 10)))
    # case GT
    expect_equal(toupper(lookup_sequence_by_index(organism_name = "toy", start = 21, stop = 40)),
                 c(rep("G", 10), rep("T", 10)))
    # case ACG
    expect_equal(toupper(lookup_sequence_by_index(organism_name = "toy", start = 1, stop = 30)),
                 c(rep("A", 10), rep("C", 10), rep("G", 10)))
    # case CGT
    expect_equal(toupper(lookup_sequence_by_index(organism_name = "toy", start = 11, stop = 40)),
                 c(rep("C", 10), rep("G", 10), rep("T", 10)))
})
