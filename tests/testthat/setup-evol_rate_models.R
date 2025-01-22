# Fixture for a very simple GTR+I+G model
get_simple_GTR_model <- function() { # nolint: object_name_linter
    model <- SEEPS::generate_rate_model(
        a2c = 1, a2g = 1,
        a2t = 1, c2g = 1,
        c2t = 1, g2t = 1,
        fa = 0.25, fc = 0.25,
        fg = 0.25, ft = 0.25,
        i = 0.5,
        alpha = 0.5, ncat = 2
    )
    return(model)
}
