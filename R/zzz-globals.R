# nocov start
# Silence R CMD check NOTEs for non-standard evaluation patterns
# (foreach loop variables and pipe-style data manipulation).
utils::globalVariables(c("payload", "marker_name", "n"))
# nocov end
