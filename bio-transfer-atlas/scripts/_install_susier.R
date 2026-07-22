dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
install.packages(c("susieR", "Rcpp", "RcppArmadillo"), repos = "https://cloud.r-project.org")
cat("susieR=", requireNamespace("susieR", quietly = TRUE), "\n")
