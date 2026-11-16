user_lib <- Sys.getenv("R_LIBS_USER", unset = file.path(Sys.getenv("HOME"), ".local", "R", "library"))
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))
if (!requireNamespace("qtl2", quietly = TRUE)) {
  install.packages("qtl2", repos = "https://cloud.r-project.org", lib = user_lib)
}
cat(as.character(utils::packageVersion("qtl2")), "\n")
