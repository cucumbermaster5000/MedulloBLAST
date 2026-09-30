# Run from the project root. No credentials or hosting deployment involved.
dir.create(".R-library",showWarnings=FALSE)
.libPaths(c(normalizePath(".R-library"),.libPaths()))
options(repos=c(CRAN="https://cloud.r-project.org"),timeout=600)
if(!requireNamespace("renv",quietly=TRUE))install.packages("renv",lib=".R-library")
renv::restore(project=".",library=normalizePath(".R-library"),prompt=FALSE)
