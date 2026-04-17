if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos="https://cloud.r-project.org")
if (!requireNamespace("coRdon", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
    BiocManager::install("coRdon", update=FALSE, ask=FALSE)
}
devtools::install_github("jlw-ecoevo/gRodon2", upgrade="never")
library(gRodon)
cat("gRodon installed OK\n")
