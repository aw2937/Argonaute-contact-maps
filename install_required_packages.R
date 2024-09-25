#Do not forget to set the working directory to the source file location.

# Read the list of required packages
required_packages <- readLines("required_packages.txt")

# Install the packages
install.packages(required_packages)

if (!require("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
BiocManager::install("ggtree")
if (!require("devtools", quietly = TRUE)) {install.packages("devtools")}
devtools::install_github("mattflor/chorddiag")