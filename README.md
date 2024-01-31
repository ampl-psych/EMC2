# EMC2
This repository is work in progress.

To install the R package, and its dependencies you can, download the tar.gz:
https://github.com/ampl-psych/EMC2/archive/refs/heads/main.tar.gz

And install it (and its dependencies) using:
install.packages("emc2-main.tar.gz", repos = NULL, type = "source",dependencies=TRUE)

Alternatively, link the Rstudio project to your github account. (Google should
help this differs per version and OS), install the remotes package, and run:

remotes::install_github("ampl-psych/EMC2",dependencies=TRUE)
