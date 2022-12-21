FROM rocker/r-ver:4	

RUN apt-get update && apt-get install -y --no-install-recommends python3-pandas python3-numpy python3-matplotlib python3-seaborn python3-sklearn
RUN R -e 'install.packages(c("missForest", "glmnet", "caret", "doParallel", "dbplyr", "randomforest"))'
RUN R -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("biomaRt"); '
