# Code for scRNAseq BRCA Atlas

To install the R packages for the R analysis done in this project run the following in a R session:
```
packages <- readLines("scripts/requirements.txt")
install.packages(packages)
```

To setup the conda environments for the python analysis done in this project run the commands found in:
```
scripts/conda_install.sh
scripts/r-sceasy.sh
```
