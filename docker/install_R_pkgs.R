options(repos=structure(c(CRAN='https://cloud.r-project.org')))
library("groundhog")


pkgs <- c(
    "optparse",
    "stringi",
    "RCurl",
    "rmdformats",
    "knitr",
    "rmarkdown",
    "DT",
    "data.table",
    "R.utils",
    "ggplot2",
    "skimr",
    "dplyr",
    "remotes"
)
groundhog.library(pkgs, "2022-04-01")
groundhog.library("Hmisc", "2022-04-01")

bioc_version = Sys.getenv("BIOC_VERSION")
BiocManager::install(version=bioc_version)
BiocManager::install(c("biovizBase", "ggbio", "Homo.sapiens"))

remotes::install_github("exaexa/scattermore@v0.8")
remotes::install_github("roman-tremmel/ggfastman@v.1.2")
