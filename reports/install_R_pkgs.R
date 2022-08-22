options(repos=structure(c(CRAN='https://cloud.r-project.org')))
library("groundhog")

pkgs <- c("BiocManager", "remotes","DT","data.table","R.utils","ggplot2","skimr","kableExtra","dplyr")
groundhog.library(pkgs, "2022-04-01")

bioc_version = Sys.getenv("BIOC_VERSION")
BiocManager::install(version=bioc_version)
BiocManager::install("bioc::ggbio")
BiocManager::install("bioc::Homo.sapiens")

remotes::install_github('exaexa/scattermore@v0.8')
remotes::install_github('roman-tremmel/ggfastman@v1.2-9-g367d419')