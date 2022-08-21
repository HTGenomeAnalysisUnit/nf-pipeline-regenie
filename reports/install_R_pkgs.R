#This can help creating an appropriate R environment from scratch
#Alternatively, you can use the renev.lock file and renv::restore()
options(repos=structure(c(CRAN='https://cloud.r-project.org')))

install.packages('remotes','BiocManager')
remotes::install_github('rstudio/renv@0.15.5')

renv::install(c("DT","data.table","R.utils","ggplot2","skimr","kableExtra","dplyr"))
renv::install("exaexa/scattermore")
renv::install("bioc::ggbio")
renv::install("bioc::Homo.sapiens")
renv::install("roman-tremmel/ggfastman")

renv::snapshot()