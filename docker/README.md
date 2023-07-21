# Pipeline containers

All the analysis steps of the pipeline depend on a main container defined in the `Dockerfile`.

An additional Dockerfile is provided to generate reports since this is more likely to change frequently as we improve graphics and plots. This is needed only when `make_report` is true. The container for the reports is defined in `Dockerfile_reports` with a supporting python scripts to install gwaslab reference data.