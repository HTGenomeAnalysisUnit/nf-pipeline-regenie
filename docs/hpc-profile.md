# Define HPC profile

To create a profile for your HPC cluster, copy the template file `templates/profile_template.config` and adapt the parameters according to your HPC cluster configuration. 

You can use this to configure the scheduler, the queue names, the number of concurrent jobs, the memory and time limits, etc. for your specific environment. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more details.

Then use the profile in the pipeline execution by providing the new config file and the profile name, something like `-c myprofile.config -profile singularity,<profile_name>`.