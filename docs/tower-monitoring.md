# Monitoring with Tower

To easily monitor pipeline execution, we suggest to use Nextflow Tower. 

1. First, register to the [Nextflow tower](https://cloud.tower.nf/) using your GitHub or Google account. 
2. Then, click on your profile (upper right corner) and select `Your tokens`. Then click add token and follow the instructions to create a new token. Make sure to copy the generated token, since you were not able to see it again.
3. Add the following to `single_project.conf` file or create a separate `tower.conf` file and then run the pipeline adding `-c tower.conf`:

```nextflow
tower {
  enabled = true
  accessToken = 'your_token_here'
}
```

Now when the pipeline is running you should be able to monitor progress in the [Nextflow tower](https://cloud.tower.nf/) under your runs.