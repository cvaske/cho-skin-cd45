# Skin sample merging

Some scripts and a dockerfile for performing 10x sample merging with
the Harmony R package.


## Files in this repository

* `docker` - Dockerfiles for various run environments
* `docker-seurat5.0.0` - runs of Harmony on Seurat 5.0.0
  * `*.R` - scripts for generating the output directories stored on Box.com
  * `*.tsv` - cutoffs and outputs used
  * `*.ipynb` - notebooks used for interactive exploration of the data and code prototyping
  
Both the input HDF5 files and output plots/tables are kept out of git
(due to size) via `.gitignore`. Check Box for those.

TODO: make a full Snakemake or Nextflow workflow so this is more reproducible

## Running the environment

The docker environment can be launched with:

```shell
docker run -p 8888:8888 --oom-kill-disable --memory 12g -v $PWD:/data -it --entrypoint /usr/local/bin/jupyter -v $PWD:/data cvaske/cho-skin-cd45:harmony-v1.2.0_0 lab --ip=0.0.0.0 --allow-root
```

and then open the URL that is pasted after Jupyter starts up.

This will give you access to an R kernel for interactive plotting with
Seurat/Harmony.

The scripts can also be run to create more permanent and reproducible
outputs.
