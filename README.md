# Skin sample merging

Some scripts and a dockerfile for performing 10x sample merging with
the Harmony R package.


## Running the environment

The docker environment can be launched with:

```shell
docker run -p 8888:8888 -v $PWD:/data -it --entrypoint --entrypoint /usr/local/bin/jupyter -v $PWD:/data cvaske/cho-skin-cd45:harmony-v1.2.0_0 lab --ip=0.0.0.0 --allow-root
```

and then open the URL that is pasted after Jupyter starts up.

This will give you access to an R kernel for interactive plotting with
Seurat/Harmony.

The scripts can also be run to create more permanent and reproducible
outputs.
