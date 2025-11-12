To download the singularity image

https://quay.io/repository/biocontainers/xengsort?tab=tags

```
module load apptainer/1.0.2
apptainer build xengsort209.sif docker://quay.io/biocontainers/xengsort:2.0.9--pyhdfd78af_0  # download protocol
```

Run the image:
```
apptainer run xengsort209.sif xengsort index --help
apptainer run xengsort209.sif xengsort classify --help
```
