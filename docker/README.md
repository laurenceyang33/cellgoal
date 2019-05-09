# Docker container providing cellgoals software (BIG-BOSS and ProtOS)

Accompanies the paper: XXX
to appear in KDD

## Run the image (to build the image see **Building the Image**)
```
docker run --rm -i -v $(pwd):/workdir -t sbrg/cellgoals:1.0 bash
```
From here you can run the test suite and reproduce figures in the paper.
You can also run the algorithms using your own data and genome-scale models.
Please refer to documentation at
https://github.com/laurenceyang33/cellgoal

## Building the Image
1. TODO
1. Build the image
```
$ docker build -t sbrg/cellgoal:1.0 .
```

## Transfer the image to another computer
1. [optional] Save the image
```
docker save -o <image_file> sbrg/cellgoal:1.0
```
2. Copy (rsync, scp) the <image_file> above to your computer. You might want to tar/gzip it first.
3. Load the image
```
docker load -i <image_file>
```
