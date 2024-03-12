#!/bin/bash


#build docker image
# docker build \
# --platform linux/amd64 \
# -t mckeowr1/prep_sims .

docker build -t mckeowr1/prep_sims_m2 --platform linux/arm64/v8 .

#check if doker image was built
docker run -ti mckeowr1/prep_sims_m2 sh

#Tag image with a version 
docker image tag mckeowr1/asess_sims:latest mckeowr1/prep_sims:1.0

#push the image to docker hub
docker push mckeowr1/prep_sims:1.0