#!/usr/bin/env bash
docker run -e PASSWORD=nopassword -v `pwd`/src:/home/rstudio/R -p 8787:8787 quay.io/comp-bio-aging/gene-modules:latest /init
