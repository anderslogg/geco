# Installation

## Downloading the code
Clone this repository using either https:

    git clone https://github.com/anderslogg/geco.git

or ssh:

    git clone git@github.com:anderslogg/geco.git

## Dependencies

GECo depends on FEniCS version 2019.2.0.

The easiest way to get the dependencies is to use the provided Docker image.
The Docker image can be easily built by entering the `docker` directory and then
running the following commands:

    ./docker-build-image
    ./docker-build-container
    ./docker-start-container

The build may take a few minutes. 

## Getting started

To run the demos for GECo, install `geco` using `pip` by running the following
command from inside the package directory:

    pip install .

Alternatively, you can issue the command

    export PYTHONPATH=`pwd`:$PYTHONPATH

This will allow you to run the demos inside the `demos` directory without
actually installing GECo. Note that these commands should be run from inside the
Docker container (unless you have the dependencies installed on your host
system).