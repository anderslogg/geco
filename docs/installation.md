# Installation

Run the command

    pip install .

from inside the package directory.

Alternatively, you can issue the command

    export PYTHONPATH=`pwd`:$PYTHONPATH

This will allow you to run the demos inside the `demos` directory without
actually installing GECo.

## Dependencies

GECo depends on FEniCS version 2019.2.0.

The easiest way to get the dependendencies is to use the provided Docker image.
The Docker image can be easily built by entering the `docker` directory and then
running the following commands:

    ./docker-build-image
    ./docker-build-container
    ./docker-start-container
