# Use Phusion base image (minimal Docker-friendly Ubuntu)
FROM phusion/baseimage:master as base

# Set basic variables
ENV USER geco
ENV HOME /home/$USER
ENV DIR geco

# Add user
RUN useradd -m $USER -G sudo && \
    echo "$USER:$USER" | chpasswd && \
    echo "$USER ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers

# Add shared volume
VOLUME $HOME/$DIR
WORKDIR $HOME/$DIR

# Add welcome message
COPY Welcome $HOME/.welcome
RUN echo "cat $HOME/.welcome" >> $HOME/.bashrc

# Start bash login shell
FROM base as dev
ENTRYPOINT ["/bin/bash", "-l", "-c"]
CMD ["/bin/bash", "-i"]

# Change to root
USER root

# Install essential packages
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    apt-utils gcc g++ openssh-server cmake build-essential gdb gdbserver rsync vim pip

# Install FEniCS
RUN apt-get install software-properties-common
RUN add-apt-repository ppa:fenics-packages/fenics
RUN apt-get update
RUN apt-get install -y fenics

# Change to user
USER $USER
