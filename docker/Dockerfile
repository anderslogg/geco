# Use Phusion base image (minimal Docker-friendly Ubuntu)
FROM phusion/baseimage:focal-1.2.0

# Set some variables
ENV USER geco
ENV HOME /home/$USER
ENV DIR geco

# Install essential packages
RUN apt-get update
RUN apt-get install -y locales sudo build-essential cmake pip

# Install FEniCS
RUN apt-get install -y software-properties-common
RUN add-apt-repository ppa:fenics-packages/fenics
RUN apt-get update
RUN apt-get install -y fenics

# Clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Add user
RUN useradd -m $USER -G sudo && \
    echo "$USER:$USER" | chpasswd && \
    echo "$USER ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers
USER $USER

# Create shared volume
VOLUME $HOME/$DIR
WORKDIR $HOME/$DIR

# Generate welcome message printed at login
COPY Welcome $HOME/.welcome
RUN echo "cat $HOME/.welcome" >> $HOME/.bashrc

# Start login shell
CMD ["/bin/bash", "-i"]
