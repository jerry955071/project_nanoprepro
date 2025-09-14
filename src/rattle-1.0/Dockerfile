FROM gcc:9

# install cmake
RUN apt-get update && apt-get install -y cmake git libz-dev libbz2-dev liblzma-dev libncurses5-dev

# git clone rattle
RUN git clone --recurse-submodules https://github.com/comprna/RATTLE --branch v1.0
WORKDIR /RATTLE
RUN ./build.sh

# configure PATH
WORKDIR /
ENTRYPOINT [ "/RATTLE/rattle" ]