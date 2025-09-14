FROM condaforge/miniforge3:25.3.1-0

# Configure the Bioconda channel
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict
    
# Install cutadapt version 5.1
RUN conda create -n cutadapt cutadapt=5.1

ENTRYPOINT [ "conda", "run", "-n", "cutadapt"]