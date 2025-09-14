FROM condaforge/miniforge3:25.3.1-0

# Configure channels
RUN conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --add channels nanoporetech && \
    conda config --set channel_priority strict

RUN conda create -n pychopper "nanoporetech::pychopper"
RUN conda run -n pychopper pip install memory-profiler==0.61.0

ENTRYPOINT [ "conda", "run", "-n", "pychopper"]