FROM nfcore/base:1.14
LABEL authors="Olga Botvinnik" \
      description="Docker image containing all software requirements for the nf-core/sicilian pipeline"

# Install the conda environment
COPY environment.yml /
# Install "mamba" for faster dependency resolution and build times
RUN conda install -c conda-forge mamba
RUN mamba env create --quiet -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-sicilian-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-sicilian-1.0dev > nf-core-sicilian-1.0dev.yml
