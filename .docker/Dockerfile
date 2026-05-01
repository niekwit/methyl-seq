FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="f705961e1d893a1b9d5c99b097fb5eea1a614f3e38a93153dc40b8891ce4d9d5"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v2.10.0/bio/preseq/lc_extrap/environment.yaml
#   prefix: /conda-envs/dbf6d01858bf329ddb05cf98bd61f111
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - preseq =3.2.0
RUN mkdir -p /conda-envs/dbf6d01858bf329ddb05cf98bd61f111
ADD https://github.com/snakemake/snakemake-wrappers/raw/v2.10.0/bio/preseq/lc_extrap/environment.yaml /conda-envs/dbf6d01858bf329ddb05cf98bd61f111/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v7.6.0/bio/fastqc/environment.yaml
#   prefix: /conda-envs/fc3f74ad3a4edd6540650576745e8ff5
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.8.0
RUN mkdir -p /conda-envs/fc3f74ad3a4edd6540650576745e8ff5
ADD https://github.com/snakemake/snakemake-wrappers/raw/v7.6.0/bio/fastqc/environment.yaml /conda-envs/fc3f74ad3a4edd6540650576745e8ff5/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v8.1.1/bio/multiqc/environment.yaml
#   prefix: /conda-envs/d644a06128ee144c2c0f210018caa0b9
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.33
#     - snakemake-wrapper-utils =0.8.0
RUN mkdir -p /conda-envs/d644a06128ee144c2c0f210018caa0b9
ADD https://github.com/snakemake/snakemake-wrappers/raw/v8.1.1/bio/multiqc/environment.yaml /conda-envs/d644a06128ee144c2c0f210018caa0b9/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v9.0.0/bio/samtools/sort/environment.yaml
#   prefix: /conda-envs/5e89b2a5d29ff7b7f49c8dc9669267b4
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.23
#     - snakemake-wrapper-utils =0.8.0
RUN mkdir -p /conda-envs/5e89b2a5d29ff7b7f49c8dc9669267b4
ADD https://github.com/snakemake/snakemake-wrappers/raw/v9.0.0/bio/samtools/sort/environment.yaml /conda-envs/5e89b2a5d29ff7b7f49c8dc9669267b4/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v9.3.0/bio/trim_galore/se/environment.yaml
#   prefix: /conda-envs/9a44d3e05abb1b2d4bbef33d61bf314d
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - trim-galore =0.6.11
RUN mkdir -p /conda-envs/9a44d3e05abb1b2d4bbef33d61bf314d
ADD https://github.com/snakemake/snakemake-wrappers/raw/v9.3.0/bio/trim_galore/se/environment.yaml /conda-envs/9a44d3e05abb1b2d4bbef33d61bf314d/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v9.4.2/bio/seqkit/environment.yaml
#   prefix: /conda-envs/f008bd698a6714e6623cc359e3d3e6b1
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - seqkit =2.13.0
#     - htslib =1.23.1
#     - snakemake-wrapper-utils =0.8.0
RUN mkdir -p /conda-envs/f008bd698a6714e6623cc359e3d3e6b1
ADD https://github.com/snakemake/snakemake-wrappers/raw/v9.4.2/bio/seqkit/environment.yaml /conda-envs/f008bd698a6714e6623cc359e3d3e6b1/environment.yaml

# Conda environment:
#   source: workflow/envs/R.yaml
#   prefix: /conda-envs/b1e5f7e50783d0a664cdf5d4cf9ff01f
#   name: R
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - r-base=4.4
#     - r-tidyverse=2.0.0
#     - r-cowplot=1.2.0
#     - r-rcolorbrewer
#     - r-ggrepel=0.9.6
#     - r-scales=1.4.0
#     - r-data.table=1.17.8
RUN mkdir -p /conda-envs/b1e5f7e50783d0a664cdf5d4cf9ff01f
COPY workflow/envs/R.yaml /conda-envs/b1e5f7e50783d0a664cdf5d4cf9ff01f/environment.yaml

# Conda environment:
#   source: workflow/envs/bismark.yaml
#   prefix: /conda-envs/e254bd4b0d07081bedb140ad007ad7ca
#   name: bismark
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - pigz=2.8
#     - bioconda:bismark=0.25.1
RUN mkdir -p /conda-envs/e254bd4b0d07081bedb140ad007ad7ca
COPY workflow/envs/bismark.yaml /conda-envs/e254bd4b0d07081bedb140ad007ad7ca/environment.yaml

# Conda environment:
#   source: workflow/envs/deeptools.yaml
#   prefix: /conda-envs/d6dfe6c4ba9fc3e1d2289533b5e5696e
#   name: deeptools
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - pigz=2.8
#     - python=3.12
#     - bioconda::deeptools=3.5.6
#     - bioconda::ucsc-bedgraphtobigwig=482
#     - bioconda::bedtools=2.31.1
#     - bioconda::snakemake=8.25.5
#     - bioconda::pyfaidx=0.9.0.3
#     - bioconda::wiggletools=1.2.11
RUN mkdir -p /conda-envs/d6dfe6c4ba9fc3e1d2289533b5e5696e
COPY workflow/envs/deeptools.yaml /conda-envs/d6dfe6c4ba9fc3e1d2289533b5e5696e/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/dbf6d01858bf329ddb05cf98bd61f111 --file /conda-envs/dbf6d01858bf329ddb05cf98bd61f111/environment.yaml && \
    conda env create --prefix /conda-envs/fc3f74ad3a4edd6540650576745e8ff5 --file /conda-envs/fc3f74ad3a4edd6540650576745e8ff5/environment.yaml && \
    conda env create --prefix /conda-envs/d644a06128ee144c2c0f210018caa0b9 --file /conda-envs/d644a06128ee144c2c0f210018caa0b9/environment.yaml && \
    conda env create --prefix /conda-envs/5e89b2a5d29ff7b7f49c8dc9669267b4 --file /conda-envs/5e89b2a5d29ff7b7f49c8dc9669267b4/environment.yaml && \
    conda env create --prefix /conda-envs/9a44d3e05abb1b2d4bbef33d61bf314d --file /conda-envs/9a44d3e05abb1b2d4bbef33d61bf314d/environment.yaml && \
    conda env create --prefix /conda-envs/f008bd698a6714e6623cc359e3d3e6b1 --file /conda-envs/f008bd698a6714e6623cc359e3d3e6b1/environment.yaml && \
    conda env create --prefix /conda-envs/b1e5f7e50783d0a664cdf5d4cf9ff01f --file /conda-envs/b1e5f7e50783d0a664cdf5d4cf9ff01f/environment.yaml && \
    conda env create --prefix /conda-envs/e254bd4b0d07081bedb140ad007ad7ca --file /conda-envs/e254bd4b0d07081bedb140ad007ad7ca/environment.yaml && \
    conda env create --prefix /conda-envs/d6dfe6c4ba9fc3e1d2289533b5e5696e --file /conda-envs/d6dfe6c4ba9fc3e1d2289533b5e5696e/environment.yaml && \
    conda clean --all -y