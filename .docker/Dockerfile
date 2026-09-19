FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="faaf398b7874c2f138bd010f28df01709a2c7798eb3f461a05afb2107584651a"

# Step 2: Retrieve conda environments

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
#   prefix: /conda-envs/2ab973151943ea4481469dc78c1becb6
#   name: R
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - r-base=4.5
#     - r-tidyverse=2.0.0
#     - r-cowplot=1.2.0
#     - r-rcolorbrewer
#     - r-ggrepel=0.9.6
#     - r-scales=1.4.0
#     - r-data.table=1.17.8
#     - r-gprofiler2=0.2.4
#     - bioconductor-methylkit=1.36.0
#     - bioconductor-genomicranges=1.62.1
#     - bioconductor-chipseeker=1.46.1
#     - bioconductor-txdb.mmusculus.ucsc.mm39.knowngene=3.22.0
#     - bioconductor-txdb.hsapiens.ucsc.hg38.knowngene=3.22.0
#     - bioconductor-org.mm.eg.db=3.22.0
#     - bioconductor-org.hs.eg.db=3.22.0
#     - bioconductor-genomeinfodb=1.46.2
#     - bioconductor-rtracklayer=1.70.1
RUN mkdir -p /conda-envs/2ab973151943ea4481469dc78c1becb6
COPY workflow/envs/R.yaml /conda-envs/2ab973151943ea4481469dc78c1becb6/environment.yaml

# Conda environment:
#   source: workflow/envs/bismark.yaml
#   prefix: /conda-envs/4f8aeead075fd4b7d6fcd7fde18e8d05
#   name: bismark
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - pigz=2.8
#     - bioconda:bismark=3.1.0
RUN mkdir -p /conda-envs/4f8aeead075fd4b7d6fcd7fde18e8d05
COPY workflow/envs/bismark.yaml /conda-envs/4f8aeead075fd4b7d6fcd7fde18e8d05/environment.yaml

# Conda environment:
#   source: workflow/envs/deeptools.yaml
#   prefix: /conda-envs/2bdd78eb3d93aba4429dda48c67f8f3c
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
#     - bioconda::ucsc-bigwigaverageoverbed=469
#     - bioconda::bedtools=2.31.1
#     - bioconda::pybedtools=0.10.0
#     - bioconda::snakemake=8.25.5
#     - conda-forge::pandas=2.2.3
#     - bioconda::pyfaidx=0.9.0.3
#     - bioconda::wiggletools=1.2.11
RUN mkdir -p /conda-envs/2bdd78eb3d93aba4429dda48c67f8f3c
COPY workflow/envs/deeptools.yaml /conda-envs/2bdd78eb3d93aba4429dda48c67f8f3c/environment.yaml

# Conda environment:
#   source: workflow/envs/trim_galore.yml
#   prefix: /conda-envs/f25143b31d35300a2e0f8a6febf76392
#   name: trim_galore
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     # Pinned to this exact build: build h9ee0642_0 was linked against a newer
#     # glibc than older hosts (e.g. Ubuntu 20.04, glibc 2.31) provide, causing
#     # `GLIBC_2.32/2.33/2.34 not found` errors. hf1b6044_1 fixes this (glibc <=2.16).
#     - trim-galore=2.2.0=hf1b6044_1
RUN mkdir -p /conda-envs/f25143b31d35300a2e0f8a6febf76392
COPY workflow/envs/trim_galore.yml /conda-envs/f25143b31d35300a2e0f8a6febf76392/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/fc3f74ad3a4edd6540650576745e8ff5 --file /conda-envs/fc3f74ad3a4edd6540650576745e8ff5/environment.yaml && \
    conda env create --prefix /conda-envs/d644a06128ee144c2c0f210018caa0b9 --file /conda-envs/d644a06128ee144c2c0f210018caa0b9/environment.yaml && \
    conda env create --prefix /conda-envs/f008bd698a6714e6623cc359e3d3e6b1 --file /conda-envs/f008bd698a6714e6623cc359e3d3e6b1/environment.yaml && \
    conda env create --prefix /conda-envs/2ab973151943ea4481469dc78c1becb6 --file /conda-envs/2ab973151943ea4481469dc78c1becb6/environment.yaml && \
    conda env create --prefix /conda-envs/4f8aeead075fd4b7d6fcd7fde18e8d05 --file /conda-envs/4f8aeead075fd4b7d6fcd7fde18e8d05/environment.yaml && \
    conda env create --prefix /conda-envs/2bdd78eb3d93aba4429dda48c67f8f3c --file /conda-envs/2bdd78eb3d93aba4429dda48c67f8f3c/environment.yaml && \
    conda env create --prefix /conda-envs/f25143b31d35300a2e0f8a6febf76392 --file /conda-envs/f25143b31d35300a2e0f8a6febf76392/environment.yaml && \
    conda clean --all -y
