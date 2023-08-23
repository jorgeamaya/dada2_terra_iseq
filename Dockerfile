FROM continuumio/miniconda3:23.3.1-0

RUN conda config \
    --add channels defaults \
    --add channels bioconda \
    --add channels conda-forge

#RUN conda install -y python=3.10.0
RUN conda install -y pandas
RUN conda install -y biopython
RUN conda install -y bbmap
RUN conda install -y trim-galore 
RUN conda install -y cutadapt

RUN apt-get update
RUN apt-get install -y curl gnupg1 r-base #trim-galore cutadapt

RUN conda install -y bioconductor-dada2
RUN conda install -y bioconductor-limma
RUN conda install -y muscle

# Specific for google cloud support
RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && apt-get update -y && apt-get install google-cloud-sdk -y

RUN curl https://sdk.cloud.google.com | bash
ENV PATH=/root/google-cloud-sdk/bin/:${PATH}

#Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod

COPY Code Code
