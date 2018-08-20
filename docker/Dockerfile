FROM ubuntu:trusty
RUN apt-get update && apt-get install --no-install-recommends -y \
		curl \
		ca-certificates \
		gcc \
		unzip

RUN curl -k -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash ./miniconda.sh -b -p /miniconda
RUN rm -f ./miniconda.sh
ENV PATH=/miniconda/bin:$PATH
ENV MPLBACKEND="PS"
RUN conda update -n base conda
RUN conda install -y cython numpy pandas scikit-learn scipy h5py
RUN conda install -y -c conda-forge pot
COPY wot /wot
RUN cd /wot && pip install -e .



