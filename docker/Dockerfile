FROM ubuntu:trusty
RUN apt-get update && apt-get install --no-install-recommends -y \
		curl \
		ca-certificates \
		git \
		gcc \
		libgl1-mesa-glx

RUN curl -k -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash ./miniconda.sh -b -p /root/miniconda
RUN rm -f ./miniconda.sh
ENV PATH=/root/miniconda/bin:$PATH
RUN conda install -y cython numpy pandas scikit-learn scipy h5py
RUN conda install -y -c conda-forge pot
RUN git clone -b develop https://github.com/broadinstitute/wot.git
WORKDIR /wot
RUN git pull
RUN pip install -e .
WORKDIR /




