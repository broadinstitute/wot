FROM ubuntu:trusty
RUN apt-get update && apt-get install --no-install-recommends -y \
		curl \
		ca-certificates \
		gcc \
		unzip

RUN curl -k -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash ./miniconda.sh -b -p /root/miniconda
RUN rm -f ./miniconda.sh
ENV PATH=/root/miniconda/bin:$PATH
ENV MPLBACKEND="PS"
RUN conda install -y cython numpy pandas scikit-learn scipy h5py
RUN conda install -y -c conda-forge pot
RUN curl -L -k -o wot.zip https://github.com/broadinstitute/wot/archive/develop.zip && unzip wot.zip
RUN cd /wot-develop && pip install -e .
RUN rm -f wot.zip


