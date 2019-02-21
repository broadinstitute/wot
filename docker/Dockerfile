FROM continuumio/anaconda3:latest
SHELL ["/bin/bash", "-c"]
ADD https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/monitor_script.sh /software/
RUN apt-get update \
&& conda update -y -n base -c defaults conda \
&& conda install -y -c conda-forge pot \
&& pip install --upgrade pip \
&& pip install wot scanpy \
&& chmod a+rx /software/*
ENV PATH=/software:$PATH
