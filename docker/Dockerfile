FROM continuumio/anaconda3:latest
SHELL ["/bin/bash", "-c"]

RUN apt-get update \
&& conda update -y -n base -c defaults conda \
&& conda install -y -c conda-forge pot \
&& pip install --upgrade pip \
&& pip install wot==$WOT_VERSION
