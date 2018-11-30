FROM continuumio/miniconda3:latest
RUN conda install -q -y \
cython \
h5py \
numpy \
pandas \
psutil \
scikit-learn \
scipy \
&& conda install -q -y -c conda-forge pot \
&& conda clean --all -y && \
pip install scanpy wot

