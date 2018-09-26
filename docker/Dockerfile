FROM continuumio/miniconda3:latest
RUN conda install -q -y \
cython \
h5py \
numpy \
pandas \
scanpy \
scikit-learn \
scipy \
&& conda install -q -y -c conda-forge pot \
&& conda clean --all -y
RUN pip install wot

