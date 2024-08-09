FROM --platform=linux/amd64 satijalab/seurat:5.0.0

RUN R -e 'install.packages("harmony")'
RUN python3 -m pip install jupyterlab
RUN R -e "install.packages('IRkernel'); IRkernel::installspec()"

EXPOSE 8888/tcp
