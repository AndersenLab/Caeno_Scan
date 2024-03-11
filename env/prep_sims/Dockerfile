FROM continuumio/miniconda

COPY prep_sims.yml .
RUN \
   conda env update -n root -f prep_sims.yml \
&& conda clean -a