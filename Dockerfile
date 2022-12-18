FROM rocker/tidyverse:4.2.2
RUN apt-get update && apt install -y libglpk40 libxt6 ghostscript
WORKDIR /home/rstudio
RUN R -e "install.packages('renv', version='0.16.0')"
# COPY renv.lock renv.lock
# RUN R -e "renv::restore()"

