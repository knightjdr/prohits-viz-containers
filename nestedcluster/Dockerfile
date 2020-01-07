FROM ubuntu:18.04

# Install prerequisites.
RUN apt-get update && apt-get install -y \
  build-essential \
  curl

# Create directories.
RUN mkdir /app /build /files

# Download nestedcluster prerequisites.
WORKDIR /build
RUN curl -O http://gnu.mirror.iweb.com/gsl/gsl-2.6.tar.gz && \
  tar -xvzf gsl-2.6.tar.gz
WORKDIR /build/gsl-2.6
RUN ./configure && \
  make && \
  make install
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

# Download and build nestedcluster.
WORKDIR /build
RUN curl -o biclust.tar.gz -L https://sourceforge.net/projects/nestedcluster/files/biclust_0.0.1.tar.gz/download && \
  tar -xvzf biclust.tar.gz
WORKDIR /build/biclust
RUN /build/biclust/compile
RUN mv /build/biclust/nestedcluster /app/

# Install and configure R.
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y r-base
RUN Rscript -e 'install.packages("gplots")'

COPY ./nestedcluster/script.bash /app/
RUN chmod +x /app/script.bash
COPY ./nestedcluster/biclustering.R /app/

WORKDIR /files
ENTRYPOINT ["/app/script.bash"]
