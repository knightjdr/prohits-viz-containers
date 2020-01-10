FROM ubuntu:18.04 as builder

# Install prerequisites.
RUN apt-get update && apt-get install -y \
  build-essential \
  curl

# Create directories.
RUN mkdir /app /build /files

# Download nestedcluster prerequisites.
WORKDIR /build
RUN curl -O http://gnu.mirror.iweb.com/gsl/gsl-2.6.tar.gz \
  && tar -xvzf gsl-2.6.tar.gz \
  && cd /build/gsl-2.6 \
  && ./configure \
  && make \
  && make install
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

# Download and build nestedcluster.
RUN curl -o biclust.tar.gz -L https://sourceforge.net/projects/nestedcluster/files/biclust_0.0.1.tar.gz/download \
  && tar -xvzf biclust.tar.gz \
  && cd /build/biclust \
  && /build/biclust/compile

# Create runtime layer.
FROM ubuntu:18.04
COPY --from=builder /build/biclust/nestedcluster /app/nestedcluster
COPY --from=builder /lib64/ld-linux-x86-64.so.2 /lib64/ld-linux-x86-64.so.2
COPY --from=builder \
  /usr/local/lib/libgsl.so.25.0.0 \
  /usr/local/lib/libgslcblas.so.0.0.0 \
  /usr/local/lib/
COPY --from=builder \
  /lib/x86_64-linux-gnu/libm-2.27.so \
  /lib/x86_64-linux-gnu/libc-2.27.so \
  /lib/x86_64-linux-gnu/
RUN ln -sf /usr/local/lib/libgsl.so.25.0.0 /usr/local/lib/libgsl.so.25 \
  && ln -sf /usr/local/lib/libgslcblas.so.0.0.0 /usr/local/lib/libgslcblas.so.0 \
  && ln -sf /lib/x86_64-linux-gnu/libm-2.27.so /lib/x86_64-linux-gnu/libm.so.6 \
  && ln -sf /lib/x86_64-linux-gnu/libc-2.27.so /lib/x86_64-linux-gnu/libc.so.6 
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"
RUN ldconfig

# Install and configure R.
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
  && apt-get install -y r-base \
  && Rscript -e 'install.packages("gplots")'

COPY ./nestedcluster/script.bash /app/
COPY ./nestedcluster/biclustering.R /app/
RUN chmod +x /app/script.bash

WORKDIR /files
ENTRYPOINT ["/app/script.bash"]
