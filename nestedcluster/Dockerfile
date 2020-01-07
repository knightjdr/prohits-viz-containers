FROM alpine:3.11.2 as builder

# Add build prerequisites.
RUN apk add build-base

# Create directories.
RUN mkdir /build
RUN mkdir /files

# Download prerequisites.
WORKDIR /build
RUN wget http://gnu.mirror.iweb.com/gsl/gsl-2.6.tar.gz
RUN tar -xvzf gsl-2.6.tar.gz
WORKDIR /build/gsl-2.6
RUN ./configure
RUN make
RUN make install
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

# Download and build nestedcluster.
WORKDIR /build
RUN wget -O biclust.tar.gz https://sourceforge.net/projects/nestedcluster/files/biclust_0.0.1.tar.gz/download
RUN tar -xvzf biclust.tar.gz
WORKDIR /build/biclust
RUN /build/biclust/compile

# Create runtime layer.
FROM alpine:3.11.2
COPY --from=builder /build/biclust/nestedcluster /app/
COPY --from=builder /usr/local/lib/libgsl* usr/local/lib/
COPY --from=builder /files /files
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"

WORKDIR /files
ENTRYPOINT ["/app/nestedcluster"]
