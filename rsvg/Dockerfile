FROM ubuntu:18.10

RUN mkdir /build /files

# Install dependecies listed at https://gitlab.gnome.org/GNOME/librsvg/blob/master/COMPILING.md.
RUN apt-get update && apt-get install -y \
  gcc make rustc cargo \
  automake autoconf libtool gettext itstool \
  libgdk-pixbuf2.0-dev libgirepository1.0-dev \
  gtk-doc-tools git \
  libxml2-dev libcairo2-dev libpango1.0-dev \
  && rm -rf /var/lib/apt/lists/*

# Install required cairo verion.
WORKDIR /build
RUN curl -O https://www.cairographics.org/releases/cairo-1.16.0.tar.xz && tar xf cairo-1.16.0.tar.xz
WORKDIR /build/cairo-1.16.0
RUN ./configure && make && make install

# Install required croco version.
WORKDIR /build
RUN curl -L https://github.com/GNOME/libcroco/archive/0.6.13.tar.gz -o libcroco-0.6.13.tar.gz \
  && tar -xvzf libcroco-0.6.13.tar.gz
WORKDIR /build/libcroco-0.6.13
RUN ./autogen.sh && make && make install

# Install librsvg.
WORKDIR /build
RUN curl -L https://github.com/GNOME/librsvg/archive/2.46.3.tar.gz -o librsvg-2.46.3.tar.gz \
  && tar -xvzf librsvg-2.46.3.tar.gz
WORKDIR /build/librsvg-2.46.3
ENV PATH="$PATH:/usr/lib/x86_64-linux-gnu/gdk-pixbuf-2.0"
RUN ./autogen.sh && make && make install

# Configure LD path.
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"
RUN ldconfig

# Clean up.
RUN \rm -rf /build
RUN apt-get clean

# Run rsvg-convert when the container launches
WORKDIR /files
ENTRYPOINT ["rsvg-convert"]
CMD ["--help"]
