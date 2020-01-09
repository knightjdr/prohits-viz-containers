FROM ubuntu:18.10 as builder

# Install dependecies listed at https://gitlab.gnome.org/GNOME/librsvg/blob/master/COMPILING.md.
RUN apt-get update && apt-get install -y \
  gcc make rustc cargo \
  automake autoconf libtool gettext itstool \
  libgdk-pixbuf2.0-dev libgirepository1.0-dev \
  gtk-doc-tools git \
  libxml2-dev libcairo2-dev libpango1.0-dev

RUN mkdir /build
WORKDIR /build

# Install required cairo verion.
RUN curl -O https://www.cairographics.org/releases/cairo-1.16.0.tar.xz \
  && tar xf cairo-1.16.0.tar.xz \
  && cd /build/cairo-1.16.0 \
  && ./configure \
  && make \
  && make install

# Install required croco version.
RUN curl -L https://github.com/GNOME/libcroco/archive/0.6.13.tar.gz -o libcroco-0.6.13.tar.gz \
  && tar -xvzf libcroco-0.6.13.tar.gz \
  && cd /build/libcroco-0.6.13 \
  && ./autogen.sh \
  && make \
  && make install

# Install librsvg.
ENV PATH="$PATH:/usr/lib/x86_64-linux-gnu/gdk-pixbuf-2.0"
RUN curl -L https://github.com/GNOME/librsvg/archive/2.46.3.tar.gz -o librsvg-2.46.3.tar.gz \
  && tar -xvzf librsvg-2.46.3.tar.gz \
  && cd /build/librsvg-2.46.3 \
  && ./autogen.sh \
  && make \
  && make install

# Create runtime layer.
FROM ubuntu:18.10
COPY --from=builder /usr/local/bin/rsvg-convert /usr/local/bin/rsvg-convert
COPY --from=builder /lib64/ld-linux-x86-64.so.2 /lib64/ld-linux-x86-64.so.2
COPY --from=builder \
  /usr/local/lib/librsvg-2.so.2 \
  /usr/local/lib/libcairo.so.2 \
  /usr/local/lib/libcroco-0.6.so.3 \
  /usr/local/lib/
COPY --from=builder \
  /usr/lib/x86_64-linux-gnu/libgio-2.0.so.0 \
  /usr/lib/x86_64-linux-gnu/libgobject-2.0.so.0 \
  /usr/lib/x86_64-linux-gnu/libglib-2.0.so.0 \
  /usr/lib/x86_64-linux-gnu/libgdk_pixbuf-2.0.so.0 \
  /usr/lib/x86_64-linux-gnu/libpangocairo-1.0.so.0 \
  /usr/lib/x86_64-linux-gnu/libpango-1.0.so.0 \
  /usr/lib/x86_64-linux-gnu/libxml2.so.2 \
  /usr/lib/x86_64-linux-gnu/libgmodule-2.0.so.0 \
  /usr/lib/x86_64-linux-gnu/libffi.so.6.0.4 \
  /usr/lib/x86_64-linux-gnu/libpixman-1.so.0 \
  /usr/lib/x86_64-linux-gnu/libfontconfig.so.1 \
  /usr/lib/x86_64-linux-gnu/libfreetype.so.6 \
  /usr/lib/x86_64-linux-gnu/libpng16.so.16 \
  /usr/lib/x86_64-linux-gnu/libxcb-shm.so.0 \
  /usr/lib/x86_64-linux-gnu/libxcb.so.1 \
  /usr/lib/x86_64-linux-gnu/libxcb-render.so.0 \
  /usr/lib/x86_64-linux-gnu/libXrender.so.1 \
  /usr/lib/x86_64-linux-gnu/libX11.so.6 \
  /usr/lib/x86_64-linux-gnu/libXext.so.6 \
  /usr/lib/x86_64-linux-gnu/libpangoft2-1.0.so.0 \
  /usr/lib/x86_64-linux-gnu/libthai.so.0 \
  /usr/lib/x86_64-linux-gnu/libfribidi.so.0 \
  /usr/lib/x86_64-linux-gnu/libicuuc.so.60 \
  /usr/lib/x86_64-linux-gnu/libXau.so.6 \
  /usr/lib/x86_64-linux-gnu/libXdmcp.so.6 \
  /usr/lib/x86_64-linux-gnu/libharfbuzz.so.0 \
  /usr/lib/x86_64-linux-gnu/libdatrie.so.1 \
  /usr/lib/x86_64-linux-gnu/libicudata.so.60 \
  /usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.25 \
  /usr/lib/x86_64-linux-gnu/libgraphite2.so.3 \
  /usr/lib/x86_64-linux-gnu/
COPY --from=builder \
  /lib/x86_64-linux-gnu/libm.so.6 \
  /lib/x86_64-linux-gnu/libpthread.so.0 \
  /lib/x86_64-linux-gnu/libc.so.6 \
  /lib/x86_64-linux-gnu/libdl.so.2 \
  /lib/x86_64-linux-gnu/libgcc_s.so.1 \
  /lib/x86_64-linux-gnu/libz.so.1.2.11 \
  /lib/x86_64-linux-gnu/libselinux.so.1 \
  /lib/x86_64-linux-gnu/libresolv.so.2 \
  /lib/x86_64-linux-gnu/libmount.so.1.1.0 \
  /lib/x86_64-linux-gnu/libpcre.so.3.13.3 \
  /lib/x86_64-linux-gnu/librt.so.1 \
  /lib/x86_64-linux-gnu/liblzma.so.5.2.2 \
  /lib/x86_64-linux-gnu/libblkid.so.1.1.0 \
  /lib/x86_64-linux-gnu/libexpat.so.1 \
  /lib/x86_64-linux-gnu/libuuid.so.1.3.0 \
  /lib/x86_64-linux-gnu/libbsd.so.0 \
  /lib/x86_64-linux-gnu/
RUN ln -sf /lib/x86_64-linux-gnu/libpcre.so.3.13.3 /lib/x86_64-linux-gnu/libpcre.so.3 \
  && ln -sf /lib/x86_64-linux-gnu/libz.so.1.2.11 /lib/x86_64-linux-gnu/libz.so.1 \
  && ln -sf /lib/x86_64-linux-gnu/libblkid.so.1.1.0 /lib/x86_64-linux-gnu/libblkid.so.1 \
  && ln -sf /lib/x86_64-linux-gnu/liblzma.so.5.2.2 /lib/x86_64-linux-gnu/liblzma.so.5 \
  && ln -sf /lib/x86_64-linux-gnu/libuuid.so.1.3.0 /lib/x86_64-linux-gnu/libuuid.so.1 \
  && ln -sf /lib/x86_64-linux-gnu/libmount.so.1.1.0 /lib/x86_64-linux-gnu/libmount.so.1 \
  && ln -sf /usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.25 /usr/lib/x86_64-linux-gnu/libstdc++.so.6 \
  && ln -sf /usr/lib/x86_64-linux-gnu/libffi.so.6.0.4 /usr/lib/x86_64-linux-gnu/libffi.so.6

COPY --from=builder \
  /etc/fonts/fonts.conf \
  /etc/fonts/conf.avail \
  /etc/fonts/conf.d \
  /etc/fonts/
COPY --from=builder /usr/share/fonts /usr/share/fonts
COPY ./rsvg/fonts/* /usr/share/fonts/truetype/arial/

# Configure LD path.
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"
RUN ldconfig

RUN mkdir /files

# Run rsvg-convert when the container launches.
WORKDIR /files
ENTRYPOINT ["rsvg-convert"]
CMD ["--help"]
