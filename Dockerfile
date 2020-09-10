FROM python:3.8-slim

LABEL MAINTAINER="Flavio Lichtenstein <flalix@gmail.com>"

# RUN adduser -D centd
RUN   adduser --system centd

WORKDIR /home/centd

# https://dev.to/faizanbashir/building-python-data-science-container-usingdocker-3f8p
# Building Python Data Science Container using Docker  Faizan Bashir 20 de jan. de 2019
# SOFTWARE PACKAGES
#   * musl: standard C library
#   * lib6-compat: compatibility libraries for glibc
#   * linux-headers: commonly needed, and an unusual package name from Alpine.
#   * build-base: used so we include the basic development packages (gcc)
#   * bash: so we can access /bin/bash
#   * git: to ease up clones of repos
#   * ca-certificates: for SSL verification during Pip and easy_install
#   * freetype: library used to render text onto bitmaps, and provides support font-related operations
#   * libgfortran: contains a Fortran shared library, needed to run Fortran
#   * libgcc: contains shared code that would be inefficient to duplicate every time as well as auxiliary helper routines and runtime support
#   * libstdc++: The GNU Standard C++ Library. This package contains an additional runtime library for C++ programs built with the GNU compiler
#   * openblas: open source implementation of the BLAS(Basic Linear Algebra Subprograms) API with many hand-crafted optimizations for specific processor types
#   * tcl: scripting language
#   * tk: GUI toolkit for the Tcl scripting language
#   * libssl1.0: SSL shared libraries

ENV PACKAGES="\
    dumb-init \
    musl \
    libc6-compat \
    linux-headers \
    build-base \
    bash \
    git \
    ca-certificates \
    freetype \
    libgfortran \
    libgcc \
    libstdc++ \
    openblas \
    tcl \
    tk "

# apt-cache search openblas
# libblas-test - Basic Linear Algebra Subroutines 3, testing programs
# libopenblas-base - Optimized BLAS (linear algebra) library based on GotoBLAS2
# libopenblas-dev - Optimized BLAS (linear algebra) library based on GotoBLAS2

RUN apt-get -y update && apt-get upgrade
# RUN apt-get add libssl1.0
RUN apt-get install -y python3 python3-dev build-essential
RUN apt-get install -y musl-dev libc6-dev
RUN apt-get install -y libopenblas-dev libopenblas-base libfreetype6-dev pkg-config
RUN apt-get install -y dumb-init musl libc6-dev
RUN apt-get install -y bash git ca-certificates freetype2-demos
RUN apt-get install -y gfortran-7 aptitude libstdc++6 libopenblas-dev tcl
# RUN apt-get install -y linux-headers-$(uname -r)


COPY requirements.txt requirements.txt
RUN python -m venv venv
RUN venv/bin/pip install --upgrade pip
RUN venv/bin/pip install -r requirements.txt
RUN venv/bin/pip install gunicorn pymysql

COPY app app
COPY migrations migrations
COPY centd.py config.py boot.sh ./

# cd /projects/src
# cp -r src/ CENTD/
# ce /projects
# mv colaboracoes/ web/CENTD

COPY src/ /home/src/
COPY colaboracoes/covid/fasta/sarscov2_202007' /colaboracoes/covid/fasta/sarscov2_202007/
COPY colaboracoes/biobanco/screening'  /colaboracoes/biobanco/screening/

RUN chmod a+x boot.sh
ENV FLASK_APP centd.py

RUN chown -R centd:centd ./
USER centd

EXPOSE 5000
ENTRYPOINT ["./boot.sh"]
