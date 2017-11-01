FROM centos

# basic packages
RUN yum update yum
RUN yum install -y \
    gcc \
    make \
    git \
    wget \
    zlib-devel \
    bzip2 \
    bzip2-devel \
    xz-devel \
    ncurses-devel

# bwa: requiring [git, gcc, make, zlib-devel]
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    make && chmod 755 ./bwa && \
    cp ./bwa /bin/bwa && cd /

# samtools: requiring [wget, gcc, make, bzip2, zlib-devel, bzip2-devel, xz-devel]
RUN wget https://github.com/samtools/samtools/releases/download/1.4.1/samtools-1.4.1.tar.bz2 && \
    tar -xvf samtools-1.4.1.tar.bz2 && \
    cd samtools-1.4.1 && \
    make && make install && \
    mv samtools /bin/ && cd /

# pip: requiring [curl, python, git]
RUN curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py"
RUN python get-pip.py

# GenomonFisher: requiring [pysam, scipy]
RUN pip install pysam scipy
RUN pip install git+https://github.com/Genomon-Project/GenomonFisher.git

# Entrypoint
ADD ./main.sh /bin
ENTRYPOINT /bin/main.sh
