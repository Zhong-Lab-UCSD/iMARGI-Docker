FROM ubuntu:18.04
ENV TIMEZONE America/Los_Angeles

RUN apt-get update && \
    apt-get install -y \
    git build-essential libz-dev libbz2-dev liblzma-dev libssl-dev libcurl4-gnutls-dev \
    autoconf automake libncurses5-dev wget gawk parallel && \
    cd /tmp && git clone https://github.com/lh3/seqtk.git && \
    cd seqtk && make && make install && \
    cd /tmp && git clone https://github.com/samtools/htslib && \
    cd htslib && autoheader && autoconf && \
    ./configure --prefix=/usr/local && make && make install && \
    cd /tmp && git clone https://github.com/samtools/samtools && \
    cd samtools && autoheader && autoconf && \
    ./configure --prefix=/usr/local && make && make install && \
    cd /tmp && git clone https://github.com/lh3/bwa.git && \
    cd bwa && make && cp bwa /usr/local/bin && \
    cd /tmp && git clone https://github.com/nh13/pbgzip && \
    cd pbgzip && sh autogen.sh && ./configure && make && make install && \
    cd /tmp && git clone https://github.com/lz4/lz4 && \
    cd lz4 && make && make install && \
    cd /tmp && wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.4/sratoolkit.2.9.4-ubuntu64.tar.gz && \
    tar zxvf sratoolkit.2.9.4-ubuntu64.tar.gz && cp -R sratoolkit.2.9.4-ubuntu64/bin/* /usr/local/bin && \
    cd /tmp && rm -rf * 

RUN apt-get install -y python3-dev libopenblas-dev python3-pip && \
    pip3 install numpy cython scipy pandas click && \
    pip3 install pairtools cooler HTSeq

RUN apt-get --purge remove -y build-essential git autoconf automake wget

RUN apt-get -y install locales && \
    sed -i -e 's/# \(en_US\.UTF-8 .*\)/\1/' /etc/locale.gen && locale-gen
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

COPY ./imargi_* /usr/local/bin/
RUN chmod +x /usr/local/bin/imargi_* && mkdir /imargi

WORKDIR /imargi
