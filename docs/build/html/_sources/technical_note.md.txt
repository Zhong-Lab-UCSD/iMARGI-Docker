# Technical Notes

## Dockerfile

Here we describe the Dockerfile for building iMARGI-Docker. 

- **Linux base:** Ubuntu (version 18.04)
- **Main tools installed:** sra-tools, seqtk, htslib, samtools, bwa, pbgzip, lz4, python3, pairtools, cooler,
  imargi_wrapper.sh, imargi_clean.sh, imargi_parse.sh, imargi_restrict.py, imargi_rsfrags.sh, imargi_stats.sh,
  imargi_convert.sh, imargi_distfilter.sh, imargi_annotate.sh.
- **Working Directory**: `/imargi`

The source code of Dockerfile is shown below:

``` Docker
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
    pip3 install pairtools cooler

RUN apt-get --purge remove -y build-essential git autoconf automake wget

RUN apt-get -y install locales && \
    sed -i -e 's/# \(en_US\.UTF-8 .*\)/\1/' /etc/locale.gen && locale-gen
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

COPY ./imargi_* /usr/local/bin/
RUN chmod +x /usr/local/bin/imargi_* && mkdir /imargi

WORKDIR /imargi
```

## Run iMARGI-Docker with Non-root User

root (id = 0) is the default user within a container. It will cause some permission problem of some files or directories
created by Docker container. So it's better to run iMARGI-Docker container using `-u (--user)`  option to override the
default root user with your own user id (UID).

You can use command `id` in your linux system to get your own UID. For example, my UID is 1043, so I can run iMARGI-Docker
with `-u 1043`, then all the output files and directories are all belong to my user ID.

``` bash
id
uid=1043(frankyan) gid=1048(frankyan)
docker run -u 1043 -v ~/imargi_example:/imargi imargi imargi_wrapper.sh \
    -r hg38 \
    -N HEK_iMARGI \
    -t 16 \
    -g ./ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
    -1 ./data/SRR8206679_1.fastq.gz \
    -2 ./data/SRR8206679_1.fastq.gz \
    -o ./output
```
