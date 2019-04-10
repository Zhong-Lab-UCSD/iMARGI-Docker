# Technical Notes

- [Technical Notes](#technical-notes)
  - [Dockerfile](#dockerfile)
  - [Install Docker on Different systems](#install-docker-on-different-systems)
  - [Change Docker Memory Settings on Windows and macOS](#change-docker-memory-settings-on-windows-and-macos)
  - [Run iMARGI-Docker with Non-root User](#run-imargi-docker-with-non-root-user)
  - [Solve `bwa index` Failure Problem](#solve-bwa-index-failure-problem)

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

## Install Docker on Different systems

There are official docs of Docker installation guides for different OS, which can be found in
[Docker official webpage](https://docs.docker.com/install/).

Here are only some essential instructions. Install Docker on Linux is the easiest.

- **Linux**: Support the most recent 64 bit stable releases of Ubuntu, Debian, Fedora and CentOS. You need `root` or `sudo`
  privileges. Generally, the following script will automatically install Docker in your system.
  
  ``` Bash
  sudo curl -fsSL https://get.docker.com |sh -
  
  # replace frank with you user name
  sudo usermod -aG docker frank
  ```

- **macOS (modern)**: Docker Desktop for macOS. Support macOS Sierra 10.12 and newer on a Apple computer after 2010.
  
  Download Docker Desktop software for macOS and install.
  [Click here to check instructions](https://docs.docker.com/docker-for-mac/install/)

- **Windows 10 (modern)**: Docker Desktop for Windows. Support the latest Windows 10 (64 bit) Pro, Enterprise or
  Education version.
  
  First, enable virtualization of your CPU (most of modern Intel CPUs support virtualization).
  [Check here to see how to enable it in BIOS.](https://www.isumsoft.com/computer/enable-virtualization-technology-vt-x-in-bios-or-uefi.html)
  Then, turn on Hyper-V. [[Check here to see how to turn on Hyper-V.](https://docs.microsoft.com/en-us/virtualization/hyper-v-on-windows/quick-start/enable-hyper-v)
  Finally, download Docker Desktop software for Windows and install,
  [Click here to check instructions](https://docs.docker.com/docker-for-windows/install/)

- **Legacy solution**: For older Mac and Windows systems that do not meet the requirements of Docker Desktop for Mac and
  Docker Desktop for Windows, you can install Docker Toolbox to use Docker.

  Download Docker Toolbox for macOS or Windows and install.
  [Click here to check instructions](https://docs.docker.com/toolbox/overview/)

## Change Docker Memory Settings on Windows and macOS

Enough memory is important to the iMARGI data processing pipeline. The required amount of memory depends on the size
of reference genome. For human genome, at least 8GB free memory are required by BWA. Hence,** the memory on the machine
needs to be more than 8 GB, which usually is 16 GB**. If the memory is not enough, BWA will generate an empty BAM file,
then it will throw out some strange error information in following steps, such as `"KeyError: 'chr1'"`.

**System memory requirement: 16 GB.**

In addition, if you are using Windows or macOS, there is a memory limit to Docker, which is 2 GB as default. Hence,
you need to increase it to more than 8 GB.

Here are simple instructions of how to change the settings.

- If you are using Docker Desktop for Windows or macOS, you can easily change the settings by right click the
  Docker icon (Whale) in the task bar, then go to Settings -> Advanced to change memory and CPU limits.
  More detail can be found in the Docker official docs of
  [Get started with Docker for Windows](https://docs.docker.com/docker-for-windows/), and
  [Get started with Docker Desktop for Mac](https://docs.docker.com/docker-for-mac/#memory).

- If you are using Docker Toolbox for Windows or macOS, which uses VirtualBox as backend, so you need to open VirtualBox,
  then stop default VM, select it and click on settings, then make changes as you want.

There isn't any limitation to Docker on Linux system, so don't worry about it.

## Run iMARGI-Docker with Non-root User

root (id = 0) is the default user within a container. It will cause some permission problem of some files or directories
created by Docker container. So it's better to run iMARGI-Docker container using `-u (--user)`  option to override the
default root user with your own user id (UID).

You can use command `id` in your linux system to get your own UID. For example, my UID is `1043`, so I can run iMARGI-Docker
with `-u 1043`, then all the output files and directories are all belong to my user ID. You need to replace the `1043`
with your own UID.

``` bash
id
# uid=1043(frankyan) gid=1048(frankyan)
docker run --rm -u 1043 -v ~/imargi_example:/imargi zhonglab/imargi imargi_wrapper.sh \
    -r hg38 \
    -N HEK_iMARGI \
    -t 16 \
    -g ./ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
    -1 ./data/sample_R1.fastq.gz\
    -2 ./data/sample_R2.fastq.gz \
    -o ./output
```

## Solve `bwa index` Failure Problem

As we created iMARGI-Docker on Linux system, so it works perfectly on Linux system. Generally, it should also work
perfectly on Windows and macOS. However, in our test, we found one critical problem which will cause `baw index`
failure. Here we explain the problem and give a solution to this problem.

_First of all, you need to check you memory. Your computer needs 16 GB memory. If you are using Docker on Windows or
macOS, you also need to set the Docker memory limit to more than 8 GB. If you are sure that the memory is enough, then
try the following solution._

- **Description**: When you use `imargi_wrapper.sh` without `-i` option, the script will generate bwa index files
  automatically. But it might fail when you run it on Windows or macOS.

- **Cause**: Different operating systems have different file system formats, such as NTFS in Windows, APFS in macOS and
  ext4 in Linux. When we use `-v` option to mount host directory to Docker container, it's a kind of map between
  different file system. Most of time, there isn't any problem. However, some tools cannot handle this kind of hybrid
  situation, such as `fsync`, which is utilized by `bwa index` when building large genome index. So `bwa index` fails on
  Windows and macOS.

- **Solution**: Only Windows and macOS users need the solution. The simplest way is use pre-built bwa index files or
  build without Docker. We provide human hg38 bwa index files on our server
  ([link to download](https://sysbio.ucsd.edu/imargi_pipeline/bwa_index.tar.gz)).
  A technical solution to the problem is to build bwa index files in a Docker volume instead of
  a mounted host directory. This solution requires some knowledge of Docker volume. Here we only show demo command
  lines. You need to read Docker documentation if you want to learn more.
  
  ``` bash
  # create Docker volume ref_vol
  docker volume create --name ref_vol
  # start a temporary container to cp ref genome FASTA to Docker volume ref_vol
  docker run -v ref_vol:/data --name helper busybox true
  docker cp ./ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta helper:/data
  docker rm helper
  # run imargi pipeline with two '-v' arguments
  docker run \
    -v /d/imargi_example:/imargi \
    -v ref_vol:/imargi/ref_data \
    zhonglab/imargi \
    imargi_wrapper.sh \
    -r hg38 \
    -N test_sample \
    -t 4 \
    -g ./ref_data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
    -1 ./data/sample_R1.fastq.gz \
    -2 ./data/sample_R2.fastq.gz \
    -o ./output
  ```
  