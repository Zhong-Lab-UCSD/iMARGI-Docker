# Tools and Installation

In the iMARGI pipeline, a bundle of tools and several bash scripts created by us are required. For convenience and
reproducibility, we built [iMARGI-Docker](https://hub.docker.com/r/zhonglab/imargi) to distribute the pipeline.
It delivers all the well configured tools.

- [Tools and Installation](#tools-and-installation)
  - [Docker Container Usage Instructions](#docker-container-usage-instructions)
  - [Dependencies Instruction](#dependencies-instruction)

## Docker Container Usage Instructions

An iMARGI Docker image is available in [Docker-Hub](https://hub.docker.com/r/zhonglab/imargi), and its source files are
hosted in [iMARGI-Docker GitHub repo](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker). It's much easier to apply the
iMARGI pipeline using the docker container than installing and configuring all the required tools.

Docker Community Edition or Enterprise Edition is required to be installed and configured on the machine (root
authority is required for installing Docker CE).
[Here is the guides of installing Docker CE on Ubuntu](https://docs.docker.com/install/linux/docker-ce/ubuntu/).

After Docker CE was installed, you can pull the latest iMARGI Docker image to your sever.

```bash
docker pull zhonglab/imargi
```

To use the tools in the iMARGI Docker image, you need to run a Docker container. Here is an example of creating a
new directory with `mkdir` command through a Docker container.

``` bash
docker run -v ~/test:/imargi zhonglab/imargi mkdir new_dir
```

![](./figures/docker_command_example.png)

You should know the `-v` or `--volume` option, which assigns the `~/test` directory on your host machine to the
working directory `/imargi` of the iMARGI container. The container can only operate the files in `~/test` directory.
In the example, a folder `new_dir` will be created in the folder `~/test/`. You can just replace the command to use
any tool in the iMARGI Docker container. Besides, there are many other options are useful. For example, you can use
`--rm` to automatically clean up the container after finished its job. For more usage information of Docker, please
refer to [Docker documentation](https://docs.docker.com/engine/reference/commandline/cli/).

## Dependencies Instruction

We strongly recommend using iMARGI-Docker instead of configuring all the dependencies of iMARGI pipeline on your Linux
server. However, if you really cannot run Docker on your machine, you might want to try to configure these tools. It
requires root access to your machine and solid experience of Linux server administration.

We cannot guarantee success of local configuration. If you encounter some problems or have suggestions, please create
issues in the [iMARGI-Docker GitHub repo](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker). If you are using Ubuntu
(18.04), the following command lines we used to configure iMARGI-Docker might help you.

``` bash
# run with root account
apt-get update
apt-get install git build-essential libz-dev libbz2-dev liblzma-dev libssl-dev libcurl4-gnutls-dev \
    autoconf automake libncurses5-dev wget gawk parallel
cd /tmp && git clone https://github.com/lh3/seqtk.git && \
    cd seqtk && make && make install
cd /tmp && git clone https://github.com/samtools/htslib && \
    cd htslib && autoheader && autoconf && \
    ./configure --prefix=/usr/local && make && make install
cd /tmp && git clone https://github.com/samtools/samtools && \
    cd samtools && autoheader && autoconf && \
    ./configure --prefix=/usr/local && make && make install
cd /tmp && git clone https://github.com/lh3/bwa.git && \
    cd bwa && make && cp bwa /usr/local/bin
cd /tmp && git clone https://github.com/nh13/pbgzip && \
    cd pbgzip && sh autogen.sh && ./configure && make && make install
cd /tmp && git clone https://github.com/lz4/lz4 && \
    cd lz4 && make && make install
cd /tmp && wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.4/sratoolkit.2.9.4-ubuntu64.tar.gz && \
    tar zxvf sratoolkit.2.9.4-ubuntu64.tar.gz && cp -R sratoolkit.2.9.4-ubuntu64/bin/* /usr/local/bin

apt-get install -y python3-dev libopenblas-dev python3-pip
pip3 install numpy cython scipy pandas click
pip3 install pairtools cooler HTSeq

cd /tmp && git clone https://github.com/Zhong-Lab-UCSD/iMARGI-Docker.git &&\
    cp iMARGI-Docker/src/imargi_* /usr/local/bin/

chmod +x /usr/local/bin/imargi_*
```

The following table shows all the required tools with simple descriptions. Some of these tools, such as `bash`, `sort`
and `zcat` are usually default installed in most of Linux distributions. Besides, you might need root access or
compiling tools on Linux system to install some of these tools.

Tool | Version  | Installation | Brief description
---------|----------|---------|-----------
Python | 3.x | [Following instruction](https://www.python.org/downloads/) | Running environment for several tools
seqtk | 1.3 | [Following instruction](https://github.com/lh3/seqtk)| Processing FASTA/FASTQ files
bwa | 0.7.17 | [Following instruction](https://github.com/lh3/bwa) | Mapping reads to reference genome
samtools | 1.9 | [Following instruction](http://www.htslib.org/download/)| Manipulating SAM/BAM files
htslib | 1.9 | [Following instruction](http://www.htslib.org/download/)| Manipulating SAM/BAM files
pairtools | 0.2.2 | [Following instruction](https://pairtools.readthedocs.io/en/latest/installation.html)| Utilities for processing interaction pairs
lz4 | 1.8.3 | [Following instruction](https://github.com/lz4/lz4) | Extremely fast compression
pbgzip | - | [Following instruction](https://github.com/nh13/pbgzip)| Compression for Genomics Data
cooler | 0.8.3 | [Following instruction](https://github.com/mirnylab/cooler)| Utilities for genomic interaction data
HTSeq | 0.11.2 | [Following instruction](https://htseq.readthedocs.io/en/master/install.html)| Utilities for annotating interactions
SRA Toolkit | 2.9.4  | [Following instruction](https://github.com/ncbi/sra-tools) | NCBI SRA tools
GNU parallel | - | Linux package "parallel" | Executing jobs in parallel
[GNU awk](https://www.gnu.org/software/gawk/manual/html_node/Quick-Installation.html)| - | Linux package "gawk", set alias awk | Text file processing tool
bash | - | Linux package "bash" | Shell environment
sort | - | Linux package "sort" | Sort text
gunzip | - | Linux package "gunzip" | Compression tool
zcat | - | Linux package "zcat" | Readout compressed text file
imargi_wrapper.sh | 0.1 | Download and `chmod +x` | All-in-one pipeline wrapper
imargi_clean.sh | 0.1 | Download and `chmod +x` | Clean iMARGI paired end fastq files
imargi_parse.sh | 0.1| Download and `chmod +x` | Parse BAM to valid RNA-DNA interaction pairs
imargi_stats.sh | 0.1 | Download and `chmod +x` | Simple stats report of .pairs file
imargi_convert.sh | 0.1| Download and `chmod +x` | Convert .pairs format to other formats
imargi_distfilter.sh | 0.1 | Download and `chmod +x` | Filter .pairs or BEDPE file with genomic distance threshold
imargi_rsfrags.sh | 0.1 | Download and `chmod +x` | Generate restriction fragment BED file
imargi_restrict.py | 0.1 | Download and `chmod +x` | Restriction site analysis of .pairs file
imargi_annotate.sh | 0.1| Download and `chmod +x` | Annotate RNA/DNA-ends with genomic annotations
imargi_ant.py | 0.1 | Download and `chmod +x` | Annotate RNA/DNA-ends with genomic annotations, used by imargi_annotate.sh