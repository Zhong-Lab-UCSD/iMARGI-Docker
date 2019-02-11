# Tools and Installation

In the iMARGI pipeline, a bundle of tools and several bash scripts created by us are required. For convenience,
we built [iMARGI-Docker](https://hub.docker.com/r/zhonglab/imargi) to distribute the pipeline. It delivers all the
configured tools.

- [Tools and Installation](#tools-and-installation)
  - [Docker Container Usage Instructions](#docker-container-usage-instructions)
  - [Dependencies Instruction](#dependencies-instruction)

## Docker Container Usage Instructions

An iMARGI Docker image is available in [Docker-Hub](https://hub.docker.com/r/zhonglab/imargi). It's
much easier to apply the iMARGI pipeline using the docker container than installing and configuring all the required
tools.

Docker Community Edition or Enterprise Edition is required to be installed and configured on the Linux server (root
authority is required for installing Docker CE).
[Here is the guides of installing Docker CE on Ubuntu](https://docs.docker.com/install/linux/docker-ce/ubuntu/).

After Docker CE was installed, you can pull the latest iMARGI Docker image to your sever.

```bash
docker pull zhonglab/imargi
```

To use the tools in the iMARGI Docker image, you need to run a Docker container. Here is an example of creating a
new directory with `mkdir` command through a Docker container.

``` bash
docker run -v ~/test:/imargi imargi mkdir new_dir
```

![](./figures/docker_command_example.png)

You should know the `-v` or `--volume` option, which assigns the `~/test` directory on your host machine to the
working directory `/imargi` of the iMARGI container. The container can only operate the files in `~/test` directory.
In the example, a folder `new_dir` will be created in the folder `~/test/`. You can just replace the command to use
any tool in the iMARGI Docker container.

## Dependencies Instruction

All the dependencies are listed here. However, we strongly recommend using iMARGI-Docker instead of configuring these
dependencies on your Linux server. The following table shows all those required tools with simple descriptions. Some
of these tools, such as `bash`, `sort` and `zcat` are usually default installed in most of Linux distributions. Besides,
you might need root access or compiling tools on Linux system to install some of these tools.

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
cooler | 0.8.2 | [Following instruction](https://github.com/mirnylab/cooler)| Utilities for genomic interaction data
SRA Toolkit | 2.9.4  | [Following instruction](https://github.com/ncbi/sra-tools) | NCBI SRA tools
GNU parallel | - | Linux package "parallel" | Executing jobs in parallel
[GNU awk](https://www.gnu.org/software/gawk/manual/html_node/Quick-Installation.html)| - | Linux package "gawk" | Text file processing tool
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
imargi_annotate.sh | 0.1| Download and `chmod +x` | Annotate RNA/DNA-ends with gene annotations