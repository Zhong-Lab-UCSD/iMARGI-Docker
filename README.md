# iMARGI-Docker

iMARGI-Docker distributes the iMARGI sequencing data processing pipeline

[![License](https://img.shields.io/badge/License-BSD%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)

- [iMARGI-Docker](#imargi-docker)
  - [Description](#description)
  - [Repo Contents](#repo-contents)
  - [Installation Guide](#installation-guide)
    - [Hardware Requirements](#hardware-requirements)
    - [Software Requirements](#software-requirements)
    - [Installation](#installation)
      - [Pull from Docker Hub](#pull-from-docker-hub)
      - [Build with Dockerfile](#build-with-dockerfile)
  - [Software Testing Demo](#software-testing-demo)
    - [Testing Data](#testing-data)
      - [iMARGI sequencing data (paired FASTQ)](#imargi-sequencing-data-paired-fastq)
      - [Reference genome data (FASTA)](#reference-genome-data-fasta)
      - [bwa index data](#bwa-index-data)
    - [Testing Command](#testing-command)
    - [Testing Results](#testing-results)
      - [Running Time Profile](#running-time-profile)
      - [Expected Result files](#expected-result-files)
  - [License](#license)

## Description

*in situ* MARGI (**iMARGI**) is a sequencing technique to genome-wide determine the potential genomic interaction loci
of Chromatin associated RNAs (caRNAs). To minimize variations in data processing, we developed a complete data
processing pipeline to improve analysis reproducibility by standardizing data processing steps. **iMARGI-Docker**, a
Docker image, was built to perform the data processing pipeline in a more convenient way.

This repo hosts the iMARGI-Docker source code with brief introductions. For more detail of performing the iMARGI data
analysis using iMARGI-Docker, please read our online comprehensive
[**documentation**](https://sysbio.ucsd.edu/imargi_pipeline).

We hope every user can perform the iMARGI pipeline with iMARGI-Docker. However, some old machines or operating systems
might not support Docker technique, so users have to configure all the tools used in the pipeline to run it locally.
It can only be done on Linux/MacOS or Windows Subsystem of Linux (WSL), and it's not easy. Please read the
[installation dependencies section of iMARGI pipeline documentation](https://sysbio.ucsd.edu/imargi_pipeline/installation.html)
for details. Users might need root access and have some system administration experience to successfully complete it.

If you encounter any problems, please creating issues in this GitHub repo.

## Repo Contents

- src: source code, such as the Dockerfile of iMARGI-Docker
- data: small chunk of data for testing
- docs: source file of [documentation](https://sysbio.ucsd.edu/imargi_pipeline)

## Installation Guide

### Hardware Requirements

There isn't specific high performance hardware requirements of running iMARGI-Docker. However, as iMARGI generates hugh
amount of sequencing data, usually more than 300 million read pairs, so a high performance computer will save you a lot
of time. Generally, a faster multi-core CPU, larger memory and hard drive storage will benefits you a lot. We suggest
the following specs:

- CPU: at least 4 core CPU
- RAM: at least 16 GB
- Hard drive storage: Depends on your data, typically at least 160 GB is required. Besides, fast IO storage is better,
  such as SSD.

### Software Requirements

iMARGI-Docker only requires Docker. You can use Docker [Community Edition (CE)](https://docs.docker.com/install/) or
[Enterprise Edition (EE)](https://docs.docker.com/ee/). Docker supports all the mainstream OS, such as Linux, Windows
and MacOS. For how to install and configure Docker, please read the
[official documentation of Docker](https://docs.docker.com/).

We recommend using Linux system, because its filesystem is better for large file processing. All the example command
lines here and in the documentation are ran on a Linux system. Most of time, the operations in MacOS is the same as in
Linux system, as it's also a Unix system. However, if you are using Windows system, some command lines need to be
modified. Besides, you need to configure the CPU and memory settings of Docker on Windows system.

### Installation

We recommend pulling the iMARGI-Docker image from Docker Hub. You can also re-build it on your own machine with the
source files in `src` folder.

If you cannot use Docker, please read the
[installation section of iMARGI pipeline documentation](https://sysbio.ucsd.edu/imargi_pipeline/installation.html) for
alternative instructions.

#### Pull from Docker Hub

When Docker was installed, it's easy to install iMARGI-Docker by pulling from
[Docker Hub](https://hub.docker.com/r/zhonglab/imargi). The latest version of iMARGI-Docker image in Docker
Hub is based on the most recent released stable version. It takes about 10 seconds to
install, which depends on your network speed.

``` bash
docker pull zhonglab/imargi
```

#### Build with Dockerfile

Instead of pulling from Docker Hub, you can also build the iMARGI-Docker on your own computer. We provided all the
source code for building iMARGI-Docker in the [`src`](./src/) folder, including Dockerfile and all the script tools.
You can download the most recent stable release or `git clone` from the master branch. So you can modify and rebuild
your own iMARGI-Docker image. It will take about several minutes to build, which depends on your computer performance
and network speed. The stable release is v1.0, which is the master branch.

## Software Testing Demo

To test whether you have successfully installed iMARGI-Docker, you can follow instructions below to do a demo test run.

### Testing Data

#### iMARGI sequencing data (paired FASTQ)

As real iMARGI sequencing data are always very big, so we randomly extracted a small chunk of real data for software
testing. The data can be downloaded from the following links.

- [R1 reads](https://sysbio.ucsd.edu/imargi_pipeline/sample_R1.fastq.gz)
- [R2 reads](https://sysbio.ucsd.edu/imargi_pipeline/sample_R2.fastq.gz)

#### Reference genome data (FASTA)

Besides, you need to download a human genome reference FASTA file.
We use the reference genome used by
[4D Nucleome](https://www.4dnucleome.org/) and
[ENCODE project](https://www.encodeproject.org/data-standards/reference-sequences/).

The FASTA file of the reference
genome is too large for us to host it in GitHub repo. You can be download it use the link:

- [GRCh38_no_alt_analysis_set_GCA_000001405.15](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)

It needs to be decompressed using `gunzip -d` command on Linux/MacOS. If your system is Windows, you can use 7Zip
software to decompress the `.gz` file. Besides, you can also use the `gunzip` tool delivered in iMARGI-Docker.

#### bwa index data

As `bwa index` process will cost a lot of time (more than 1 hour), we suggest to download our pre-built index files for the reference
genome. Please download the following gzip compressed `bwa_index` folder and decompress it (`tar zxvf`) on your machine.

- [bwa index files](https://sysbio.ucsd.edu/imargi_pipeline/bwa_index.tar.gz)

*We assume that you put the data and reference files in the following directory structure.*

``` bash
~/imargi_example
    ├── data
    │   ├── sample_R1.fastq.gz
    │   └── sample_R2.fastq.gz
    ├── output
    └── ref
        ├── GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
        └── bwa_index
            ├── bwa_index_hg38.amb
            ├── bwa_index_hg38.ann
            ├── bwa_index_hg38.bwt
            ├── sample_R1.fastq.pac
            └── sample_R2.fastq.sa
```

### Testing Command

We can use one command line to perform the whole pipeline to the testing data.

``` bash
docker run --rm -u 1043 -v ~/imargi_example:/imargi zhonglab/imargi \
    imargi_wrapper.sh \
    -r hg38 \
    -N test_sample \
    -t 4 \
    -g ./ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
    -i ./ref/bwa_index/bwa_index_hg38 \
    -1 ./data/sample_R1.fastq.gz \
    -2 ./data/sample_R2.fastq.gz \
    -o ./output
```

*Tips:*

- `-u 1043`: Run docker with your own UID of your Linux system (use `id` command to check your UID) to avoid file/dir
  permission problem.

- `-v ~/imargi_example:/imargi`: It mounts the `~/imargi_example` directory in your host machine to workspace of the
  running docker container. The path must be a full path. The example was ran on a Linux computer. If you ran it on a
  Windows  computer, the path is a little different. For example, Windows path `D:\test\imargi_example` needs to be
  rewritten as `/d/test/imargi_example`, so the `-v` argument needs to be `-v /d/test/imargi_example:/imargi`. When you
  executed it on Windows, a window might pop up to verify that you want to share the folder.

- The command line is long, so `\` was used for splitting it into multiple lines in the example. It's a Linux or MacOS
  style. However, in Windows, you need to replace `\` with `^`.

- `-i`: Building bwa index will cost a lot time, so we used the pre-built index files with `-i` argument. There
  are some other arguments can be used for pre-generated files, such as `-R` for restriction fragment BED file and
  `-c` for chromsize file.See more details in the
  [documentation of command line API section](https://sysbio.ucsd.edu/imargi_pipeline/commandline_api.html#imargi-wrapper-sh)

- `-i`: If you don't supply bwa index files, the `imargi_wrapper.sh` will generated     it automatically. It works
  perfectly on Linux system. However, it doesn't work on Windows and MacOS because `bwa index` use `fsync` when build
  large genome index, which cannot handle different driver formats (`-v` mount Windows/MacOS driver to Linux container).
  So it's better to build it in advance. In fact, there's a solution to the problem if you are familiar with Docker
  volume. Please read the
  [technical note of iMARGI pipeline documentation](https://sysbio.ucsd.edu/imargi_pipeline/technical_note.html#solve-bwa-index-failure-problem) for
  detail.

### Testing Results

#### Running Time Profile

It took about 10 minutes to perform the pipeline (with `-i` bwa index argument).

Step | Time | Speed up suggestion
---------|----------|----------
Generating chromosome size file | 10 sec | It's fast, but you can also supply with `-c` once you've generated it before.
Generating bwa index (skipped) | 75 min | Supply with `-i` if you've pre-built index files.
Generating restriction fragment file | 4 min | Supply with `-R` when you've already created it before.
cleaning | 10 sec | It's fast and not parallelization.
bwa mapping | 2 min | More CPU cores with `-t`.
interaction pair parsing | 1 min | More CPU cores with `-t`.

#### Expected Result files

The output result files are in the folder assign with `-o` argument. The final output `.pairs` format file for further
analysis is `final_test_sample.pairs.gz`. Besides, multiple intermediate output files of each step are in the
`clean_fastq`, `bwa_output`, and `parse_temp` sub-directories of the `output` directory. In addition, the generated
chromosome size file, bwa index folder and restriction fragment BED file are all in the `ref` directory, in which the
reference genome FASTA file is. Here is the final directory structure after completing the pipeline.

``` bash
~/imargi_example/
    ├── data
    │   ├── sample_R1.fastq.gz
    │   └── sample_R2.fastq.gz
    ├── ref
    │   ├── GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    │   ├── GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.fai
    │   ├── chromsize.hg38.txt
    │   ├── AluI_frags.bed.gz
    │   ├── AluI_frags.bed.gz.tbi
    │   └── bwa_index
    │       ├── bwa_index_hg38.amb
    │       ├── bwa_index_hg38.ann
    │       ├── bwa_index_hg38.bwt
    │       ├── bwa_index_hg38.pac
    │       └── bwa_index_hg38.sa
    └── output
        ├── bwa_output
        │   ├── bwa_log_test_sample.txt
        │   └── test_sample.bam
        ├── clean_fastq
        │   ├── clean_test_sample_R1.fastq.gz
        │   └── clean_test_sample_R2.fastq.gz
        ├── parse_temp
        │   ├── dedup_test_sample.pairs.gz
        │   ├── drop_test_sample.pairs.gz
        │   ├── duplication_test_sample.pairs.gz
        │   ├── sorted_all_test_sample.pairs.gz
        │   ├── stats_dedup_test_sample.txt
        │   ├── stats_final_test_sample.txt
        │   └── unmapped_test_sample.pairs.gz
        ├── final_test_sample.pairs.gz
        └── pipelineStats_test_sample.log
```

## License

iMARGI-Docker source code is licensed under the [BSD 2 license](./src/LICENSE).