# iMARGI-Docker

iMARGI-Docker distributes the iMARGI sequencing data processing pipeline

[![License](https://img.shields.io/badge/License-BSD%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)

- [iMARGI-Docker](#imargi-docker)
  - [1. Description](#1-description)
  - [2. Repository Contents](#2-repository-contents)
  - [3. Installation Guide](#3-installation-guide)
    - [3.1. Hardware Requirements](#31-hardware-requirements)
    - [3.2. Software Requirements](#32-software-requirements)
      - [3.2.1. Docker Installation](#321-docker-installation)
      - [3.2.2. Docker settings (macOS or Windows)](#322-docker-settings-macos-or-windows)
    - [3.3. iMARGI-Docker Installation](#33-imargi-docker-installation)
      - [3.3.1. Pull from Docker Hub](#331-pull-from-docker-hub)
      - [3.3.2. Build with Dockerfile](#332-build-with-dockerfile)
  - [4. Software Testing Demo](#4-software-testing-demo)
    - [4.1. Testing Data](#41-testing-data)
      - [4.1.1. iMARGI sequencing data (paired FASTQ)](#411-imargi-sequencing-data-paired-fastq)
      - [4.1.2. Reference genome data (FASTA)](#412-reference-genome-data-fasta)
      - [4.1.3. bwa index data](#413-bwa-index-data)
    - [4.2. Testing Command](#42-testing-command)
    - [4.3. Testing Results](#43-testing-results)
      - [4.3.1. Running Time Profile](#431-running-time-profile)
      - [4.3.2. Expected Result files](#432-expected-result-files)
  - [License](#license)

## 1. Description

*in situ* MARGI (**iMARGI**) is a sequencing technique to genome-wide determine the potential genomic interaction loci
of Chromatin associated RNAs (caRNAs). To minimize variations in data processing, we developed a complete data
processing pipeline to improve analysis reproducibility by standardizing data processing steps. **iMARGI-Docker**, a
Docker image, was built to perform the data processing pipeline in a more convenient way.

This repo hosts the iMARGI-Docker source code with brief introductions. For more detail of performing the iMARGI data
analysis using iMARGI-Docker, please read our online comprehensive
[**documentation**](https://sysbio.ucsd.edu/imargi_pipeline).

We hope every user can perform the iMARGI pipeline with iMARGI-Docker. However, some old machines or operating systems
might not support Docker technique, so users have to configure all the tools used in the pipeline to run it locally.
It can only be done on Linux/macOS, and it requires solid experience in Linux system configuration. Please read the
[installation dependencies section of iMARGI pipeline documentation](https://sysbio.ucsd.edu/imargi_pipeline/installation.html)
for details.

If you encounter any problems, please create issues in this GitHub repo.

## 2. Repository Contents

- src: source code, such as the Dockerfile of iMARGI-Docker
- data: small chunk of data for testing
- docs: source file of [documentation](https://sysbio.ucsd.edu/imargi_pipeline)

## 3. Installation Guide

### 3.1. Hardware Requirements

There isn't specific high performance hardware requirements of running iMARGI-Docker. However, as iMARGI generates hugh
amount of sequencing data, usually more than 300 million read pairs, so a high performance computer will save you a lot
of time. Generally, a faster multi-core CPU, larger memory and hard drive storage will benefits you a lot. We suggest
the following specs:

- CPU: At least dual core CPU. More CPU cores will speed up the processing.
  
- RAM: 16 GB. Depends on the size of reference genome. For human genome, at least 8GB free memory are required by BWA,
  so the memory on the machine needs to be more than 8 GB, which usually is 16 GB. **Out of memory will cause ERROR.**

- Hard drive storage: Depends on your data, typically at least 160 GB free space is required for 300M 2x100 read pairs.
  Besides, fast IO storage is better, such as SSD.

### 3.2. Software Requirements

iMARGI-Docker only requires Docker installed on your computer.
You can use [Docker Community Edition (CE)](https://docs.docker.com/install/).

Although Docker supports all the mainstream OS, such as Linux, Windows and macOS,
**we strongly recommend using Linux system**, because it's much easier to setup and its filesystem is better for large
file processing. You can install Docker CE with only two commands on well supported 64-bit Linux distributions, including
Ubuntu, Debian, Fedora, and CentOS.

Keep in mind, all the example command lines here and in the documentation are ran on a Linux system.
Most of time, the operations in macOS is the same as in Linux system, as it's also a Unix system. However, if you are
using Windows system, some command lines need to be modified. Besides, you need to do additional configurations
of Docker on Windows or macOS system.

#### 3.2.1. Docker Installation

Here are some essential instructions for installing Docker on different systems. Install Docker on Linux is the easiest.

- **Linux**: Support the most recent 64 bit stable releases of Ubuntu, Debian, Fedora and CentOS. You need `root` or `sudo`
  privileges. Generally, the following commands will automatically install Docker in your system. The second command
  will allow you to run `docker` commands without `sudo` privileges.
  [Learn more from the official documentation.](https://docs.docker.com/install/linux/docker-ce/ubuntu/)
  
  ``` Bash
  # install Docker on Ubuntu, Debian, Fedora, and CentOS
  sudo curl -fsSL https://get.docker.com |sh -
  
  # set Docker user, replace demo_user with you own user name,
  # then you can use docker command without sudo
  sudo usermod -aG docker demo_user
  ```

- **macOS (modern)**: Docker Desktop for macOS. Support macOS Sierra 10.12 and newer on a Apple computer after 2010.
  
  Download Docker Desktop software for macOS and install.
  [Click here to check instructions](https://docs.docker.com/docker-for-mac/install/)

- **Windows 10 (modern)**: Docker Desktop for Windows. Support the latest Windows 10 (64 bit) Pro, Enterprise or
  Education version.
  
  First, enable virtualization of your CPU (most of modern Intel CPUs support virtualization).
  [Check here to see how to enable it in BIOS.](https://www.isumsoft.com/computer/enable-virtualization-technology-vt-x-in-bios-or-uefi.html)
  
  Then, turn on Hyper-V. [Check here to see how to turn on Hyper-V.](https://docs.microsoft.com/en-us/virtualization/hyper-v-on-windows/quick-start/enable-hyper-v)
  
  Finally, download Docker Desktop software for Windows and install,
  [Click here to check instructions](https://docs.docker.com/docker-for-windows/install/)

- **Legacy solution**: For older Mac and Windows systems that do not meet the requirements of Docker Desktop for Mac and
  Docker Desktop for Windows, you can install Docker Toolbox to use Docker.

  First, enable virtualization of your CPU (most of modern Intel CPUs support virtualization).
  [Check here to see how to enable it in BIOS.](https://www.isumsoft.com/computer/enable-virtualization-technology-vt-x-in-bios-or-uefi.html)
  
  [Download the latest version of Docker Toolbox from GitHub repo](https://github.com/docker/toolbox/releases)
  
  [Instructions of Docker Toolbox for Windows](https://docs.docker.com/toolbox/toolbox_install_windows/)
  
  [Instructions of Docker Toolbox for macOS](https://docs.docker.com/toolbox/toolbox_install_mac/)

If you are using macOS or Windows, you can check the
[Technical Notes of installing Docker on different systems](https://sysbio.ucsd.edu/imargi_pipeline/technical_note.html#install-docker-on-different-systems)
to learn how to install Docker on other systems.

#### 3.2.2. Docker settings (macOS or Windows)

For macOS and Windows, there are CPU and memory limitations to Docker, which are 1 CPU core and 2 GB memory as default.
The memory is far from the requirement of BWA for human genome, which will cause ERROR. So it must be changed to more
than 8 GB memory. If you have 4 CPU cores, it's better to increase the CPU limitation.

Here are simple instructions of how to change the settings.

- If you are using Docker Desktop for Windows or macOS, you can easily change the settings by right click the
  Docker icon (Whale) in the task bar, then go to Settings -> Advanced to change memory and CPU limits.
  More detail can be found in the Docker official docs of
  [Get started with Docker for Windows](https://docs.docker.com/docker-for-windows/), and
  [Get started with Docker Desktop for Mac](https://docs.docker.com/docker-for-mac/#memory).

- If you are using Docker Toolbox for Windows or macOS, which uses VirtualBox as backend, so you need to open VirtualBox,
  then stop default VM, Select it and click on settings, then make changes as you want.

There isn't any limitation to Docker on Linux system, so don't worry about it.

### 3.3. iMARGI-Docker Installation

We recommend pulling the iMARGI-Docker image from Docker Hub. You can also re-build it on your own machine with the
source files in `src` folder.

If you cannot use Docker, please read the
[installation section of iMARGI pipeline documentation](https://sysbio.ucsd.edu/imargi_pipeline/installation.html) for
alternative instructions.

#### 3.3.1. Pull from Docker Hub

When Docker was installed, it's easy to install iMARGI-Docker by pulling from
[Docker Hub](https://hub.docker.com/r/zhonglab/imargi). The latest version of iMARGI-Docker image in Docker
Hub is based on the most recent released stable version. It takes about 10 seconds to
install, which depends on your network speed.

``` bash
docker pull zhonglab/imargi
```

#### 3.3.2. Build with Dockerfile

Instead of pulling from Docker Hub, you can also build the iMARGI-Docker on your own computer. We provided all the
source code for building iMARGI-Docker in the [`src`](./src/) folder, including Dockerfile and all the script tools.
You can download the most recent stable release or `git clone` from the master branch. So you can modify and rebuild
your own iMARGI-Docker image. It will take about several minutes to build, which depends on your computer performance
and network speed. Currently, the stable release is v1.0, which is the master branch.

## 4. Software Testing Demo

To test whether you have successfully deployed iMARGI-Docker, you can follow instructions below to do a demo test run.

### 4.1. Testing Data

#### 4.1.1. iMARGI sequencing data (paired FASTQ)

As real iMARGI sequencing data are always very big, so we randomly extracted a small chunk of real data for software
testing. The data can be downloaded from the following links.

- [R1 reads](https://sysbio.ucsd.edu/imargi_pipeline/sample_R1.fastq.gz)
- [R2 reads](https://sysbio.ucsd.edu/imargi_pipeline/sample_R2.fastq.gz)

#### 4.1.2. Reference genome data (FASTA)

Besides, you need to download a human genome reference FASTA file.
We use the reference genome used by
[4D Nucleome](https://www.4dnucleome.org/) and
[ENCODE project](https://www.encodeproject.org/data-standards/reference-sequences/).

The FASTA file of the reference
genome is too large for us to host it in GitHub repo. You can be download it use the link:

- [GRCh38_no_alt_analysis_set_GCA_000001405.15](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz)

**It needs to be decompressed** using `gunzip -d` or `gzip -d` command on Linux/macOS. If your system is Windows, you can
use 7Zip or other software to decompress the `.gz` file. Besides, you can also use the `gunzip` tool delivered in iMARGI-Docker.

#### 4.1.3. bwa index data

As `bwa index` process will cost a lot of time (more than 1 hour), we suggest to download our pre-built index files for
the reference genome. Please download the following gzip compressed `bwa_index` folder and decompress it
(`tar zxvf`) on your machine.

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
            ├── bwa_index_hg38.pac
            └── bwa_index_hg38.sa
```

### 4.2. Testing Command

We can use one command line to perform the whole pipeline to the testing data.

``` bash
cd ~/imargi_example

# replace "-u 1043" with your own UID, see the tips below
# replace "-v ~/imargi_example:/imargi" with your working directory if not ~/imargi_example

docker run --rm -t -u 1043 -v ~/imargi_example:/imargi zhonglab/imargi \
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
- `--rm`: By default a container’s file system persists even after the container exits. Hence, the container file
  systems can really pile up. `--rm` option will automatically clean up the container after the container exits.

- `-t`: Allocate a pseudo-TTY. With `-t`, you can use `Ctrl + c` to stop the run in terminal. However, without `-t`, you
  have to use `docker ps -a` to check your run container id and then use `docker stop <container_id>` to stop the run.

- `-u 1043`: Run docker with your own UID of your Linux system (use `id` command to check your UID) to avoid file/dir
  permission problem.

- `-v ~/imargi_example:/imargi`: It mounts the `~/imargi_example` directory in your host machine to workspace of the
  running docker container. The path must be a full path. The example was ran on a Linux computer. If you ran it on a
  Windows  computer, the path is a little different. For example, Windows path `D:\test\imargi_example` needs to be
  rewritten as `/d/test/imargi_example`, so the `-v` argument needs to be `-v /d/test/imargi_example:/imargi`. When you
  executed it on Windows, a window might pop up to verify that you want to share the folder.

- The command line is long, so `\` was used for splitting it into multiple lines in the example. It's a Linux or macOS
  style. However, in Windows, you need to replace `\` with `^`.

- `-i`: Building BWA index will cost a lot time, so we used the pre-built index files with `-i` argument. If you don't
  supply BWA index files, the `imargi_wrapper.sh` will generated it automatically based on the reference genome sequence
  supplied by `-g` parameter. Building BWA index needs large memory as we required (16 GB). There
  are some other arguments can be used for pre-generated files, such as `-R` for restriction fragment BED file
  (the automatically generated file is named as `AluI_frags.bed.gz`) and `-c` for chromsize file. See more details in the
  [documentation of command line API section](https://sysbio.ucsd.edu/imargi_pipeline/commandline_api.html#imargi-wrapper-sh)

### 4.3. Testing Results

#### 4.3.1. Running Time Profile

It took about 10 minutes to perform the pipeline on our computer (with `-i` bwa index argument).

Step | Time | Speed up suggestion
---------|----------|----------
Generating chromosome size file | 10 sec | It's fast, but you can also supply with `-c` once you've generated it before.
Generating bwa index (skipped) | 75 min | Supply with `-i` if you've pre-built index files.
Generating restriction fragment file | 4 min | Supply with `-R` when you've already created it before.
cleaning | 10 sec | It's fast and not parallelization.
bwa mapping | 2 min | More CPU cores with `-t`.
interaction pair parsing | 1 min | More CPU cores with `-t`.

#### 4.3.2. Expected Result files

The output result files are in the folder assign with `-o` argument. The final output `.pairs` format file for further
analysis is `final_test_sample.pairs.gz`. Besides, multiple intermediate output files of each step are in the
`clean_fastq`, `bwa_output`, and `parse_temp` sub-directories of the `output` directory. In addition, the generated
chromosome size file, bwa index folder and restriction fragment BED file are all in the `ref` directory, in which the
reference genome FASTA file is. Besides, there is also a simple stats file, `pipelineStats_test_sample.log`, which reports the
total processed read pairs number, BWA mapping stats and number of valid RNA-DNA interaction in the
final `.pairs.gz` file. Here is the final directory structure after completing the pipeline.

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