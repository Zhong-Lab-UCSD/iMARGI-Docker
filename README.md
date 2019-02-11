# iMARGI-Docker

iMARGI-Docker distributes the iMARGI sequencing data processing pipeline

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
    - [Testing Command](#testing-command)
    - [Testing results](#testing-results)
  - [License](#license)
  - [Citation](#citation)

## Description

*in situ* MARGI (**iMARGI**) is a sequencing technique to genome-wide determine the potential genomic interaction loci
of Chromatin associated RNAs (caRNAs). To minimize variations in data processing, we developed a complete data
processing pipeline to improve analysis reproducibility by standardizing data processing steps. **iMARGI-Docker**, a
Docker image, was built to perform the data processing pipeline in a more convenient way.

This repo hosts the iMARGI-Docker source code with brief introductions. For more detail of performing the iMARGI data
analysis using iMARGI-Docker, please read our online comprehensive
[**documentation**](https://sysbio.ucsd.edu/imargi_pipeline).

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
and Mac. For how to install and configure Docker, please read the
[official documentation of Docker](https://docs.docker.com/).

### Installation

#### Pull from Docker Hub

When Docker was installed, it's easy to install iMARGI-Docker by pulling from
[Docker Hub](https://hub.docker.com/r/zhonglab/imargi).

``` bash
docker pull zhonglab/imargi
```

#### Build with Dockerfile

We provided all the source files for building iMARGI-Docker in the [`src`](./src/) folder, including Dockerfile and all
the script tools. So you can modify and rebuild your own Docker image.

## Software Testing Demo

### Testing Data

As real iMARGI sequencing data are always very big, so we randomly extracted a small chunk of real data for software
testing. The data can be found in [data](./data/) folder.

### Testing Command

### Testing results

## License

iMARGI-Docker is licensed under the [BSD 2 license](./src/LICENSE).

## Citation