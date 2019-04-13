# Quick Start Example

This quick start example shows how to apply the iMARGI Pipeline using iMARGI-Docker. The all-in-one wrapper script
`imargi_wrapper.sh` is used here to automate the whole pipeline in one command line.
Please go to the [Step-by-step Illustration](./step_by_step_illustration.md) section to learn the details of
the pipeline.

- [Quick Start Example](#quick-start-example)
  - [Pull the iMARGI Docker Image](#pull-the-imargi-docker-image)
  - [Example Data Preparation](#example-data-preparation)
    - [Set Working Directory](#set-working-directory)
    - [Prepare FASTQ Files of Paired-end Sequencing](#prepare-fastq-files-of-paired-end-sequencing)
    - [Download Reference Genome](#download-reference-genome)
    - [Download BWA Index Files](#download-bwa-index-files)
    - [Checklist of Preparation](#checklist-of-preparation)
  - [Perform iMARGI Pipeline in One Command](#perform-imargi-pipeline-in-one-command)
  - [Output of iMARGI Pipeline](#output-of-imargi-pipeline)

## Pull the iMARGI Docker Image

```bash
docker pull zhonglab/imargi
```

## Example Data Preparation

### Set Working Directory

We need to mount a directory on host machine to Docker container as working directory. Here we set the
`~/imargi_example` directory as the mounted working directory, which means we need to use `-v ~/imargi_example:/imargi`
option when we run iMARGI-Docker. We created some sub-directories for organizing data and results.

``` bash
cd ~/imargi_example
mkdir ref data output
```

Then the structure of working directory is:

``` bash
~/imargi_example
    ├── ref
    │   └──
    ├── data
    │   └──
    └── output
        └──
```

### Prepare FASTQ Files of Paired-end Sequencing

Here we take a small demo dataset as the example. The demo dataset can be processed in about 15 mins.
Here are the download links:

- [Click Here to Download: R1 reads fastq.gz](https://sysbio.ucsd.edu/imargi_pipeline/sample_R1.fastq.gz)

- [Click Here to Download: R2 reads fastq.gz](https://sysbio.ucsd.edu/imargi_pipeline/sample_R2.fastq.gz)

The demo dataset is randomly sampled from a real HiSeq iMARGI data in HEK cells generated in our previous published
paper <a id="a1">[[1]](#f1)</a>. If you are interested in the real dataset, you can get it in NCBI SRA database with
accession number [SRR8206679](https://www.ncbi.nlm.nih.gov/sra?term=SRX5026004). It's a very large dataset
(more than 40 GB) and processing it will take a lot of time (more than 10 hours).

### Download Reference Genome

The minimum required reference file is only a reference genome sequence FASTA file.

We use the same reference genome used by
[4D Nucleome](https://www.4dnucleome.org/) and
[ENCODE project](https://www.encodeproject.org/data-standards/reference-sequences/). The FASTA file of the reference
genome can be downloaded use the link:

[Click Here to Download: GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz).

It needs to be decompressed using `gunzip -d` command. The command lines for preparing the reference genome files:

``` bash
cd ~/imargi_example

wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
gunzip -d GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
mv GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta ./ref
```

### Download BWA Index Files

Although `imargi_wrapper.sh` can automatically build bwa index files, it will cost a lot of time (more than 1 hour
for human genome). Hence, we suggest to download our pre-built BWA index files for the human hg38 reference
genome and supply it to `imargi_wrapper.sh` with `-i` option. Please download the following gzip compressed `bwa_index`
folder and decompress it (`tar zxvf`) in your machine.

- [Click Here to Download: bwa index files](https://sysbio.ucsd.edu/imargi_pipeline/bwa_index.tar.gz)

``` bash
cd ~/imargi_example/ref
wget  https://sysbio.ucsd.edu/imargi_pipeline/bwa_index.tar.gz
tar zxvf bwa_index.tar.gz
```

### Checklist of Preparation

After preparation, you will get all the required data and reference files in the working directory. It's the minimum
requirements, which are sequencing reads data in FASTQ format and reference genome sequence in FASTA format.

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

## Perform iMARGI Pipeline in One Command

The easiest way of is using the all-in-one wrapper script `imargi_wrapper.sh` with default parameters.
We used `-N test_sample` argument to set base identifier for all the result files.

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

> *Tips:*
> 
> - `--rm`: By default a container’s file system persists even after the container exits. Hence, the container file
> systems can really pile up. `--rm` option will automatically clean up the container after the container exits.
> 
> - `-t`: Allocate a pseudo-TTY. Without `-t`, you cannot use `Ctrl + c` to stop the run.
> 
> - `-u 1043`: Run docker with your own UID of your Linux system (use `id` command to check your own UID and replace
>   `1043` with it) to avoid file/dir permission problem.
> 
> - `-v ~/imargi_example:/imargi`: It mounts the `~/imargi_example` directory in your host machine to workspace of the
>   running docker container. The path must be a full path. The example was ran on a Linux computer. If you ran it on a
>   Windows  computer, the path is a little different. For example, Windows path `D:\test\imargi_example` needs to be
>   rewritten as `/d/test/imargi_example`, so the `-v` argument needs to be `-v /d/test/imargi_example:/imargi`. When you
>   executed it on Windows, a window might pop up to verify that you want to share the folder.
> 
> - The command line is long, so `\` was used for splitting it into multiple lines in the example. It's a Linux or MacOS
>   style. However, in Windows, you might need to replace `\` with `^`.
> - `-i`: Building bwa index will cost a lot time, so we used the pre-built index files with `-i` argument. There
>   are also some other arguments can be used for pre-generated files, such as `-R` for restriction fragment BED file and
>   `-c` for chromsize file.See more details in the
>   [documentation of command line API section](https://sysbio.ucsd.edu/imargi_pipeline/commandline_api.html#imargi-wrapper-sh)
> 
> - `-i`: If you don't supply bwa index files, the `imargi_wrapper.sh` will generated     it automatically. It works
>   perfectly on Linux system. However, it doesn't work on Windows and MacOS because `bwa index` use `fsync` when build
>   large genome index, which cannot handle different driver formats (`-v` mount Windows/MacOS driver to Linux container).
>   So it's better to build it in advance. In fact, there's a solution to the problem if you are familiar with Docker
>   volume. Please read the
>   [technical note of iMARGI pipeline documentation](https://sysbio.ucsd.edu/imargi_pipeline/technical_note.html#solve-bwa-index-failure-problem) for
>   detail.

It will take about 10 minutes to complete the whole processing pipeline to the demo dataset on our computer
(with `-i` bwa index argument).

Step | Time | Speed up suggestion
---------|----------|----------
Generating chromosome size file | 10 sec | It's fast, but you can also supply with `-c` once you've generated it before.
Generating bwa index (skipped) | 75 min | Supply with `-i` if you've pre-built index files.
Generating restriction fragment file | 4 min | Supply with `-R` when you've already created it before.
cleaning | 10 sec | It's fast and not parallelization.
bwa mapping | 2 min | More CPU cores with `-t`.
interaction pair parsing | 1 min | More CPU cores with `-t`.

## Output of iMARGI Pipeline

Once the iMARGI Pipeline has completed, all the result files are in the output directory. The final tree structure of the
whole working directory is shown below. The final RNA-DNA interaction map is in the `final_test_sample.pairs.gz` file,
which is a compressed `.pairs` format file and can be used for further analysis. All the intermediate result files are
also kept in several sub-directories. As we used the minimum input requirements, so the pipeline automatically
generated several new reference files in the same directory of reference genome, including chromosome sizes, bwa index
and AluI digestion fragments. When processing new dataset, you can reuse these new generated reference files with
corresponding arguments instead of only using `-g` argument, which will save you some time and disk space.
Besides, there is a simple stats file, `pipelineStats_test_sample.log`, which reports the total processed read pairs
number, BWA mapping stats and number of valid RNA-DNA interaction in the final `.pairs.gz` file. For more detail, please
check the [Step-by-step Illustration](./step_by_step_illustration.md) and
[Command-line API](./commandline_api.md) sections.

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

**Reference:**

<small>[[1]](#a1) <span id="f1"></span> Yan, Z. et al. Genome-wide co-localization of RNA-DNA interactions and fusion RNA pairs. PNAS February 19, 2019, 116 (8) 3328-3337. https://doi.org/10.1073/pnas.1819788116</small>