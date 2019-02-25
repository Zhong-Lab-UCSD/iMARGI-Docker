# Quick Start Example

This quick start example shows how to apply the iMARGI Pipeline to a real iMARGI sequencing dataset. The
all-in-one wrapper script `imargi_wrapper.sh` is used here to automate the whole pipeline in one command line.
Please go to the [Step-by-step Illustration](./step_by_step_illustration.md) section to learn the details of
the pipeline.

- [Quick Start Example](#quick-start-example)
  - [Pull the iMARGI Docker Image](#pull-the-imargi-docker-image)
  - [Example Data Preparation](#example-data-preparation)
    - [Set Working Directory](#set-working-directory)
    - [Prepare FASTQ Files of Paired-end Sequencing](#prepare-fastq-files-of-paired-end-sequencing)
    - [Download Reference Genome](#download-reference-genome)
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
option when we run docker container. We created some sub-directories for organizing data and results.

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

Here we take the HEK (HEK293T) iMARGI data set generated in our previous published paper <a id="a1">[[1]](#f1)</a> as the
example. The relatd data can be accessed in GEO database by
[GSM3478205](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3478205). The raw sequencing reads data is released
in NCBI SRA database with accession number [SRR8206679](https://www.ncbi.nlm.nih.gov/sra?term=SRX5026004). So it can be
downloaded and prepared by fastq-dump tool. Then we can get the iMARGI sequencing data in a pair of compressed FASTQ
format files, `SRR8206679_1.fastq.gz` and `SRR8206679_2.fastq.gz`.

``` bash
docker run -v ~/imargi_example:/imargi imargi fastq-dump --gzip --split-3 SRR8206679

cd ~/imargi_example
mv ./SRR8206679_*.fastq.gz ./data
```

### Download Reference Genome

The required reference files include a reference genome sequence FASTA file and a chromosome size text file.

We use the same reference genome used by
[4D Nucleome](https://www.4dnucleome.org/) and
[ENCODE project](https://www.encodeproject.org/data-standards/reference-sequences/). The FASTA file of the reference
genome can be downloaded use the link:
[GRCh38_no_alt_analysis_set_GCA_000001405.15](https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz).
It needs to be decompressed using `gunzip -d` command.

The command lines for preparing the reference genome files:

``` bash
cd ~/imargi_example

wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
gunzip -d GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
mv GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta ./ref
```

### Checklist of Preparation

After preparation, you will get all the required data and reference files in the working directory. It's the minimum
requirements, which are sequencing reads data in FASTQ format and reference genome sequence in FASTA format.

``` bash
~/imargi_example
    ├── data
    │   ├── SRR8206679_1.fastq.gz
    │   └── SRR8206679_2.fastq.gz
    ├── output
    └── ref
        └── GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
```

## Perform iMARGI Pipeline in One Command

The easiest way of is using the all-in-one wrapper script `imargi_wrapper.sh` with default parameters.
We used ‘-N HEK_iMARGI’ argument to set base name for all the result files.

``` bash
cd ~/imargi_example
mkdir ./output
docker run -v ~/imargi_example:/imargi imargi imargi_wrapper.sh \
    -r hg38 \
    -N HEK_iMARGI \
    -t 16 \
    -g ./ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
    -1 ./data/SRR8206679_1.fastq.gz \
    -2 ./data/SRR8206679_1.fastq.gz \
    -o ./output
```

As the size of sequencing data is very large (more than 350 million read pairs), it will cost about 12 hours to process.

## Output of iMARGI Pipeline

Once the iMARGI Pipeline has completed, all the result files are in the output directory. The tree structure of the
whole working directory is shown below. The final RNA-DNA interaction map is in the `final_HEK_iMARGI.pairs` file,
which is a compressed .pairs format file and can be used for further analysis. All the intermediate result files are
also kept in several sub-directories. As we used the minimum input requirements, so the pipeline automatically
generated several new reference files in the same directory of reference genome, including chromosome sizes, bwa index
and AluI digestion fragments. When processing new dataset, you can reuse these new generated reference files with
corresponding arguments instead of only using `-g` argument, which will save you some time and disk space. For more
detail, please check the [Step-by-step Illustration](./step_by_step_illustration.md) and
[Command-line API](./commandline_api.md) sections.

``` bash
~/imargi_example/
    ├── data
    │   ├── SRR8206679_1.fastq.gz
    │   └── SRR8206679_1.fastq.gz
    ├── ref
    │   ├── GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    │   ├── GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.fai
    │   ├── chrom.sizes.hg38.txt
    │   ├── AluI_frags.bed.gz
    │   └── bwa_index
    │       ├── bwa_index_hg38.amb
    │       ├── bwa_index_hg38.ann
    │       ├── bwa_index_hg38.bwt
    │       ├── bwa_index_hg38.pac
    │       └── bwa_index_hg38.sa
    └── output
        ├── bwa_output
        │   └── HEK_iMARGI.bam
        ├── clean_fastq
        │   ├── clean_SRR8206679_1.fastq.gz
        │   └── clean_SRR8206679_2.fastq.gz
        ├── parse_temp
        │   ├── dedup_HEK_iMARGI.pairs.gz
        │   ├── drop_HEK_iMARGI.pairs.gz
        │   ├── duplication_HEK_iMARGI.pairs.gz
        │   ├── sorted_all_HEK_iMARGI.pairs.gz
        │   ├── stats_dedup_HEK_iMARGI.txt
        |   ├── stats_final_HEK_iMARGI.txt
        │   └── unmapped_HEK_iMARGI.pairs.gz
        ├── final_HEK_iMARGI.pairs.gz
        └── pipelineStats_HEK_iMARGI.log
```

**Reference:**

<small>[[1]](#a1) <span id="f1"></span> Yan, Z. et al. Genome-wide co-localization of RNA-DNA interactions and fusion RNA pairs. bioRxiv 472019 (2018). doi:10.1101/472019</small>