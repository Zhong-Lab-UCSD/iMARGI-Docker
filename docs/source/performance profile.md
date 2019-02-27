## bwa index

```
docker run -v ~/research/MARGI_4DN/imargi_example:/imargi zhonglab/imargi bwa index \
    -p ./ref/bwa_index/bwa_index_GRCh38 ./ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
```

 Version: 0.7.17-r1194-dirty
 CMD: bwa index -p ./ref/bwa_index/bwa_index_GRCh38 ./ref/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
 Real time: 4839.265 sec; CPU: 4822.200 sec

## clean

``` bash
time bash tools/imargi_clean_fastq.sh -1  raw_fastq/SRR8206679_1.fastq.gz -2 raw_fastq/SRR8206679_2.fastq.gz \
    -o output -t 16 -b 2000000
```

6101.94user 2459.60system 11:32.13elapsed 1236%CPU (0avgtext+0avgdata 32000maxresident)k
75975280inputs+85582640outputs (1major+417576335minor)pagefaults 0swaps

## mapping

``` bash
time bwa mem -t 16 -SP5M ./ref/bwa_index/bwa_index_GRCh38 \
    ./clean_fastq/clean_SRR8206679_1.fastq.gz ./clean_fastq/clean_SRR8206679_2.fastq.gz | \
    samtools view -@ 15 -Shb - >./output/HEK_iMARGI.bam
```

229676.05user 10590.82system 4:12:19elapsed 1586%CPU (0avgtext+0avgdata 9430528maxresident)k
96178192inputs+0outputs (0major+4909367983minor)pagefaults 0swaps

## parsing and filtering

Start parsing and deduplication ...
Sun Jan 20 19:28:20 PST 2019
55948.08user 5881.36system 5:27:54elapsed 314%CPU (0avgtext+0avgdata 1130424maxresident)k
122919624inputs+154081216outputs (0major+657876360minor)pagefaults 0swaps

Mon Jan 21 00:59:47 PST 2019
Start filtering ... 
Mon Jan 21 01:15:54 PST 2019
Parsing and filtering finished.
2541.66user 1929.22system 16:07.00elapsed 462%CPU (0avgtext+0avgdata 2099000maxresident)k
10272424inputs+9512520outputs (89major+68513040minor)pagefaults 0swaps
