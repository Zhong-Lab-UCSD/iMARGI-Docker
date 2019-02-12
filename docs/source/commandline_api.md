# Command-line API

We created several script tools. Here we show the usage and source code of all these tools.

## imargi_wrapper.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_wrapper.sh)

```
    Usage: $PROGNAME [-r <ref_name>] [-N <base_name>] [-g <ref_fasta>]
                     [-c <chromSize_file>] [-i <bwa_index>] [-R <restrict_sites>]
                     [-G <max_inter_align_gap>] [-O <offset_restriction_site>] [-M <max_ligation_size>]
                     [-t <threads>] [-1 <fastq.gz_R1>] [-2 <fastq.gz_R2>] [-o <output_dir>] 
    
    Dependency: seqtk, samtools, bwa, pairtools, pbgzip

    This is an all-in-one wrapper script of iMARGI Pipeline.

    -r : Reference assembly name, such as "hg38"
    -N : Base name for ouput result. Such as -N HEK_iMARGI, then the output file will be renamed using the base name,
         such as final_HEK_iMARGI.paris.gz, HEK_iMARGI.bam.
    -g : Reference genome fasta file. It will be used to generate ref files of '-c', '-i' and '-R'. If any of them are
         not supplied. It will be generated from the reference fasta file.
    -c : Chromosome size file.
    -i : bwa index
    -R : DNA restriction enzyme digestion fragments bed file.
    -G : Max inter align gap for pairtools parsing. Default 20. It will allow R1 5' end clipping.
    -O : Max offset bases for filtering pairs based on R2 5' end positions to restriction sites. Default 3.
    -M : Max size of ligation fragment for sequencing. It's used for filtering unligated DNA sequence.
    -t : Max CPU threads for parallelized processing, at least 4. (Default 8)
    -1 : R1 fastq.gz file, if there are multiple files, just separated with space, such as -1 lane1_R1.fq lane2_R1.fq
    -2 : R2 fastq.gz file, if there are multiple files, just separated with space, such as -1 lane1_R2.fq lane2_R2.fq
    -o : Output directoy
    -h : Show usage help
```

## imargi_clean.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_clean.sh)

``` 
    Usage: imargi_clean.sh [-1 <fastq.gz_R1>] [-2 <fastq.gz_R2>] [-o <output_dir>] [-f <filter_CT>] [-d <drop>] 
                [-t <threads>] [-b <block_size>]
    
    Dependencies: seqtk, gzip, zcat, awk, parallel

    This script will clean the paired reads (R1 and R2) of iMARGI sequencing Fastq data. According to the iMARGI design,
    RNA end reads (R1) start with 2 random based, and DNA end reads (R2) of successful ligation fragments start
    with "CT". We need to remove the first 2 bases of R1 for better mapping. For R2, We can strictly clean the data by
    filtering out those R2 reads not starting with "CT" in this step by setting "-f CT". Alternatively, 
    you can skip the "CT" filtering without "-f" parameter. You can apply the "CT" filtering in interaction pairs
    filtering step.
    If you choose to do "CT" filtering, the script also fixes the paired reads in R1. If "-d" was set as "true",
    it will drop all the non "CT" started R2 reads and paired R1 reads, which outputs two fastq files with prefix
    "clean_". If "-d" was "false", the filtered read pairs would also be outputed in a pair of fastq files with prefix
    "drop_". "-d" only works when "-f" is set, and the default setting of "-d" is "false".
    The input fastq files must be gzip files, i.e., fastq.gz or fq.gz. The output files are also gzipped files fastq.gz.

    -1 : R1 fastq.gz file, if there are multiple files, just separated with space, such as -1 lane1_R1.fq lane2_R1.fq
    -2 : R2 fastq.gz file, if there are multiple files, just separated with space, such as -1 lane1_R2.fq lane2_R2.fq
    -o : Output directoy
    -f : Filtering sequence by 5' start of R2. If not set, no filtering applied. "CT" filtering can be set as "-f CT"
    -d : Flag of dropping, working with "-f". Default is false, i.e., output drop_*fastq.gz files of dropped read pairs.
    -t : Max CPU threads for parallelized processing, at least 4. (Default 8)
    -b : Fastq data block size (number of reads) for each thread. Default 2000000.
    -h : Show usage help
```

## imargi_rsfrags.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_rsfrags.sh)

```
    Usage: imargi_rsfrags.sh [-r <ref_fasta>] [-c <chromSize_file>] [-e <enzyme_name>] [-C <cut_position>] [-o <output_dir>] 
                    [-g <max_inter_align_gap>] [-O offset_restriction_site] [-d <drop>] [-D <intermediate_dir>] 
                    [-s <stats_flag>] [-t <threads>] 
    
    Dependency: cooler

    This script use cooler digest to generate the restriction Enzyme digested fragments bed file for iMARGI

    -r : Reference genome fasta file
    -c : Chromosome size file.
    -e : Enzyme name, we use AluI in iMARGI.
    -C : Cut position right offset of the Enzyme sequence. In iMARGI, the enzyme restriction sequence is AGCT and it
         cuts between G and C, so the right offset is 2, then -C 2.
    -o : Output file name with directory. It's a .gz file (tabix indexed bgzip file).
    -h : Show usage help
```

## imargi_parse.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_parse.sh)

```
    Usage: imargi_parse.sh [-r <ref_name>] [-c <chromSize_file>] [-R <restrict_sites>] [-b <bam_file>] [-o <output_dir>] 
                    [-Q <min_mapq>] [-G <max_inter_align_gap>] [-O <offset_restriction_site>] [-M <max_ligation_size>]
                    [-d <drop>] [-D <intermediate_dir>] [-t <threads>] 
    
    Dependency: pairtools pbgzip

    This script will use pairtools to parse the BAM alignments to interaction read pairs in .pairs format, and apply
    de-duplication and filtering.

    -r : Reference assembly name, such as "hg38"
    -c : Chromosome size file.
    -R : DNA restriction enzyme digestion sites bed file.
    -b : BAM file generated by "bwa mem -SP5M" mapping of iMARGI data.
    -o : Output directoy
    -Q : Min MAPQ value, default 1.
    -G : Max inter align gap for pairtools parsing. Default 20. It will allow R1 5' end clipping.
    -O : Max mis-offset bases for filtering pairs based on R2 5' end positions to restriction sites. Default 0.
    -M : Max size of ligation fragment for sequencing. It's used for filtering unligated DNA sequence.
    -d : Flag of dropping. Default is false, i.e., output all the intermediate results.
    -D : Directory for intermediate results. Works when -d false. Default is a sub-folder "intermediate_results" 
         in output directory.
    -t : Max CPU threads for parallelized processing, at least 4. (Default 8)
    -h : Show usage help
```

## imargi_stats.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_stats.sh)

```
    Usage: imargi_stats.sh [-D <distance_type>] [-d <distance_threshold>] [-F <deal_with_filter>]
                [-i <input_file>] [-o <output_file>]

    Dependency: gzip, awk

    This script can be used to filter out short-range intra-chromosomal interactions with a threshold genomic distance.

    -D : Distance type. The default genomic position in .pairs file is the 5' end position, so the default distance
         type is the distance between 5' end of R1 and R2, which meand '-D 5end'. The alternative choices are 'outer'
         and 'inner'. 'outer' means the outer-side distance between R1 and R2 read pairs. 'inner' means the inner-side
         distance between R1 and R2 read pairs. Default is '5end'.
    -d : sets the genomic distance threshold of defining proximal and distal interactions (intra-chromosomal). It will
         report the number of proximal and distal interactions. Besides, you can set a list of values separated with
         comma ',', such as '-d 1000,2000,10000,20000,100000,1000000', then the report will include the statistics
         number with different thresholds (space is not allowed).
    -i : Input .pairs.gz file.
    -o : Output .pairs.gz file.
    -h : Show usage help
```

## imargi_distfilter.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_distfilter.sh)

```
    Usage: imargi_distfilter.sh [-D <distance_type>] [-d <distance_threshold>] [-F <deal_with_filter>]
                [-i <input_file>] [-o <output_file>]

    Dependency: gzip, awk

    This script can be used to filter out short-range intra-chromosomal interactions with a threshold genomic distance.

    -D : Distance type. The default genomic position in .pairs file is the 5' end position, so the default distance
         type is the distance between 5' end of R1 and R2, which meand '-D 5end'. The alternative choices are 'outer'
         and 'inner'. 'outer' means the outer-side distance between R1 and R2 read pairs. 'inner' means the inner-side
         distance between R1 and R2 read pairs. Default is '5end'.
    -d : The distance threshold for filtering. Default is 200000 (distance <200000 will be filtered out).
    -F : How to deal with the interactions need to be filtered out? '-F' accepts 'drop' and 'output'. 'drop' means
         drop those interactions. 'output' means output an new file of all the filtered out interactions with prefix
         'filterOut_'. Default is 'drop'.
    -i : Input .pairs.gz file.
    -o : Output .pairs.gz file.
    -h : Show usage help
```

## imargi_convert.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_convert.sh)

```
    Usage: imargi_convert.sh [-f <file_format>] [-k <keep_cols>] [-b <bin_size>] [-i <input_file>] [-o <output_file>]

    Dependency: gzip, awk, cool

    This script can convert .pairs format to BEDPE, .cool, and GIVE interaction format.

    -f : The target format, only accept 'cool', 'bedpe' and 'give'. For 'cool', it will generate
         a ".cool" file with defined resolution of -b option and a multi-resolution ".mcool" file
         based on the ".cool" file.
    -k : Keep extra information column in BEDPE. Columns ids in .pairs file you want to keep.
         For example, 'cigar1,cigar2'. Default value is "", i.e., drop all extra cols.
    -b : bin size for cool format. Default is 5000.
    -i : Input file.
    -o : Output file.
    -h : Show usage help
```

## imargi_annotate.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_annotate.sh)

```
    Usage: imargi_annotate.sh [-f <file_format>] [-A <annotation_type>] [-a <annotation_file>]
                [-C <end_for_annotate>] [-c <add_col_names>] [-m <min_overlap>]
                [-i <input_file>] [-o <output_file>]

    Dependency: gzip, awk, cool, BEDOPS

    This script can annotate both RNA and DNA ends with gene annotations in GTF/GFF format or any other genomic
    features in a simple BED file (each line is a named genomic feature).

    -f : The input file format, only accept 'pairs', or 'bedpe'.
    -A : Only accept '-A gtf' for annotating gene with GTF file or '-A bed' for annotating any other genomic features
         in a simple BED format file. Default is '-A gtf'.
    -a : The annotation information file.
    -C : Determine Which end to be annotated. It accepts '-C RNA', '-C DNA' or '-C both'. Default is '-C both',
         then both RNA and DNA ends will be annotated.
    -c : The names of annotation column in the output file. If '-C both', then you need set two column names separated
         with comma ',', i.e., '-c gene1,gene2'.
    -m:  Minimum bases of overlapping for annotation. '-m 1' means at least one base overlapping. Default '-m 1'.
         '-m 0' means it must be inside of the genomic feature.
    -i : Input file.
    -o : Output file.
    -h : Show usage help
```