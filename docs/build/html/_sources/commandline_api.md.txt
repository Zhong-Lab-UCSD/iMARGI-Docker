# Command-line API

We created several script tools. Here we show the usage and source code of all these tools.

- [Command-line API](#command-line-api)
  - [imargi_wrapper.sh](#imargi_wrappersh)
  - [imargi_clean.sh](#imargi_cleansh)
  - [imargi_rsfrags.sh](#imargi_rsfragssh)
  - [imargi_parse.sh](#imargi_parsesh)
  - [imargi_stats.sh](#imargi_statssh)
  - [imargi_distfilter.sh](#imargi_distfiltersh)
  - [imargi_convert.sh](#imargi_convertsh)
  - [imargi_annotate.sh](#imargi_annotatesh)

## imargi_wrapper.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_wrapper.sh)

``` bash
    Usage: $PROGNAME [-r <ref_name>] [-N <base_name>] [-g <ref_fasta>]
                     [-c <chromSize_file>] [-i <bwa_index>] [-R <restrict_sites>] 
                     [-Q <min_mapq>] [-G <max_inter_align_gap>]
                     [-O <offset_restriction_site>] [-M <max_ligation_size>] [-t <threads>]
                     [-1 <fastq.gz_R1>] [-2 <fastq.gz_R2>]
                     [-o <output_dir>]
    
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
    -Q : Min MAPQ value for parsing as unique mapping. Default 1.
    -G : Max inter align gap for pairtools parsing. Default 20. It will allow R1 5' end clipping.
    -O : Max offset bases for filtering pairs based on R2 5' end positions to restriction sites. Default 3.
    -M : Max size of ligation fragment for sequencing. It's used for filtering unligated DNA sequence.
    -t : Max CPU threads for parallelized processing, at least 4. (Default 8)
    -1 : R1 fastq.gz file, if there are multiple files, just separated with space or use wildcard,
         such as '-1 lane1_R1.fq.gz lane2_R1.fq.gz', or '-1  lane*_R1.fq.gz'.
    -2 : R2 fastq.gz file, if there are multiple files, just separated with space or use wildcard,
         such as '-2 lane1_R2.fq.gz lane2_R2.fq.gz', or '-2  lane*_R2.fq.gz'.
    -o : Output directoy
    -h : Show usage help
```

## imargi_clean.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_clean.sh)

``` bash
    Usage: $PROGNAME [-1 <fastq.gz_R1>] [-2 <fastq.gz_R2>] [-N <base_name>] [-o <output_dir>] [-t <threads>]
    
    Dependencies: seqtk, gzip, zcat, awk, parallel

    This script will clean the paired reads (R1 and R2) of iMARGI sequencing Fastq data. According to the iMARGI design,
    RNA end reads (R1) start with 2 random based. We need to remove the first 2 bases of R1 for better mapping. 
    If you provided multiple input files (different lanes) in '-1' and '-2' with ',' separator or contains wildcard,
    then the output will merge multi-lanes fastq files to one clean fastq file.

    The input fastq files must be gzip files, i.e., fastq.gz or fq.gz. The output files are also gzipped files fastq.gz.

    -1 : R1 fastq.gz file, if there are multiple files, just separated with space or use wildcard,
         such as '-1 lane1_R1.fq.gz lane2_R1.fq.gz', or '-1  lane*_R1.fq.gz'.
    -2 : R2 fastq.gz file, if there are multiple files, just separated with space or use wildcard,
         such as '-2 lane1_R2.fq.gz lane2_R2.fq.gz', or '-2  lane*_R2.fq.gz'.
    -N : Base name for ouput result. Such as -N HEK_iMARGI, then output cleaned and merged fastq.gz file will be
         renamed using the base name.
    -o : Output directoy
    -t : Max CPU threads for parallelized processing, at least 4. (Default 8)
    -h : Show usage help
```

## imargi_rsfrags.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_rsfrags.sh)

``` bash
    Usage: $PROGNAME [-r <ref_fasta>] [-c <chromSize_file>] [-e <enzyme_name>] [-C <cut_position>] [-o <output_dir>] 

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

``` bash
    Usage: $PROGNAME [-r <ref_name>] [-c <chromSize_file>] [-R <restrict_sites>] [-b <bam_file>] [-o <output_dir>] 
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
    -O : Max mis-offset bases for filtering pairs based on R2 5' end positions to restriction sites. Default 3.
    -M : Max size of ligation fragment for sequencing. It's used for filtering unligated DNA sequence. Default 1000.
    -d : Flag of dropping. Default is false, i.e., output all the intermediate results.
    -D : Directory for intermediate results. Works when -d false. Default is a sub-folder "intermediate_results" 
         in output directory.
    -t : Max CPU threads for parallelized processing, at least 4. (Default 8)
    -h : Show usage help
```

## imargi_stats.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_stats.sh)

``` bash
    Usage: $PROGNAME [-D <distance_type>] [-d <distance_threshold>] [-i <input_file>] [-o <output_file>]

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
    -o : Output stats text file.
    -h : Show usage help
```

## imargi_distfilter.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_distfilter.sh)

``` bash
    Usage: $PROGNAME [-D <distance_type>] [-d <distance_threshold>] [-F <deal_with_filter>] [-i <input_file>] 
                     [-o <output_file>]

    Dependency: gzip, awk
    This script can be used to filter out short-range intra-chromosomal interactions with a threshold genomic distance.

    -D : Distance type. The default genomic position in .pairs file is the 5' end position, so the default distance
         type is the distance between 5' end of R1 and R2, which meand '-D 5end'. The alternative choices are 'outer'
         and 'inner'. 'outer' means the outer-side distance between R1 and R2 read pairs. 'inner' means the inner-side
         distance between R1 and R2 read pairs. Default is '5end'.
    -d : The distance threshold for filtering. Default is 200000 (distance <200000 will be filtered out).
    -F : How to deal with the interactions need to be filtered out? '-F' accepts 'drop' and 'output'. 'drop' means
         drop those interactions. 'output' means output an new file of all the filtered out interactions with prefix
         'filterOut_'. Default is 'output'.
    -i : Input .pairs.gz file.
    -o : Output .pairs.gz file.
    -h : Show usage help
```

## imargi_convert.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_convert.sh)

``` bash
    Usage: $PROGNAME [-f <file_format>] [-k <keep_cols>] 
                     [-b <bin_size>] [-r <resolution>] [-T <transpose>] 
                     [-i <input_file>] [-o <output_file>] 

    Dependency: gzip, awk, cool
    This script can convert .pairs format to BEDPE, .cool, and GIVE interaction format.
    -f : The target format, only accept 'cool', 'bedpe' and 'give'. For 'cool', it will generate a ".cool" file
         with defined resolution of -b option and a -r defined multi-resolution ".mcool" file based on the ".cool" file.
         For 'bedpe', the output will be pbgzip compressed file. So keep in mind to name the output_file '-o' with
         '.gz' extesion. For 'give', the output is a normal text file.
    -k : (Only for BEDPE) Keep extra information column in BEDPE. Columns ids in .pairs file you want to keep.
         For example, 'cigar1,cigar2'. Default value is "", i.e., drop all extra cols.
    -b : (Only for cool/mcool) bin size for cool format. Default is 1000.
    -r : (Only for cool/mcool) resolution for cool/mcool format. Integers separated by comma. The values of resolution
         must be integer multiples of the bin size defined by -b option.
         Default is 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000
    -T : (Only for cool/mcool) mcool file can be visualized by HiGlass. Currently the heatmap orientation cannot be
         set in HiGlass control panel. So if you want to transpose the interaction map in HiGlass, you need to generate
         a transposed mcool file. The default value is 'fasle', i.e., no transpose, the RNA-DNA interactions will be
         mapped to a X-Y system as a DNA x RNA contact matrix. If set '-T true', then the generated cool/mcool file is 
         transposed, which is RNA x DNA contact matrix.
    -i : Input file.
    -o : Output file.
         BEDPE output is gzip compressed file, so it's better to have a .gz file extension.
         cool output are two files, .cool and .mcool. The -o option assigns the name of .cool file, it must use .cool as
         extension. The .mcool file will be generated based on the .cool file with .mcool extension.
    -h : Show usage help
```

## imargi_annotate.sh

[*Source Code*](https://github.com/Zhong-Lab-UCSD/iMARGI-Docker/blob/master/src/imargi_annotate.sh)

``` bash
    Usage: imargi_annotate.sh [-A <annotation_format>] [-a <annotation_file>] 
                [-l <annotation_level>] [-f <feature_attribute>]
                [-C <end_for_annotate>] [-c <add_col_names>] [-s <strand_specific>]
                [-m <min_overlap>] [-G <cigar>]
                [-t <threads>] [-i <input_file>] [-o <output_file>]

    Dependency: gzip, pairtools, lz4 pbgzip

    This script can annotate both RNA and DNA ends with gene annotations in GTF/GFF format or any other genomic
    features in a simple BED file (each line is a named genomic feature). Multiple overlapped annotation features are
    separated by ','.

    -A : Only accept '-A gtf' for annotating gene with GTF file or '-A bed' for annotating any other genomic features
         in a simple BED format file, whose 4th column is the feature name. Default is '-A gtf'.
    -a : The annotation information file.
    -l : Only work with '-A gtf'. Gene annotation level, 'exon' or 'gene'. Default is '-l gene'.
    -f : Only work with '-A gtf'. There are many gene feature attributes can be added, such as gene_id, gene_name,
         or gene_type. Use ',' to seperate the features. In the output results, these features are seperated with
         '|'. Default is '-f gene_id,gene_name,gene_type'. Blank space is not allowed!
    -C : Determine Which end to be annotated. It accepts '-C RNA', '-C DNA' or '-C both'. Default is '-C both',
         then both RNA and DNA ends will be annotated. When use '-C RNA' or '-C DNA', then '-c', '-s', '-m' and '-G'
         only need one value, if you provided two, then only the first one will be used. If use '-C both', then you must
         provide two values for '-c', '-s', '-m' and '-G'.
    -c : The names of annotation column in the output file. If '-C both', then you need set two column names separated
         with comma ',', i.e., '-c gene1,gene2'. Blank space is not allowed!
    -s : Strand specific annotation option, 's', 'r' and 'n'. 's' means strand specific, 'r' means
         reverse-strand-specific, and 'n' means none strand specific, ignoring it. As RNA end of iMARGI is
         reverse-strand-specific, and DNA end is none strand specific, so the default is '-s rn'.
    -m : Minimum bases of overlapping for annotation. '-m 1' means at least one base overlapping. Default '-m 1,1'.
         '-m 0' means it must be totally inside of the genomic feature.
    -G : Use CIGAR columsn information . Default value is use standard CIGAR info columns, i.e., '-G cigar1,cigar2'.
         Blank space is not allowed! If you just want to use the 5' end genomic coordinate reported in .pairs file for
         annotation, set '-G false,false' and the '-m' will be automatically set as '0,0'.
    -t : Max CPU threads for parallelized processing, at least 4. (Default 8)
    -i : Input file.
    -o : Output file.
    -h : Show usage help
```
