#!/usr/bin/env bash
PROGNAME=$0

usage() {
    cat << EOF >&2
    Usage: $PROGNAME [-1 <fastq.gz_R1>] [-2 <fastq.gz_R2>] [-o <output_dir>] [-f <filter_CT>] [-d <drop>] 
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
EOF
    exit 1
}

while getopts :1:2:o:f:d:t:b:h opt; do
    case $opt in
        1) R1=${OPTARG};;
        2) R2=${OPTARG};;
        o) output_dir=${OPTARG};;
        f) filter_seq=${OPTARG};;
        d) dflag=${OPTARG};;
        t) threads=${OPTARG};;
        b) block=${OPTARG};;
        h) usage;;
    esac
done

[ ! -f "$R1" ] && echo "Error!! Fastq R1 not exist: "$R1 && usage
[ ! -f "$R2" ] && echo "Error!! Fastq R2 not exist: "$R2 && usage
[ ! -d "$output_dir" ] && echo "Error!! Output directory not exist: "$output_dir && usage

if [  -z "$filter_seq" ]; then
    echo "No -f parameter. No filtering applied."
else
    if ! [[ "$filter_seq" =~ ^[ATGCatgc]+$ ]]; then
        echo "Error!! Only string of ATGC is acceptable for -f, such as CT" && usage
    else
        [  -z "$dflag" ] && echo "Use default setting '-d false'." && dflag="false"
        if [[ "$dflag" != "false" ]] && [[ "$dflag" != "true" ]]; then
            echo "Error!! Only true or false is acceptable for -d." && usage
        fi
    fi
fi

[  -z "$threads" ] && echo "Use default thread number 8'." && threads=2
if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -t" && usage 
fi

if (( $threads < 4)); then
    echo "Error!! Threads must be >= 4" && usage 
fi

[  -z "$block" ] && echo "Use default block size of reads for each thread 2000000'." && block=2000000
if ! [[ "$block" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -b" && usage 
fi

clean_R1=$output_dir"/clean_"$(basename $R1)
clean_R2=$output_dir"/clean_"$(basename $R2)
drop_R1=$output_dir"/drop_"$(basename $R1)
drop_R2=$output_dir"/drop_"$(basename $R2)
[ -f "$clean_R1" ] && echo "Error!! Output clean fastq file exists: "$clean_R1 && usage
[ -f "$clean_R2" ] && echo "Error!! Output clean fastq file exists: "$clean_R2 && usage
[ -f "$drop_R1" ] && echo "Error!! Output drop fastq file exists: "$drop_R1 && usage
[ -f "$drop_R2" ] && echo "Error!! Output drop fastq file exists: "$drop_R2 && usage

tmpdir=$output_dir"/tmp_cleanfq_"$RANDOM""$RANDOM

if [  -z "$filter_seq" ]; then
    t_pbgzip=$(( $threads - 1 ))
    echo "Start: only remove first 2 bases of reads $R1 using $threads CPU threads, 
          $R2 will be copied to output dir as $clean_R2 ..."
    seqtk trimfq -b 2 $R1 | pbgzip -n $t_pbgzip -t 0 -c  > $clean_R1
    cp $R2 $clean_R2
else
    t_pipe=$(( $threads / 2 ))
    t_pbgzip=$(( $threads - $t_pipe ))
    mkdir $tmpdir
    echo "Start: Remove first 2 bases of reads $R1, 
          filtering $R2 by $filter_seq 
          using $threads CPU threads with $block reads per single-thread per run ..."
    seqtk mergepe $R1 $R2 | parallel -j $t_pipe --pipe -L8 -N$block \
        "awk -v dflag=\"$dflag\" -v filter_seq=\"$filter_seq\" 'BEGIN{keep=0; 
                    srand(systime()\"\"{%});
                    randname=int(1000000*rand());
                    OFS=\"\n\";
                    IGNORECASE=1;
                }{
                    idx=NR%8;
                    if(idx==6){
                        if(\$0 ~ \"^\"filter_seq){
                            keep=1;
                        }else{
                            keep=0;
                        }
                    };
                    if(idx==2 || idx==4){
                        tmp=substr(\$0, 3);
                        seqinfo[idx]=tmp;
                    }else{
                        seqinfo[idx]=\$0;
                    };
                    if(idx==0){
                        if(keep==1){
                            print seqinfo[1],seqinfo[2],seqinfo[3],seqinfo[4] | \
                                \"pbgzip -n $t_pbgzip -t 0 -c > $tmpdir/clean_\"randname\"_R1.fastq.gz\";
                            print seqinfo[5],seqinfo[6],seqinfo[7],seqinfo[0] | \
                                \"pbgzip -n $t_pbgzip -t 0 -c > $tmpdir/clean_\"randname\"_R2.fastq.gz\";
                        }else{
                            if(dflag==\"false\"){
                                print seqinfo[1],seqinfo[2],seqinfo[3],seqinfo[4] | \
                                    \"pbgzip -n $t_pbgzip -t 0 -c > $tmpdir/drop_\"randname\"_R1.fastq.gz\";
                                print seqinfo[5],seqinfo[6],seqinfo[7],seqinfo[0] | \
                                    \"pbgzip -n $t_pbgzip -t 0 -c > $tmpdir/drop_\"randname\"_R2.fastq.gz\";
                            }
                        }            
                    };
                }'"
    if [ "$dflag" == "false" ]; then   
        cat $tmpdir/drop_*_R1.fastq.gz >$drop_R1
        cat $tmpdir/drop_*_R2.fastq.gz >$drop_R2
    fi
    cat $tmpdir/clean_*_R1.fastq.gz >$clean_R1
    cat $tmpdir/clean_*_R2.fastq.gz >$clean_R2
    rm -r $tmpdir
fi
