#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
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
EOF
    exit 1
}

while getopts :1:2:N:o:t:h opt; do
    case $opt in
        1) R1=("${OPTARG}")
            until [[ $(eval "echo \${$OPTIND}") =~ ^-.* ]] || [ -z $(eval "echo \${$OPTIND}") ]; do
                R1+=($(eval "echo \${$OPTIND}"))
                OPTIND=$((OPTIND + 1))
            done;;
        2) R2=("${OPTARG}")
            until [[ $(eval "echo \${$OPTIND}") =~ ^-.* ]] || [ -z $(eval "echo \${$OPTIND}") ]; do
                R2+=($(eval "echo \${$OPTIND}"))
                OPTIND=$((OPTIND + 1))
            done;;
        N) base_name=${OPTARG};;
        o) output_dir=${OPTARG};;
        t) threads=${OPTARG};;
        h) usage;;
    esac
done

for i in ${R1[@]}; do
    [ ! -f "$i" ] && echo "Error!! Fastq R1 not exist: "$i && usage
done

for i in ${R2[@]}; do
    [ ! -f "$R2" ] && echo "Error!! Fastq R2 not exist: "$i && usage
done

R1_str=''
for i in ${R1[@]};do
    R1_str=$R1_str' '$i
done

R2_str=''
for i in ${R2[@]};do
    R2_str=$R2_str' '$i
done

[ ! -d "$output_dir" ] && echo "Error!! Output directory not exist: "$output_dir && usage

[  -z "$threads" ] && echo "Use default thread number 8'." && threads=2
if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -t" && usage 
fi

clean_R1=$output_dir"/clean_"$base_name"_R1.fastq.gz"
clean_R2=$output_dir"/clean_"$base_name"_R2.fastq.gz"

[ -f "$clean_R1" ] && echo "Error!! Output clean fastq file exists: "$clean_R1 && usage
[ -f "$clean_R2" ] && echo "Error!! Output clean fastq file exists: "$clean_R2 && usage

t_pbgzip=$(( $threads - 1 ))

echo "Start cleaning:"

echo "    Remove first 2 bases of R1 reads and merge: " $R1_str
zcat $R1_str | seqtk trimfq -b 2 - | pbgzip -n $t_pbgzip -t 0 -c  > $clean_R1

echo "    Copy and merge (if needed) R2 reads: "$R2_str
cat $R2_str > $clean_R2

echo "Finished: cleaned fastq files are:"
echo "    $output_dir/$clean_R1 and $output_dir/$clean_R2"