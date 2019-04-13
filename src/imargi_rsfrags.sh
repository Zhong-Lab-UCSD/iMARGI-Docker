#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
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
EOF
    exit 1
}

while getopts :r:c:e:C:o:h opt; do
    case $opt in
        r) ref_fasta=${OPTARG};;
        c) chromsize=${OPTARG};;
        e) enzyme=${OPTARG};;
        C) cut_offset=${OPTARG};;
        o) output_file=${OPTARG};;
        h) usage;;
    esac
done

[ ! -f "$ref_fasta" ] && echo "Error!! Please provide reference genome fasta with -r" && usage
[ ! -f "$chromsize" ] && echo "Error!! Chomosome size file not exist: "$chromsize && usage
[ -z "$enzyme" ] && echo "Enzyme name is required, use AluI for iMARGI." && usage
[ -z "$cut_offset" ] && echo "Cut position is required, use AluI for iMARGI." && usage
[ -f "$output_file" ] && echo "Error!! Outputfile already exists: "$output_file && usage

cooler digest $chromsize $ref_fasta $enzyme |\
    awk -v offset="$cut_offset" \
       'BEGIN{OFS="\t"}
        NR==FNR{chrsize[$1]=$2}
        NR!=FNR{if($2==0){start=0}else{start=$2-offset}; if($3==chrsize[$1]){end=$3}else{end=$3-offset};
        print $1, start, end}' \
       $chromsize - | bgzip > $output_file

tabix -p bed $output_file
