#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
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
    -1 : R1 fastq.gz file, if there are multiple files, just separated with space or use wildcard,
         such as '-1 lane1_R1.fq.gz lane2_R1.fq.gz', or '-1  lane*_R1.fq.gz'.
    -2 : R2 fastq.gz file, if there are multiple files, just separated with space or use wildcard,
         such as '-2 lane1_R2.fq.gz lane2_R2.fq.gz', or '-2  lane*_R2.fq.gz'.
    -o : Output directoy
    -h : Show usage help
EOF
    exit 1
}

# imargi_wrapper.sh -r hg38 -N HEK_iMARGI -c ../imargi_example/ref/hg38.chrom.sizes  -i ../imargi_example/ref/bwa_index/bwa_index_GRCh38 -R ../ref/iMARGI_AluI_rsites.bed.gz -t 16 -1 test_R1.fastq.gz  -2 test_R2.fastq.gz  -o ./wrapper_output

while getopts :r:N:g:c:i:R:G:O:M:t:1:2:o:h opt; do
    case $opt in
        r) ref_name=${OPTARG};;
        N) base_name=${OPTARG};;
        g) ref_fa=${OPTARG};;
        c) chromsize=${OPTARG};;
        i) bwa_index=${OPTARG};;
        R) rsites=${OPTARG};;
        G) gap=${OPTARG};;
        O) offset=${OPTARG};;
        M) max_ligation_size=${OPTARG};;
        t) threads=${OPTARG};;
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
        o) output_dir=${OPTARG};;
        h) usage;;
    esac
done

[ -z "$ref_name" ] && echo "Error!! Please provide reference genome name with -r" && usage
[ -z "$base_name" ] && echo "Error!! Please provide base name for result files with -N" && usage

[ ! -d "$output_dir" ] && echo "Error!! Output directory not exist: "$output_dir && usage

[  -z "$gap" ] && echo "Use default max inter align gap for pairtools parsing." && gap=20
if ! [[ "$gap" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -G" && usage 
fi

[  -z "$offset" ] && echo "Use default offset 3'." && offset=3
if ! [[ "$offset" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -O" && usage
fi

[  -z "$max_ligation_size" ] && echo "Use default max ligation size 1000'." && max_ligation_size=1000
if ! [[ "$max_ligation_size" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -M" && usage 
fi

[  -z "$threads" ] && echo "Use default thread number 8'." && threads=8
if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -t" && usage 
fi

[ ! -d "$output_dir/clean_fastq" ] && mkdir $output_dir/clean_fastq
[ ! -d "$output_dir/bwa_output" ] && mkdir $output_dir/bwa_output
[ ! -d "$output_dir/parse_temp" ] && mkdir $output_dir/parse_temp

echo "================ iMARGI Pipeline BEGIN"

if [ ! -z "$ref_fa" ]; then
    [ ! -f "$ref_fa" ] && echo "Error!! Reference genome fasta file not exist: "$ref_fa && usage
else
    [ ! -f "$chromsize" ] && echo "Error!! Chomosome size file not exist: "$chromsize && usage
    [ -z "$bwa_index" ] && echo "Error!! Please provide correct bwa index with -i" && usage
    [ ! -f "$rsites" ] && echo "Error!! Resitriction sites file not exist: "$rsites && usage
fi

if [ ! -z "$ref_fa" ]; then
    if [ -z "$chromsize" ]; then
        date
        echo ">>>>>>>>>>>>>>>> Ref: no chromsize argument. Start generating chromsize file from "$ref_fa" ..."
        [ ! -f "$ref_fa" ] && echo "Error!! Reference genome fasta file not exist: "$ref_fa && usage
        ref_dir=$(dirname "$ref_fa")
        cd $ref_dir && samtools faidx $(basename "$ref_fa") && cd -
        awk 'BEGIN{OFS="\t"}{print $1,$2}' "$ref_fa".fai > $ref_dir/chromsize.txt
        chromsize=$ref_dir"/chromsize.txt"
    else
        [ ! -f "$chromsize" ] && echo "Error!! Chomosome size file not exist: "$chromsize && usage
        date
        echo ">>>>>>>>>>>>>>>> Ref: chromsize argument exist. Directly use the chromsize file "$chromsize
    fi

    if [ -z "$bwa_index" ]; then
        date
        echo ">>>>>>>>>>>>>>>> Ref: no bwa index argument. Start generating bwa index from "$ref_fa" ..."
        [ ! -f "$ref_fa" ] && echo "Error!! Reference genome fasta file not exist: "$ref_fa && usage
        ref_dir=$(dirname "$ref_fa")
        mkdir $ref_dir/bwa_index
        bwa index -p $ref_dir/bwa_index/bwa_index_$ref_name $ref_fa
        bwa_index=$ref_dir"/bwa_index/bwa_index_"$ref_name
    else
        date
        echo ">>>>>>>>>>>>>>>> Ref: bwa index argument exist. Directly use the bwa index "$chromsize      
    fi

    if [ -z "$rsites" ]; then
        date
        echo ">>>>>>>>>>>>>>>> Ref: no restriction fragment argument. Start generating from "$ref_fa" ..."
        [ ! -f "$ref_fa" ] && echo "Error!! Reference genome fasta file not exist: "$ref_fa && usage
        ref_dir=$(dirname "$ref_fa")
        imargi_rsfrags.sh -r $ref_fa -c $chromsize -e AluI -C 2 -o $ref_dir/AluI_frags.bed.gz
        rsites=$ref_dir"/AluI_frags.bed.gz"
    else
        [ ! -f "$rsites" ] && echo "Error!! Ref: Restriction fragment file not exist: "$rsites && usage
        date
        echo ">>>>>>>>>>>>>>>> Ref: Restriction fragment file argument exist. Directly use "$rsites
    fi
fi

date
echo ">>>>>>>>>>>>>>>> Start cleaning fastq ..."
R1_str=''
for i in ${R1[@]};do
    R1_str=$R1_str' '$i
done

R2_str=''
for i in ${R2[@]};do
    R2_str=$R2_str' '$i
done
echo ">>>>>>>>>>>>>>>> Start bwa mem mapping ..."
# imargi_clean.sh -1 $R1_str -2 $R2_str -N $base_name -o $output_dir/clean_fastq -t $threads

# date
# echo ">>>>>>>>>>>>>>>> Start bwa mem mapping ..."

# bwa mem -t $threads -SP5M $bwa_index \
#     $output_dir/clean_fastq/clean_${base_name}_R1.fastq.gz \
#     $output_dir/clean_fastq/clean_${base_name}_R2.fastq.gz | \
#     samtools view -@ $threads -Shb - > $output_dir/bwa_output/$base_name.bam

date
echo ">>>>>>>>>>>>>>>> Start parsing valid RNA-DNA interactions ..."
imargi_parse.sh -r $ref_name -c $chromsize -R $rsites -b $output_dir/bwa_output/$base_name.bam \
    -o $output_dir -G 20 -O 3 -M 1000 -d false -D $output_dir/parse_temp -t $threads

date
echo "================ iMARGI Pipeline END"
