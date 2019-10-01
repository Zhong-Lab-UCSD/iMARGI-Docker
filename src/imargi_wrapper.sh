#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2

    ---------------------------------------------------------------------------
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
    ---------------------------------------------------------------------------
EOF
    exit 1
}

parameter_error() {
    echo "!!!!!!!!!!!!!!!! Pipeline exit with parameter ERROR at $(date '+%Y-%m-%d %H:%M:%S %Z') !!!!!!!!!!!!!!!!"
    usage
}

while getopts :r:N:g:c:i:R:Q:G:O:M:t:1:2:o:h opt; do
    case $opt in
        r) ref_name=${OPTARG};;
        N) base_name=${OPTARG};;
        g) ref_fa=${OPTARG};;
        c) chromsize=${OPTARG};;
        i) bwa_index=${OPTARG};;
        R) rsites=${OPTARG};;
        Q) mapq=${OPTARG};;
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


echo "==================== iMARGI Pipeline BEGIN at $(date '+%Y-%m-%d %H:%M:%S %Z') ===================="

# check parameters
[ -z "$ref_name" ] &&  echo "Error!! Please provide reference genome name with -r" && parameter_error

[ -z "$base_name" ] && echo "Error!! Please provide base name for result files with -N" && parameter_error

[ ! -d "$output_dir" ] && echo "Error!! Output directory not exist: "$output_dir && parameter_error

[  -z "$mapq" ] && echo "Use default min mapq 1 for pairtools parsing." && mapq=1
if ! [[ "$mapq" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -Q" && parameter_error 
fi

[  -z "$gap" ] && echo "Use default max inter align gap for pairtools parsing." && gap=20
if ! [[ "$gap" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -G" && parameter_error
fi

[  -z "$offset" ] && echo "Use default offset 3'." && offset=3
if ! [[ "$offset" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -O" && parameter_error
fi

[  -z "$max_ligation_size" ] && echo "Use default max ligation size 1000'." && max_ligation_size=1000
if ! [[ "$max_ligation_size" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -M" && parameter_error
fi

[  -z "$threads" ] && echo "Use default thread number 8'." && threads=8
if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -t" && parameter_error 
fi

# check input parameter of fastq files
[ -z "$R1" ] && echo "Error!! Please provide fastq.gz files of R1 with -1" && parameter_error
[ -z "$R2" ] && echo "Error!! Please provide fastq.gz files of R2 with -2" && parameter_error

R1_str=" ${R1[*]} "
R2_str=" ${R2[*]} "
echo $R1_str
echo $R2_str

for item in ${R1[@]}; do
    if [[ $R2_str =~ " $item " ]] ; then
        echo "Error!! R1 and R2 fastq.gz files supplied by -1 and -2 are the same one:"
        echo "R1: "$R1_str
        echo "R2: "$R2_str
        parameter_error
    fi
done

for i in ${R1[@]}; do
    [ ! -f "$i" ] && echo "Error!! Fastq R1 not exist: "$i && parameter_error
done

for i in ${R2[@]}; do
    [ ! -f "$R2" ] && echo "Error!! Fastq R2 not exist: "$i && parameter_error
done

# check reference files
if [ ! -z "$ref_fa" ]; then
    [ ! -f "$ref_fa" ] && echo "Error!! Reference genome fasta file not exist: "$ref_fa && parameter_error
else
    [ ! -f "$chromsize" ] && echo "Error!! Chomosome size file not exist: "$chromsize && parameter_error 
    [ -z "$bwa_index" ] && echo "Error!! Please provide correct bwa index with -i" && parameter_error 
    [ ! -f "$rsites" ] && echo "Error!! Resitriction sites file not exist: "$rsites && parameter_error 
fi

if [ ! -z "$ref_fa" ]; then
    if [ -z "$chromsize" ]; then
        
        echo ">>>>>>>>>>>>>>>> [$(date '+%m-%d %H:%M:%S')] Ref: no chromsize parameter. Start generating chromsize file from "$ref_fa" ..."
        [ ! -f "$ref_fa" ] && echo "Error!! Reference genome fasta file not exist: "$ref_fa && parameter_error 
        ref_dir=$(dirname "$ref_fa")
        samtools faidx $ref_fa
        awk 'BEGIN{OFS="\t"}{print $1,$2}' "$ref_fa".fai > $ref_dir/chromsize.$ref_name.txt
        chromsize=$ref_dir"/chromsize.$ref_name.txt"
    else
        [ ! -f "$chromsize" ] && echo "Error!! Chomosome size file not exist: "$chromsize && parameter_error 
        
        echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d %H:%M:%S')] Ref: chromsize parameter exist. Directly use the chromsize file "$chromsize
    fi

    if [ -z "$bwa_index" ]; then
        date
        echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d %H:%M:%S')] Ref: no bwa index parameter. Start generating bwa index from "$ref_fa" ..."
        [ ! -f "$ref_fa" ] && echo "Error!! Reference genome fasta file not exist: "$ref_fa && parameter_error 
        ref_dir=$(dirname "$ref_fa")
        mkdir $ref_dir/bwa_index
        bwa index -p $ref_dir/bwa_index/bwa_index_$ref_name $ref_fa
        bwa_index=$ref_dir"/bwa_index/bwa_index_"$ref_name
    else
        [ ! -f "${bwa_index}.amb" ] && echo "Error!! BWA index files not correct: "$bwa_index && parameter_error 
        [ ! -f "${bwa_index}.ann" ] && echo "Error!! BWA index files not correct: "$bwa_index && parameter_error 
        [ ! -f "${bwa_index}.bwt" ] && echo "Error!! BWA index files not correct: "$bwa_index && parameter_error 
        [ ! -f "${bwa_index}.pac" ] && echo "Error!! BWA index files not correct: "$bwa_index && parameter_error 
        [ ! -f "${bwa_index}.sa" ] && echo "Error!! BWA index files not correct: "$bwa_index && parameter_error

        echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d %H:%M:%S')] Ref: bwa index parameter exist. Directly use the bwa index "$bwa_index      
    fi

    if [ -z "$rsites" ]; then
        echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d %H:%M:%S')] Ref: no restriction fragment parameter. Start generating from "$ref_fa" ..."
        [ ! -f "$ref_fa" ] && echo "Error!! Reference genome fasta file not exist: "$ref_fa && parameter_error 
        ref_dir=$(dirname "$ref_fa")
        imargi_rsfrags.sh -r $ref_fa -c $chromsize -e AluI -C 2 -o $ref_dir/AluI_frags.bed.gz
        rsites=$ref_dir"/AluI_frags.bed.gz"
    else
        [ ! -f "$rsites" ] && echo "Error!! Ref: Restriction fragment file not exist: "$rsites && parameter_error 

        echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d %H:%M:%S')] Ref: Restriction fragment file parameter exist. Directly use "$rsites
    fi
fi

[ ! -d "$output_dir/clean_fastq" ] && mkdir $output_dir/clean_fastq
[ ! -d "$output_dir/bwa_output" ] && mkdir $output_dir/bwa_output
[ ! -d "$output_dir/parse_temp" ] && mkdir $output_dir/parse_temp

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d %H:%M:%S')] Start cleaning fastq ..."

# R1_str=''
# for i in ${R1[@]};do
#     R1_str=$R1_str' '$i
# done

# R2_str=''
# for i in ${R2[@]};do
#     R2_str=$R2_str' '$i
# done

imargi_clean.sh -1 $R1_str -2 $R2_str -N $base_name -o $output_dir/clean_fastq -t $threads

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d %H:%M:%S')] Start BWA mem mapping ..."

bwa mem -t $threads -SP5M $bwa_index \
    $output_dir/clean_fastq/clean_${base_name}_R1.fastq.gz \
    $output_dir/clean_fastq/clean_${base_name}_R2.fastq.gz \
    2>$output_dir/bwa_output/bwa_log_$base_name.txt | \
    samtools view -@ $threads -Shb - > $output_dir/bwa_output/$base_name.bam

bwa_flag=$( tail -n 3 $output_dir/bwa_output/bwa_log_$base_name.txt | 
    awk 'BEGIN{flag="ERROR"}{if($0~ /^\[main\] Real time/){flag="OK";}}END{print flag;}' )

if [[ "$bwa_flag" == "OK" ]]; then
    tail -n 3 $output_dir/bwa_output/bwa_log_$base_name.txt
    echo "BWA mapping finished."
else
    echo "BWA mapping ERROR."
    tail -n 3 $output_dir/bwa_output/bwa_log_$base_name.txt
    echo "Check more error log info in  "$output_dir/bwa_output/bwa_log_$base_name.txt
    echo "It's usually caused by out of memory or incorrect index files."
    echo "!!!!!!!!!!!!!!!! Pipeline exit with BWA Error at $(date '+%Y-%m-%d %H:%M:%S %Z') !!!!!!!!!!!!!!!!" && usage
fi

echo ">>>>>>>>>>>>>>>>[$(date '+%m-%d %H:%M:%S')] Start parsing, deduplicating, and filtering RNA-DNA interactions ..."
imargi_parse.sh -r $ref_name -c $chromsize -R $rsites -b $output_dir/bwa_output/$base_name.bam \
    -o $output_dir -Q $mapq -G $gap -O $offset -M $max_ligation_size -d false -D $output_dir/parse_temp -t $threads

echo "==================== iMARGI Pipeline END at at $(date '+%Y-%m-%d %H:%M:%S %Z') ===================="
