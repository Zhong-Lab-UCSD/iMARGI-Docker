#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
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

EOF
    exit 1
}

while getopts :A:a:l:f:C:c:s:m:G:t:i:o:h opt; do
    case $opt in
        A) ant_format=${OPTARG};;
        a) ant_file=${OPTARG};;
        l) ant_level=${OPTARG};;
        f) ant_attr=${OPTARG};;
        C) ant_mode=${OPTARG};;
        c) ant_col=${OPTARG};;
        s) strand_type=${OPTARG};;
        m) min_over=${OPTARG};;
        G) cigar_col=${OPTARG};;
        t) threads=${OPTARG};;
        i) input_file=${OPTARG};;
        o) output_file=${OPTARG};;
        h) usage;;
    esac
done

[ -z "$ant_format" ] && echo "Error!! Please provide annotation file format with -A" && usage
if ! [[ "${ant_format,,}" == "bed" ||  "${ant_format,,}" == "gtf" ]]; then
     echo "Error!! The format $ant_format is not acceptable, only gtf or bed are allowed." && usage
fi

[ ! -f "$ant_file" ] && echo "Error!! Annotation file not exist: "$ant_file && usage
[ ! -f "$input_file" ] && echo "Error!! Input .pairs file not exist: "$input_file && usage
[ -f "$output_file" ] && echo "Error!! Output file already exist: "$output_file && usage

[ -z "$ant_level" ] && echo "Use default annotation level of GTF file: gene body" && ant_level="gene"
[ -z "$ant_attr" ] && echo "Use default annotation attributes of GTF file: gene_id,gene_name,gene_type" \
     && ant_attr="gene_id,gene_name,gene_type"

if ! [[ "${ant_format,,}" == "gtf" ]]; then
     echo "GTF annotation level: "$ant_level
     echo "GTF annotation attributes: "$ant_attr
fi

[ -z "$ant_mode" ] && ant_mode="both" 
if ! [[ "${ant_mode,,}" == "both" ||  "${ant_mode,,}" == "rna" || "${ant_mode,,}" == "dna"  ]]; then
     echo "Error!! Wrong annotation mode $ant_mode! Only both, rna or dna are allowed." && usage
fi
echo "Annotation mode: "$ant_mode

[ -z "$cigar_col" ] && cigar_col="cigar1,cigar2"
[ -z "$ant_col" ] && ant_col="gene1,gene2"
[ -z "$min_over" ] && min_over="1,1"
[ -z "$strand_type" ] && strand_type="rn"
[ -z "$threads" ] && threads=8
if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -t" && usage 
fi

if [[ "${ant_mode,,}" == "both" ]]; then
     if ! [[ "$ant_col" =~ ^[^,]+,[^,]+$ ]]; then
          echo "Error!! Two annotation col names are required for 'both' mode with '-c', separated by ','." && usage 
     fi
     if ! [[ "$min_over" =~ ^[0-9]+,[0-9]+$ ]]; then
          echo "Error!! Two minimum overlap thresholds are required for 'both' mode with '-m'." && usage 
     fi
     if ! [[ "$cigar_col" =~ ^[^,]+,[^,]+$ ]]; then
          echo "Error!! Two cigar info columnss are required for 'both' mode with '-G'." && usage 
     fi
     if ! [[ "$strand_type" =~ ^[srn][srn]$ ]]; then
          echo "Error!! Two characters for strand specific options are required for 'both' mode, such as 'rn'." && usage 
     fi
fi

echo ">>>>>>>>>> Start annotating: ..."
date

imargi_ant.py \
     --ant_format $ant_format \
     --ant_file $ant_file \
     --ant_level $ant_level \
     --ant_attr $ant_attr \
     --ant_mode $ant_mode \
     --ant_col $ant_col \
     --strand_type $strand_type \
     --min_over $min_over \
     --cigar_col $cigar_col \
     --nproc-in $(($threads/3+1)) \
     --nproc-out $(($threads-$threads/3-1)) \
     --output $output_file \
     $input_file

date
echo "<<<<<<<<<< Finished: annotation. "
