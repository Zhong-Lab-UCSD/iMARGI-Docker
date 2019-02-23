#!/usr/bin/env bash
PROGNAME=$0

usage() {
    cat << EOF >&2
    Usage: imargi_annotate.sh [-A <annotation_format>] [-a <annotation_file>] 
                [-l <annotation_level>] [-f <feature_attribute>]
                [-C <end_for_annotate>] [-c <add_col_names>] [-s <strand_specific>]
                [-m <min_overlap>] [-G <cigar>]
                [-t <threads>] [-i <input_file>] [-o <output_file>]

    Dependency: gzip, awk, cool, BEDOPS

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

while getopts :A:a:l:f:C:c:s:m:G:i:o:h opt; do
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
        i) input_file=${OPTARG};;
        o) output_file=${OPTARG};;
        h) usage;;
    esac
done

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
     --output $output_file \
     $input_file
