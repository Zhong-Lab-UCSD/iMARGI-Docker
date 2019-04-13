#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
    Usage: $PROGNAME [-D <distance_type>] [-d <distance_threshold>] [-F <deal_with_filter>]
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
         'filterOut_'. Default is 'output'.
    -i : Input .pairs.gz file.
    -o : Output .pairs.gz file.
    -h : Show usage help
EOF
    exit 1
}

while getopts :D:d:F:i:o:h opt; do
    case $opt in
        D) dis_type=${OPTARG};;
        d) dis_thres=${OPTARG};;
        F) filter_flag=${OPTARG};;
        i) input_file=${OPTARG};;
        o) output_file=${OPTARG};;
        h) usage;;
    esac
done

[ -z "$dis_type" ] && dis_type="5end"
if [[ "$dis_type" != "5end" ]] && [[ "$dis_type" != "inner" ]] && [[ "$dis_type" != "outer" ]]; then
    echo "Error!! Only '5end', 'inner' or 'outer' is acceptable for -D." && usage
fi
[ -z "$dis_thres" ] && dis_thres=200000
if ! [[ "$dis_thres" =~ ^[0-9]+$ ]]; then
    echo "Error!! Only integer number is acceptable for -d" && usage 
fi
[ -z "$filter_flag" ] && filter_flag="output"
if [[ "$filter_flag" != "output" ]] && [[ "$filter_flag" != "drop" ]]; then
    echo "Error!! Only 'output' or 'drop' is acceptable for -F." && usage
fi
[ ! -f "$input_file" ] && echo "Error!! Input file not exist: "$input_file && usage
[ -f "$output_file" ] && echo "Error!! Output file already exist: "$output_file && usage

filterOut_file=$(dirname "$output_file")"/filterOut_"$(basename "$output_file")
tmp_output=$output_file"_"$RANDOM
tmp_filterOut=$filterOut_file"_"$RANDOM
commandline=$PROGNAME" -D "$dis_type" -d "$dis_thres" -F "$filter_flag" -i "$input_file" - o "$output_file

echo ">>>>>>>>>> Start: Filter "$input_file" based on "$dis_type" genomic distance ..."
date

zcat $input_file |\
    awk  -v dis_type="$dis_type" -v filter_flag="$filter_flag" -v dis_thres="$dis_thres" \
    -v commandline="$commandline" -v output_file="$tmp_output" -v filterOut_file="$tmp_filterOut" \
    'BEGIN{OFS="\t"}{
        if(NR % 1000000 == 0){print NR" records processed ..."}
        if($0 ~ /^#/){
            if($0 ~/^#columns/){
                print "#samheader: @PG ID:imargi_distfilter.sh\tPN:imargi_distfilter.sh\tCL:"commandline > output_file;
                print $0 > output_file;
                if(filter_flag=="output"){
                    print "#samheader: @PG ID:imargi_distfilter.sh\tPN:imargi_distfilter.sh\tCL:"commandline > filterOut_file;
                    print $0 > filterOut_file;
                }
            }else{
                print $0 > output_file;
                if(filter_flag=="output"){
                    print $0 > filterOut_file;
                }
            }
        }else{
            if($2!=$4){
                print $0 > output_file;
                next;
            };
            if(dis_type=="5end"){
                if($3-$5 >= dis_thres || $3-$5 <= -dis_thres){
                    print $0 > output_file;
                }else{
                    if(filter_flag=="output"){
                        print $0 > filterOut_file;
                    }
                };
                next;
            };
            type_n1=split(gensub("[0-9]+", "", "g", $11), type1, "");
            val_n1=split($11, val1, "[A-Z=]"); 
            m_flag=0; m_index=0;
            for(i=1;i<=type_n1;i++){
                if(type1[i] ~ /[M=XND]/){
                    if(m_flag==0){
                        m_flag=1; 
                        m_index+=1;
                        mval1[m_index]=val1[i];
                    }else{
                        mval1[m_index]+=val1[i];
                    }
                }else{
                    if(type1[i] ~ /[SH]/){
                        m_flag=0;
                    }
                }
            };
            type_n2=split(gensub("[0-9]+", "", "g", $12), type2, "");
            val_n2=split($12, val2, "[A-Z=]");
            m_flag=0; m_index=0;
            for(i=1;i<=type_n2;i++){
                if(type2[i] ~ /[M=XND]/){
                    if(m_flag==0){
                        m_flag=1;
                        m_index+=1;
                        mval2[m_index]=val2[i];
                    }else{
                        mval2[m_index]+=val2[i];
                    }
                }else{
                    if(type2[i] ~ /[SH]/){
                        m_flag=0;
                    }
                }
            };
            if($6=="+"){
                start1 = $3; end1 = start1 + mval1[1] - 1;
            }else{
                end1 = $3; start1 = end1 - mval1[length(mval1)] + 1;
            };
            if($7=="+"){
                start2 = $5; end2 = start2 + mval2[1] - 1;
            }else{
                end2=$5; start2 = end2 - mval2[length(mval2)] + 1;
            };
            if(start1 <= start2){
                dist_outer=end2-start1;
                dist_inner=start2-end1;
            }else{
                dist_outer=end1-start2;
                dist_inner=start1-end2;
            }
            if(dis_type=="inner"){
                if(dist_inner >= dis_thres){
                    print $0 > output_file;
                }else{
                    if(filter_flag=="output"){
                        print $0 > filterOut_file;
                    }
                };
                next;
            };
            if(dis_type=="outer"){
                if(dist_outer >= dis_thres){
                    print $0 > output_file;
                }else{
                    if(filter_flag=="output"){
                        print $0 > filterOut_file;
                    }
                };
                next;
            };                
        };
    }'

pbgzip -t 0 -c $tmp_output > $output_file
rm $tmp_output

if [[ "$filter_flag" == "output" ]] ; then
    pbgzip -t 0 -c $tmp_filterOut > $filterOut_file
    rm $tmp_filterOut
fi

date
echo "<<<<<<<<<< Finished: Genomic distance filtering."
