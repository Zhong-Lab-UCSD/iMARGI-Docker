#!/usr/bin/env bash
set -e
PROGNAME=$0

usage() {
    cat << EOF >&2
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
EOF
    exit 1
}

while getopts :D:d:F:i:o:h opt; do
    case $opt in
        D) dis_type=${OPTARG};;
        d) dis_thres=${OPTARG};;
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
if ! [[ "$dis_thres" =~ ^[1-9][0-9,]+[0-9]$ ]]; then
    echo "Error!! Only integer number is acceptable for -d" && usage 
fi

[ ! -f "$input_file" ] && echo "Error!! Input file not exist: "$input_file && usage
[ -f "$output_file" ] && echo "Error!! Output file already exist: "$output_file && usage

echo ">>>>>>>>>> Start: Stats "$input_file" based on "$dis_type" genomic distance"
date

pbgzip -c -d $input_file |\
    gawk  -v dis_type="$dis_type" -v dis_thres="$dis_thres" \
    -v input_file="$input_file"  \
    'BEGIN{
        OFS="\t"; 
        total_count=0;
        inter_count=0;
        split(dis_thres, dis_arr, ",");
        for(i in dis_arr){
            prox_count[dis_arr[i]]=0;
            distal_count[dis_arr[i]]=0;
        };
    }{
        if(NR % 1000000 == 0){print NR" records processed ..." > "/dev/stderr" }
        if($0 !~ /^#/){
            total_count+=1;
            if($2!=$4){
                inter_count+=1;
                next;
            }
            if(dis_type=="5end"){
                if($3-$5 >= 0){
                    distance=$3-$5;
                }else{
                    distance=$5-$3;
                };
                for(i in distal_count){
                    if(distance>=int(i)){
                        distal_count[i]+=1;
                    }else{
                        prox_count[i]+=1;
                    };
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
                for(i in distal_count){
                    if(dist_inner>=int(i)){
                        distal_count[i]+=1;
                    }else{
                        prox_count[i]+=1;
                    }
                };
                next;
            };
            if(dis_type=="outer"){
                for(i in distal_count){
                    if(dist_outer>=int(i)){
                        distal_count[i]+=1;
                    }else{
                        prox_count[i]+=1;
                    }
                };
                next;
            };                
        };
    }END{
        print "Stats Report of "input_file" with distal threshold "dis_thres":\n";
        print "Total number of interactions", total_count;
        print "Inter-chromosome interactions", inter_count;
        print "Intra-chromosome interactions", total_count-inter_count;
        print "\n==========: Table of proximal and distal interactions based on thresholds "dis_thres;
        printf "Type";
        for(i in dis_arr){printf "\t"dis_arr[i]};
        printf "\n"; 
        printf "Proximal";
        for(i in prox_count){printf "\t"prox_count[i]};
        printf "\n";
        printf "Distal";
        for(i in distal_count){printf "\t"distal_count[i]}; 
        print "\n==========";
    }' >$output_file


date
echo "<<<<<<<<<< Finished: Stats."
