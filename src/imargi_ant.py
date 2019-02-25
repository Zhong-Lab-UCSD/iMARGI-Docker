#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# modified from original pairtools_restrict
import io
import sys
import click
import re
import HTSeq
# import numpy as np
# use pandas.read_csv to accelerate loading rsite file, much faster than np.genfromtxt
# import pandas as pd

from pairtools import _fileio, _pairsam_format, cli, _headerops, common_io_options

UTIL_NAME = 'imargi_annotate'

@cli.command()

@click.argument(
    'pairs_path', 
    type=str,
    required=False)

@click.option(
    '-A', '--ant_format',
    type=str,
    required=True,
    help='annotation file format, \'gtf\' or \'bed\'')

@click.option(
    '-a', '--ant_file',
    type=str,
    required=True,
    help='annotation file')

@click.option(
    '-l', '--ant_level',
    type=str,
    default='gene',
    help='gene annotation level for GTF annotation, default is gene')

@click.option(
    '-f', '--ant_attr',
    type=str,
    default='gene_id,gene_name,gene_type',
    help='gene annotation attributes, default is gene_id, gene_name, and gene_type seperated by | in output result')

@click.option(
    '-C', '--ant_mode',
    type=str,
    required=True,
    help='which end to be annotated, \'RNA\', \'DNA\' \'both\'')

@click.option(
    '-c', '--ant_col',
    type=str,
    required=True,
    help='colnames for anntations in output .pairs file')

@click.option(
    '-s', '--strand_type',
    type=str,
    default='rn',
    help='strand-specific settings for RNA and DNA end annotation, default is \'rn\', i.e., reverse-strand-specific \
        for RNA end and ignoring strand info for DNA end.')

@click.option(
    '-m', '--min_over',
    type=str,
    default=1,
    help='minimum bases of overlapping for annotation, 0 means it must be totally inside of the genomic feature')

@click.option(
    '-G', '--cigar_col',
    type=str,
    default='cigar1,cigar2',
    help='CIGAR information columns for annotating, default \'cigar1\' and \'cigar2\'. \'FALSE\' means ignore CIGAR')

@click.option(
    '-o', "--output", 
    type=str, 
    default="", 
    help='output .pairs/.pairsam file.'
        ' If the path ends with .gz/.lz4, the output is compressed by pbgzip/lz4c.'
        ' By default, the output is printed into stdout.')

@click.option(
    '-h', "--help", 
    type=str, 
    default="", 
    help='output .pairs/.pairsam file.'
        ' If the path ends with .gz/.lz4, the output is compressed by pbgzip/lz4c.'
        ' By default, the output is printed into stdout.')

@common_io_options
def annotate(pairs_path, ant_format, ant_file, ant_level, ant_attr, ant_mode, ant_col, strand_type,
    min_over, cigar_col, output, **kwargs):
    '''
    Annotate both RNA, DNA ends with gene annotations in GTF/GFF format or any other genomic
    features in a simple BED file (each line is a named genomic feature).
    '''
    if strand_type == 'n' or strand_type == 'nn':
        stranded = False
    else:
        stranded = True

    if ant_format.lower() == "gtf":
        ant = read_gtf(ant_file, ant_level, ant_attr, stranded)
    else:
        ant = read_bed(ant_file, stranded)

    annotate_pairs(pairs_path, ant, ant_mode, ant_col, strand_type, min_over, cigar_col, output, **kwargs)

def read_gtf(ant_file, ant_level, ant_attr, stranded = True):
    gtf = HTSeq.GFF_Reader(ant_file)
    ant = HTSeq.GenomicArrayOfSets("auto", stranded = stranded)
    ant_attr = ant_attr.split(',')
    for feature in gtf:
        if feature.type == ant_level:
            attrs = []
            for i in ant_attr:
                attrs += [feature.attr[i]]
            ant[feature.iv] += "|".join(attrs)
    return ant

def read_bed(ant_file, stranded = True):
    bed = HTSeq.BED_Reader(ant_file)
    ant = HTSeq.GenomicArrayOfSets("auto", stranded = stranded)
    for feature in bed:
        ant[ feature.iv ] += feature.name
    return ant

def annotate_region(ant, chr_str, pos, strand, match_length, min_over, strand_type):  
    ant_region = set()
    # change 1-based to zero-based genomic coordinate
    if strand == "+":
        start = pos - 1
        end = pos + match_length
    else:
        start = pos - match_length
        end = pos

    min_over = int(min_over)
    if min_over == 0:
        min_over = match_length

    ant_stranded = ant.stranded

    if ant_stranded == False:
        region_iv = HTSeq.GenomicInterval(chr_str, start, end, strand)
        for iv, val in ant[region_iv].steps():
            if iv.length >= min_over:
                ant_region |= val
    else:
        if strand_type == 'r':
            strand = '-' if strand == '+' else '+'
        else:
            strand = '.'

        if strand == '.':
            region_iv = HTSeq.GenomicInterval(chr_str, start, end, '+')
            for iv, val in ant[region_iv].steps():
                if iv.length >= min_over:
                    ant_region |= val
            region_iv = HTSeq.GenomicInterval(chr_str, start, end, '-')
            for iv, val in ant[region_iv].steps():
                if iv.length >= min_over:
                    ant_region |= val
        else:
            region_iv = HTSeq.GenomicInterval(chr_str, start, end, strand)
            for iv, val in ant[region_iv].steps():
                if iv.length >= min_over:
                    ant_region |= val
    if len(ant_region) > 0:
        return ','.join(ant_region)
    else:
        return '.'

def annotate_pairs(pairs_path, ant, ant_mode, ant_col, strand_type, min_over, cigar_col, output, **kwargs):
    instream = (_fileio.auto_open(pairs_path, mode='r', 
                                  nproc=kwargs.get('nproc_in'),
                                  command=kwargs.get('cmd_in', None)) 
                if pairs_path else sys.stdin)

    outstream = (_fileio.auto_open(output, mode='w', 
                                   nproc=kwargs.get('nproc_out'),
                                   command=kwargs.get('cmd_out', None)) 
                 if output else sys.stdout)

    header, body_stream = _headerops.get_header(instream)
    header = _headerops.append_new_pg(header, ID=UTIL_NAME, PN=UTIL_NAME)

    if len(header) == 0:
        sys.stderr.write('.pairs file doesn\'t have header rows!\n')
        raise SystemExit(1)

    col_names = header[-1].split(' ')
    if col_names[0] != '#columns:':
        sys.stderr.write('The last tow of .pairs header is not a valid col_names row (start with \'#columns:\')!\n')
        raise SystemExit(1)       
    col_names.pop(0)

    ant_col = ant_col.split(',')
    for i in ant_col:
        if i in col_names:
            sys.stderr.write('Annotation col names already exist in .pairs file!\n')
            raise SystemExit(1)

    for i in strand_type:
        if i not in ['s', 'r', 'n']:
            sys.stderr.write('Invalid strand specific type for annotation!\n')
            raise SystemExit(1)
    if ant_mode.lower == 'both':
        header[-1] = header[-1] + ' ' + ' '.join(ant_col)
    else:
        header[-1] = header[-1] + ' ' + ant_col[0]

    min_over = [int(i) for i in min_over.split(',')]

    cigar_col = cigar_col.split(',')
    cigar_idx = []
    for i in cigar_col:
        if i not in col_names and i.lower() != 'false':
            sys.stderr.write('Cigar col names doesn\'t exist in .pairs file!\n')
            raise SystemExit(1)
        else:
            cigar_idx += [col_names.index(i)]

    outstream.writelines(l + '\n' for l in header)

    for line in body_stream:
        cols = line.rstrip().split(_pairsam_format.PAIRSAM_SEP)
        
        if ant_mode.lower() == 'rna':
            chrom1, pos1, strand1, cigar1 = cols[_pairsam_format.COL_C1], int(cols[_pairsam_format.COL_P1]), \
                cols[_pairsam_format.COL_S1], cols[cigar_idx[0]]
            if cigar_idx[0] == 'false':
                match_length1 = 1
            else:
                match_length1 = sum([i.ref_iv.length for i in HTSeq.parse_cigar(cigar1)])            
            ant_str = annotate_region(ant, chrom1, pos1, strand1, match_length1, min_over[0], strand_type[0])
        elif ant_mode.lower() == 'dna':
            chrom2, pos2, strand2, cigar2 = cols[_pairsam_format.COL_C2], int(cols[_pairsam_format.COL_P2]), \
                cols[_pairsam_format.COL_S2], cols[cigar_idx[0]]
            if cigar_idx[0] == 'false':
                match_length2 = 1
            else:
                match_length2 = sum([i.ref_iv.length for i in HTSeq.parse_cigar(cigar2)])    
            ant_str = annotate_region(ant, chrom2, pos2, strand2, match_length2, min_over[0], strand_type[0])
        else:
            chrom1, pos1, strand1, cigar1 = cols[_pairsam_format.COL_C1], int(cols[_pairsam_format.COL_P1]), \
                cols[_pairsam_format.COL_S1], cols[cigar_idx[0]]
            if cigar_idx[0] == 'false':
                match_length1 = 1
            else:
                match_length1 = sum([i.ref_iv.length for i in HTSeq.parse_cigar(cigar1)])            
            ant_str = annotate_region(ant, chrom1, pos1, strand1, match_length1, min_over[0], strand_type[0])
            chrom2, pos2, strand2, cigar2 = cols[_pairsam_format.COL_C2], int(cols[_pairsam_format.COL_P2]), \
                cols[_pairsam_format.COL_S2], cols[cigar_idx[1]]
            if cigar_idx[1] == 'false':
                match_length2 = 1
            else:
                match_length2 = sum([i.ref_iv.length for i in HTSeq.parse_cigar(cigar2)])    
            ant_str += _pairsam_format.PAIRSAM_SEP + \
                annotate_region(ant, chrom2, pos2, strand2, match_length2, min_over[1], strand_type[1])
        
        outstream.write(_pairsam_format.PAIRSAM_SEP.join([line.rstrip(), ant_str]))
        outstream.write('\n')
    
    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()


if __name__ == '__main__':
    annotate()