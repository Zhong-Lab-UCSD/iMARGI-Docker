#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# modified from original pairtools_restrict
import io
import sys
import click
import subprocess
import re

import numpy as np
# use pandas.read_csv to accelerate loading rsite file, much faster than np.genfromtxt
import pandas as pd

from pairtools import _fileio, _pairsam_format, cli, _headerops, common_io_options

UTIL_NAME = 'imargi_restrict'

@cli.command()

@click.argument(
    'pairs_path', 
    type=str,
    required=False)

@click.option(
    '-f', '--frags',
    type=str,
    required=True,
    help='a tab-separated BED file with the positions of restriction fragments '
         '(chrom, start, end). Can be generated using cooler digest.')

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
def restrict(pairs_path, frags, output, **kwargs):
    '''Check and assign restriction fragments to R2 (DNA ends).
    Identify the successfully ligated RNA-DNA molecule.
    New columns, frag2_start, frag2_end, dist2_rsite will be added to the output .pairs file.
    frag2_start and frag2_end are the coordinates of the assigned restriction fragments of R2.
    dist2_rsite is the distance between the 5' end of R2 to the nearest restriction site.
    PAIRS_PATH : input .pairs/.pairsam file. If the path ends with .gz/.lz4, the 
    input is decompressed by pbgzip/lz4c. By default, the input is read from stdin.
    '''
    restrict_py(pairs_path, frags, output, **kwargs)

def restrict_py(pairs_path, frags, output, **kwargs):
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
    if len(header) > 0:
        header[-1] = header[-1] + ' frag1_start frag1_end dist1_rsite frag2_start frag2_end dist2_rsite'
    outstream.writelines((l+'\n' for l in header))

    rfrags=pd.read_csv(frags, delimiter="\t", dtype=None, comment="#", 
        names=['chrom', 'start', 'end'], encoding='utf-8')
    rfrags = rfrags.to_records()

    chrom_borders = np.r_[0,
                          1+np.where(rfrags['chrom'][:-1] != rfrags['chrom'][1:])[0],
                          rfrags.shape[0]]
    rfrags = {rfrags['chrom'][i]:np.insert(rfrags['end'][i:j]+1, 0, 1)
              for i, j in zip(chrom_borders[:-1], chrom_borders[1:])}

    for line in body_stream:
        cols = line.rstrip().split(_pairsam_format.PAIRSAM_SEP)
        # chrom1, pos1 = cols[_pairsam_format.COL_C1], int(cols[_pairsam_format.COL_P1])
        # rfrag_idx1, rfrag_start1, rfrag_end1 = find_rfrag(rfrags, chrom1, pos1)
        chrom1, pos1, strand1, cigar1 = cols[_pairsam_format.COL_C1], int(cols[_pairsam_format.COL_P1]), \
            cols[_pairsam_format.COL_S1], cols[10]
        rfrag_start1, rfrag_end1, dist1_rsite = find_rfrag(rfrags, chrom1, pos1, strand1, cigar1)
        cols += [str(rfrag_start1), str(rfrag_end1), str(dist1_rsite)]
        chrom2, pos2, strand2, cigar2 = cols[_pairsam_format.COL_C2], int(cols[_pairsam_format.COL_P2]), \
            cols[_pairsam_format.COL_S2], cols[11]
        rfrag_start2, rfrag_end2, dist2_rsite = find_rfrag(rfrags, chrom2, pos2, strand2, cigar2)
        cols += [str(rfrag_start2), str(rfrag_end2), str(dist2_rsite)]
        outstream.write(_pairsam_format.PAIRSAM_SEP.join(cols))
        outstream.write('\n')
    
    if instream != sys.stdin:
        instream.close()
    if outstream != sys.stdout:
        outstream.close()

def find_rfrag(rfrags, chrom, pos, strand, cigar):
    if chrom not in rfrags:
         return '!', '!', '!'
    rsites_chrom = rfrags[chrom]
    # idx = min(rsites_chrom.searchsorted(pos, 'right') - 1, len(rsites_chrom) - 2)
    idx = rsites_chrom.searchsorted(pos, 'right') - 1
    dist_left = pos - rsites_chrom[idx]
    dist_right = rsites_chrom[idx + 1] - 1 - pos
    if strand == "+":
        if dist_left <= dist_right:
            return rsites_chrom[idx], rsites_chrom[idx + 1] - 1, dist_left
        else:
            if idx == len(rsites_chrom) - 2:
                return rsites_chrom[idx], rsites_chrom[idx + 1] - 1, dist_left
            else:
                return rsites_chrom[idx + 1], rsites_chrom[idx + 2] - 1, - dist_right
    else:
        if dist_left < dist_right:
            if idx == 0:
                return rsites_chrom[idx], rsites_chrom[idx + 1] - 1, dist_right
            else:
                return rsites_chrom[idx-1], rsites_chrom[idx] - 1,  - dist_left
        else:
            return rsites_chrom[idx], rsites_chrom[idx + 1] - 1, dist_right
        
if __name__ == '__main__':
    restrict()