#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import re
import sys

####################################################################################################

parser = argparse.ArgumentParser(description='Convert IMGT V-REGION reference files into IMSEQ reference files')
parser.add_argument('infile', metavar='<input file>', type=argparse.FileType('r'), help='Input file (IMGT with gaps)')
parser.add_argument('pos', metavar='<triplet pos>', type=int, help='Cys / Phe triplet position (within gaps)')

opts = parser.parse_args()

####################################################################################################

def parse_fasta_id(s):
    new_id = s.split('|')[1]
    imseq_id = new_id[0:3] + '|' + new_id[3] + '|';
    new_id = new_id[4:].split('*')
    imseq_id = imseq_id + new_id[0] + '|' + new_id[1].zfill(2)
    return imseq_id

for faRec in SeqIO.parse(opts.infile, 'fasta'):
    if len(faRec.seq) < opts.pos:
        print('Skipping "' + faRec.id + '", triplet not included', file = sys.stderr)
        continue
    prefix = str(faRec.seq[0:opts.pos-1])
    suffix = str(faRec.seq[(opts.pos-1):])

    n_gaps = int(0)

    # Count gaps at beginning of suffix (shouldn't be any)
    m = re.match('^[.]+', suffix)
    if m:
        n_gaps = n_gaps + len(m[0])

    # Count gaps throughout prefix
    n_gaps = n_gaps + len(re.sub("[^.]", "", prefix))

    print('>' + parse_fasta_id(faRec.id) + '|' + str(opts.pos - n_gaps - 1))
    print(re.sub('[.]', '', str(faRec.seq)))

