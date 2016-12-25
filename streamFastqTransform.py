#!/usr/bin/env python3

# uses python version >3.2
# transform a FASTQ into a tab delimited format
# barcode seq quality read_type read_name
# ATACCAGAGT      GTTTTGTCCATCTGAAATTCTAATTTTTCTTATCTTTGTTTTAAACTGATGCTTTTTCAAATTCATTCTTCTATCTTATTTTTAATATGCTTT IDJDGEHHDIJHFJDJGDIDFJJFFDEEEHIHIHGDDGFJDFEGHDFEJDEDDDFHIHIJDHGEGEIIHGGIEHIIEGIHEDFEHFDHEJFEGF        1 @NS500770:1:H5VNJAFXX:3:21501:12626:18995

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
import gzip
import argparse
import os
import sys
import re


def is_gzip(filename):
    magic_bytes = b'\x1f\x8b\x08'
    with open(filename, "rb") as infile:
        file_start = infile.read(len(magic_bytes))
    if file_start.startswith(magic_bytes):
        return True
    return False


def interleave(iter1, iter2):
    for (forward, reverse) in zip(iter1, iter2):
        assert forward.id == reverse.id
        barcode_f = forward.seq[0:10]
        barcode_r = reverse.seq[0:10]
        forward_qual = ''.join([chr(x + 33) for x in forward.letter_annotations['phred_quality']])
        reverse_qual = ''.join([chr(x + 33) for x in reverse.letter_annotations['phred_quality']])
        if barcode_f > barcode_r:
            yield (str(barcode_f+barcode_r), str(forward.seq[12:].rstrip()), str(forward_qual[12:]), 'ab', '1', forward.description)
            yield (str(barcode_f+barcode_r), str(reverse.seq[12:].rstrip()), str(reverse_qual[12:]), 'ab', '2', reverse.description)
        else:
            yield (str(barcode_r+barcode_f), str(reverse.seq[12:].rstrip()), str(reverse_qual[12:]), 'ba', '2', reverse.description)
            yield (str(barcode_r+barcode_f), str(forward.seq[12:].rstrip()), str(forward_qual[12:]), 'ba', '1', forward.description)

if __name__ == "__main__":
    assert sys.version_info >= (3, 2), "Python version 3.2 or greater required."
    parser = argparse.ArgumentParser(description="Transform FASTQ files into SeqConsensus input")
    parser.add_argument('-s1', '-seq1', help="Sequence file 1", required=True)
    parser.add_argument('-s2', '-seq2', help="Sequence file 2", required=True)
    parser.add_argument('-o', '-output', help="Output file", required=True)
    parser.add_argument('-ot', '-output_type', help="Output Type <txt or gz>", default='txt')
    args = vars(parser.parse_args())

    if ('txt' not in args['ot'] and 'gz' not in args['ot']):
        print("Invalid output type.  Only supports 'txt' and 'gz'")
        sys.exit()

    if ('gz' in args['ot'] and ".gz$" not in args['o']):
        args['o'] = args['o'] + '.gz'
    if os.path.isfile(args['o']):
        raise FileExistsError('File exists ' + args['o'])
    if is_gzip(args['s1']):
        records_f = SeqIO.parse(gzip.open(args['s1'], "rt"), "fastq")
        records_r = SeqIO.parse(gzip.open(args['s2'], "rt"), "fastq")
    else:
        records_f = SeqIO.parse(open(args['s1'], "rt"), "fastq")
        records_r = SeqIO.parse(open(args['s2'], "rt"), "fastq")
    if '^gz$' in args['ot']:
        out = gzip.open(args['o'], "wt")
    else:
        out = open(args['o'], "wt")
    for r in interleave(records_f, records_r):
        out.writelines("\t".join(r)+"\n")
    out.close()
