#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
#import pandas as pd
from utils import *


def main():
        # 1. Use argparse to read in the path to the fasta file
    parser = ArgumentParser(description="Takes a set of ONT-primer prefixes dddRADseq files and a set of barcodes and demultiplexes them.")
    parser.add_argument('-f','--fastq', help='Path to the input fq file. Required.', type=str, required=True)
    parser.add_argument('-b','--barcodes', help='Path to the barcodes. Required.', type=str, required=True)
    args = parser.parse_args()

    print_args(args)

    fq_list = parse_seqfile(args.fastq, 'fastq')
    #seqfile_dict =  SeqIO.to_dict(SeqIO.parse(f, format))
    print(type(fq_list))

    write_seqfile('newname.fq.gz', fq_list, 'fastq')


# If this is being imported
if __name__=="__main__":
    main()