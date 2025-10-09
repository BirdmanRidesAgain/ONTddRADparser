#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
import os
#import pandas as pd
from src.utils import *

def main():
    parser = ArgumentParser(description="Takes a set of ONT-primer prefixes dddRADseq files and a set of barcodes and demultiplexes them.")
    parser.add_argument('-f','--fastq', help='Path to the input fq file. Required.', type=str, required=True)
    parser.add_argument('-d','--demux', help='Path to the demux file. Required.', type=str, required=True)
    parser.add_argument('-p','--prefix', help='Output prefix. Defaults to \'ONTddRADparse_out\' if not set.', type=str, default=os.path.basename(__file__).split(sep='.')[0])
    parser.add_argument('-e1','--enzyme1', help='The restriction enzyme associated with the index. Valid choices come from BioPython.Restriction', type=str)
    parser.add_argument('-e2','--enzyme2', help='The restriction enzyme associated with the barcode. Valid choices come from BioPython.Restriction', type=str)
    args = parser.parse_args()
    print_args(args)

    if not (args.enzyme1 in enzyme_list and args.enzyme2 in enzyme_list and args.enzyme1!=args.enzyme2):
        raise ValueError(f"Invalid enzymes. Valid options are: \n{enzyme_list}")

    # Parse in files
    fq_list = parse_seqfile(args.fastq, 'fastq')
    barcode_dict=parse_barcodes(args.demux)
    print(barcode_dict)

    ## Now, we demultiplex all reads by barcode
    print(barcode_dict.keys())
    #for i in fq_list:
    #    print(len(i.seq))
    #    #print(i.seq)

    outdir=make_outdir(args.prefix)

    # Build output directory
    filename="newname.fq.gz"
    fq_path=f"{outdir}/{filename}"
    write_seqfile(fq_path, fq_list, 'fastq')


# If this is being imported
if __name__=="__main__":
    main()