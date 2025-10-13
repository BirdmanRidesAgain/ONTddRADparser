#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
import os
from src.utils import *

def main():
    parser = ArgumentParser(description="Takes a set of ONT-primer prefixes dddRADseq files and a set of barcodes and demultiplexes them.")
    parser.add_argument('-f','--fastq', help='Path to the input fq file. Required.', type=str, required=True)
    parser.add_argument('-d','--demux', help='Path to the demux file. Required.', type=str, required=True)
    parser.add_argument('-p','--prefix', help='Output prefix. Defaults to \'ONTddRADparse_out\' if not set.', type=str, default=os.path.basename(__file__).split(sep='.')[0])
    parser.add_argument('-e1','--enzyme1', help='The restriction enzyme associated with the index. Valid choices come from BioPython.Restriction', type=str)
    parser.add_argument('-e2','--enzyme2', help='The restriction enzyme associated with the barcode. Valid choices come from BioPython.Restriction', type=str)
    args = parser.parse_args()
    #print_args(args)

    if not (args.enzyme1 in enzyme_lst and args.enzyme2 in enzyme_lst and args.enzyme1!=args.enzyme2):
        raise ValueError(f"Invalid enzymes. Valid options are: \n{enzyme_lst}")

    # Parse in files
    fq_lst = parse_seqfile(args.fastq, 'fastq')
    demux_df = parse_ONT_demux_file(args.demux)


    # initialize a data frame to store information about sequences to be removed/kept
    fq_column_lst=['seq_name', 'full_index_loc_f','full_index_loc_r']
    fq_info_df=initialize_data_frame(len(fq_lst), fq_column_lst)


    # remove all reads where an index cannot be found
    unique_full_indices=demux_df['index_full'].unique()
    full_idx_alignments_lst=get_full_index_alignments(fq_lst, unique_full_indices)
    print(full_idx_alignments_lst[0][3].score)
    print(full_idx_alignments_lst[0][3][0].aligned)
    #alignment[0].aligned
    exit(0)



    ## OUTPUT FILES TO DIRECTORY
    outdir=make_outdir(args.prefix)


    # Build output directory
    filename="newname.fq.gz"
    fq_path=f"{outdir}/{filename}"
    write_seqfile(fq_path, fq_lst, 'fastq')


# If this is being imported
if __name__=="__main__":
    main()