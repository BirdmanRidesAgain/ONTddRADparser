#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
import os
from src.utils import *

def main():
    ### DEFINE AND CHECK ARGS
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

    ### PARSE IN FILES
    fq_lst = parse_seqfile(args.fastq, 'fastq')
    demux_df = parse_ONT_demux_file(args.demux)


    ### POPULATE DATA FRAME WITH READS TO BE KEPT/REMOVED
    # initialize a data frame to store information about sequences to be removed/kept
    fq_column_lst=['seq_name', 'full_index_loc_f','full_index_loc_r']
    fq_info_df=initialize_data_frame(len(fq_lst), fq_column_lst)

    ## CODE COMMENTED OUT FOR TESTING PURPOSES BUT OTHERWISE NEEDS TO BE KEPT
    percent_max_aln_score_for_confident_assignment=.95
    # Get alignments for all combinations of reads and full_indices
    unique_full_indices=demux_df['index_full'].unique()
    full_idx_alns_lst_lst=get_fuzzy_alignments(fq_lst, unique_full_indices, percent_max_aln_score_for_confident_assignment)
    print(full_idx_alns_lst_lst[0])

    # Get alignments for all combinations of reads and full_barcodes
    unique_full_barcodes=demux_df['index_full'].unique()
    print(unique_full_barcodes)
    #full_barcode_alns_lst=get_fuzzy_alignments(fq_lst, unique_full_barcodes, percent_max_aln_score_for_confident_assignment)
    #print(full_idx_alignments_lst[0])
    exit(0)



    ### WRITE BINNED READ
    ## OUTPUT FILES TO DIRECTORY
    outdir=make_outdir(args.prefix)


    # Build output directory
    filename="newname.fq.gz"
    fq_path=f"{outdir}/{filename}"
    write_seqfile(fq_path, fq_lst, 'fastq')


# If this is being imported
if __name__=="__main__":
    main()