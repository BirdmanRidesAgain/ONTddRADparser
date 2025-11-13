#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
import os
from src.utils import *
from src.housekeeping import *
from src.classes import *

def main():
    ### DEFINE AND CHECK ARGS
    parser = ArgumentParser(description="Takes a set of ONT-primer prefixes dddRADseq files and a set of barcodes and demultiplexes them.")
    parser.add_argument('-f','--fastq', help='Path to the input fq file. Required.', type=str, required=True)
    parser.add_argument('-d','--demux', help='Path to the demux file. Required.', type=str, required=True)
    parser.add_argument('-p','--prefix', help='Output prefix. Defaults to \'ONTddRADparse_out\' if not set.', type=str, default=os.path.basename(__file__).split(sep='.')[0])
    parser.add_argument('-e1','--enzyme1', help='The restriction enzyme associated with the index. Valid choices come from BioPython.Restriction', type=str)
    parser.add_argument('-e2','--enzyme2', help='The restriction enzyme associated with the barcode. Valid choices come from BioPython.Restriction', type=str)
    parser.add_argument('-fm','--fuzzy_match_percent', help='The minimum percent identity needed to fuzzy-match a full index to a sequence.', default=.9, type=float)
    parser.add_argument('-em','--exact_match_percent', help='The minimum percent identity needed to exact-match a short index to a sequence.', default=1, type=float)

    args = parser.parse_args()
    #print_args(args)

    if not (args.enzyme1 in enzyme_lst and args.enzyme2 in enzyme_lst and args.enzyme1!=args.enzyme2):
        raise ValueError(f"Invalid enzymes. Enzymes must differ, and valid options are: \n{enzyme_lst}")

    ### PARSE IN FILES
    fq_lst = parse_seqfile(args.fastq, 'fastq')
    demux_construct_list = parse_ONT_demux_file(args.demux)

    for i in demux_construct_list:
        print(i)










    ### POPULATE DATA FRAME WITH READS TO BE KEPT/REMOVED
    # initialize a data frame to store information about sequences to be removed/kept
    fq_column_lst=['seq_name', 'full_index_loc_f','full_index_loc_r']
    fq_info_df=initialize_df(len(fq_lst), fq_column_lst)


    fuzzy_subseqs=['index_full', 'barcode_full']
    testlist=get_valid_seq_subseq_aln_boundaries(fq_lst, demux_construct_list, fuzzy_subseqs, args.fuzzy_match_percent)

    # Code block is always executed; if statement there to explicitly spell out what we're doing
    REMOVE_SEQS_WO_FULL_IDX = REMOVE_SEQS_WO_FULL_BARCODE = True
    if (REMOVE_SEQS_WO_FULL_IDX & REMOVE_SEQS_WO_FULL_BARCODE):
        boundaries_lst=[]
        for col in demux_construct_list[['index_full', 'barcode_full']]:
            unique_subseq=demux_construct_list[col].unique()
            fq_lst, fq_subseq_aln_boundaries_lst = filter_seqs_by_single_subseq_validity(fq_lst, unique_subseq, percent_max_aln_score_for_fuzzy_match)
            boundaries_lst.append(fq_subseq_aln_boundaries_lst)
        boundaries_lst=list(zip(boundaries_lst[0],boundaries_lst[1]))
        print(boundaries_lst[0])


    ###
    # We now have a list of valid seqs and the boundaries of their subseq alignments
    # We need to use the boundaries to search for barcode_





    exit(0)   


    ### WRITE BINNED READ
    ## OUTPUT FILES TO DIRECTORY
    outdir=make_outdir(args.prefix)


    ### SCAN FOR CONCATAMERS

    # Build output directory
    filename="newname.fq.gz"
    fq_path=f"{outdir}/{filename}"
    write_seqfile(fq_path, fq_lst, 'fastq')


# If this is being imported
if __name__=="__main__":
    main()