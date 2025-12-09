#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
from os import path
from src.utils import *
from src.classes import *

def main():
    ### DEFINE AND CHECK ARGS
    default_prefix=path.basename(__file__).split(sep='.')[0]
    parser = ArgumentParser(description="Takes a set of ONT-primer prefixes dddRADseq files and a set of barcodes and demultiplexes them.")
    parser.add_argument('-f','--fastq', help='Path to the input fq file. Required.', type=str, required=True)
    parser.add_argument('-d','--demux', help='Path to the demux file. Required.', type=str, required=True)
    parser.add_argument('-b','--buffer', help='The integer number of base pairs outside of the long element the short element can align to. Defaults to 0 (internal matching only).', type=int, default=0)
    parser.add_argument('-p','--prefix', help='Output prefix. Defaults to \'ONTddRADparse_out\' if not set.', type=str, default=f'{default_prefix}_out')
    parser.add_argument('-e1','--enzyme1', help='The restriction enzyme associated with the index. Valid choices come from BioPython.Restriction', type=str)
    parser.add_argument('-e2','--enzyme2', help='The restriction enzyme associated with the barcode. Valid choices come from BioPython.Restriction', type=str)
    parser.add_argument('-fa','--fuzzy_aln_percent', help='The minimum percent identity needed to fuzzy-match a full index to a sequence.', default=.9, type=float)
    parser.add_argument('-ea','--exact_aln_percent', help='The minimum percent identity needed to exact-match a short index to a sequence.', default=1, type=float)

    args = parser.parse_args()
    print_args(args)

    if not (args.enzyme1 in enzyme_lst and args.enzyme2 in enzyme_lst and args.enzyme1!=args.enzyme2):
        raise ValueError(f"Invalid enzymes. Enzymes must differ, and valid options are: \n{enzyme_lst}")

    ### PARSE IN FILES
    seq_record_lst = parse_seqfile(args.fastq, 'fastq')
    demux_df = parse_demux_file(args.demux)
    demux_construct_list = convert_demux_df_to_DemuxConstruct_lst(demux_df, args.fuzzy_aln_percent, args.exact_aln_percent, args.buffer)
    # fixme - ensure that all DemuxConstruct.sample_ids in this list are unique

    ### init aligner to avoid having to recreate it every time we call DemuxAlignment
    aligner=init_aligner()

    # generate all-against-all alignments
    valid_DCA_lst = []
    for seq_record in seq_record_lst:
        for DC in demux_construct_list:

            # these three remove all reads which are missing elements, or have individual elements present
            DCA=DemuxConstructAlignment(seq_record, DC, aligner)
            DCA.check_all_ConstructElementAlignments_validity()
            if DCA.valid:
                DCA.check_all_ConstructElementAlignments_concatamer_validity()
            # if the individual elements are fine, check them collectively
                if DCA.valid:
                    DCA.check_all_ConstructElementAlignmentPairs_validity()
                    if DCA.valid:
                        DCA.check_DemuxConstructAlignment_validity()
            
            # if de DCA is still valid, append the seq to the list.
            if (DCA.valid):
                valid_DCA_lst.append(DCA)



    # Create one DemuxxedSample for each unique sample_id
    demuxxed_sample_lst = []
    demuxxed_sample_dict = {}
    for sample_id in demux_df['sample_id'].unique():
        demuxxed_sample = DemuxxedSample(sample_id)
        demuxxed_sample_lst.append(demuxxed_sample)
        demuxxed_sample_dict[sample_id] = demuxxed_sample
    
    # Scan through all DemuxConstructAlignment objects and gather SeqRecords
    for alignment in valid_DCA_lst:
        sample_id = alignment.DemuxConstruct.sample_id
        if sample_id in demuxxed_sample_dict:
            demuxxed_sample_dict[sample_id].gather_SeqRecords_from_DemuxConstructAlignment(alignment)

    # Create output directory and write FastqFiles for each demuxxed sample
    outdir = make_outdir(args.prefix)
    for demuxxed_sample in demuxxed_sample_lst:
        fastq_file = demuxxed_sample.init_FastqFile_from_Demuxxed_Sample(outdir=outdir)
        fastq_file.write_FastqFile_to_outdir()


# If this is being imported
if __name__=="__main__":
    main()