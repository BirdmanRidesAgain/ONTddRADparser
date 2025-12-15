#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
from os import path
from tqdm import tqdm
import multiprocessing as multi
from src.utils import *
from src.classes import *
#import time

def main():
    ### DEFINE AND CHECK ARGS
    default_prefix=path.basename(__file__).split(sep='.')[0]
    parser = ArgumentParser(description="Takes a set of ONT-primer prefixes dddRADseq files and a set of barcodes and demultiplexes them.")
    parser.add_argument('-f','--fastq', help='Path to the input fq file. Required.', type=str, required=True)
    parser.add_argument('-d','--demux', help='Path to the demux file. Required.', type=str, required=True)
    parser.add_argument('-b','--buffer', help='The integer number of base pairs outside of the long element the short element can align to. Defaults to 0 (internal matching only).', type=int, default=0)
    parser.add_argument('-t','--threads', help='The number of threads (actually CPU cores) this script can use to parallelize processes. Defaults to half the capacity of the host machine.', type=int, default=multi.cpu_count() // 2)    
    parser.add_argument('-p','--prefix', help='Output prefix. Defaults to \'ONTddRADparse_out\' if not set.', type=str, default=f'{default_prefix}_out')
    parser.add_argument('-fa','--fuzzy_aln_percent', help='The minimum percent identity needed to fuzzy-match a full index to a sequence.', default=.9, type=float)
    parser.add_argument('-ea','--exact_aln_percent', help='The minimum percent identity needed to exact-match a short index to a sequence.', default=1, type=float)

    args = parser.parse_args()
    print_args(args)

    ### PARSE IN FILES
    seq_record_lst = parse_seqfile(args.fastq) # uses FastqGeneralIterator to read big FAs cheaply
    Demux_df = parse_demux_file(args.demux)
    DC_lst = convert_demux_df_to_DemuxConstruct_lst(Demux_df, args.fuzzy_aln_percent, args.exact_aln_percent, args.buffer)
    # fixme - ensure that all DemuxConstruct.sample_ids in this list are unique

    ### init aligner to avoid having to recreate it every time we call DemuxAlignment
    aligner=init_aligner()
    filter_lst = [
    'check_all_ConstructElementAlignments_validity',
    'check_all_ConstructElementAlignments_concatamer_validity',
    'check_all_ConstructElementAlignmentPairs_validity',
    'check_DemuxConstructAlignment_validity',
    ]


    # init empty lists to keep track of fates of various metrics
    # separate this list into two.
    # one loop creates all the DCAs and runs in parallel
    # the other assesses them for valid/invalid, and can run sequentially
    seq_record_fate_lst = []
    DCA_lst_valid = []
    DCA_lst_invalid = []

    for seq_record in tqdm(seq_record_lst):
        for DC in DC_lst:
            # these three remove all reads which are missing elements, or have individual elements present
            DCA=DemuxConstructAlignment(seq_record, DC, aligner)
            for filter in filter_lst:
                getattr(DCA, filter)()
                if not DCA.valid:
                    DCA_lst_invalid.append(DCA)
                    seq_record_fate_lst.append(['fail', DC.sample_id, seq_record.id, filter])
                    break
            # trim DCA I guess
            if DCA.valid:
                DCA.trim_ConstructElements_from_SimpleSeqRecord()
                DCA_lst_valid.append(DCA)
                seq_record_fate_lst.append(['success', DC.sample_id, seq_record.id, 'all_checks_valid'])
                break

    # Create one DemuxxedSample for each unique sample_id
    DS_lst = []
    DS_dict = {}
    for sample_id in Demux_df['sample_id'].unique():
        DS = DemuxxedSample(sample_id)
        DS_lst.append(DS)
        DS_dict[sample_id] = DS

    # Scan through all DemuxConstructAlignment objects and gather SimpleSeqRecords
    for DCA in DCA_lst_valid:
        sample_id = DCA.DemuxConstruct.sample_id
        if sample_id in DS_dict:
            DS_dict[sample_id].gather_SimpleSeqRecords_from_DemuxConstructAlignment(DCA)

    # Create output directory to write outputs to
    outdir = make_outdir(args.prefix)

    # write the fates of all sequences + a plot of them to 'outdir'
    barplot = calc_SimpleSeqRecordFates_stats(seq_record_fate_lst, outdir)


    # and write FastqFiles for each demuxxed sample
    for DS in DS_lst:
        fastq_file = DS.init_FastqFile_from_Demuxxed_Sample(outdir=outdir)
        fastq_file.write_FastqFile_to_outdir()


# If this is being imported
if __name__=="__main__":
    main()