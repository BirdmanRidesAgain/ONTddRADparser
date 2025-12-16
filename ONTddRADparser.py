#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
from os import path
from tqdm import tqdm
import multiprocessing as multi
from itertools import product
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
    print("Reading in sequences")
    SimpleSeqRecord_lst = parse_seqfile(args.fastq) # uses FastqGeneralIterator to read big FAs cheaply
    # we want to chunk up the seq_record_lst into smaller items.

    DC_dict = get_DC_dict(args.demux, args.fuzzy_aln_percent, args.exact_aln_percent, args.buffer)

    ### init aligner to avoid having to recreate it every time we call DemuxAlignment
    aligner=init_aligner()

    # create the inputs to make DemuxConstructAlignments in parallel

    print("Making alignments")
    input_lst = list(product(SimpleSeqRecord_lst, [DC_dict], [aligner]))
    pool = multi.Pool(processes = args.threads)
    DCA_lst = pool.map(make_DCA, tqdm(input_lst))
    pool.close()
    pool.join()

    SimpleSeqRecord_fate_lst = []
    DCA_lst_valid = []
    DCA_lst_invalid = []

    print("Checking alignment validity")
    for DCA in tqdm(DCA_lst):
        if DCA.valid:
            # trim DCA if valid I guess
            DCA.trim_ConstructElements_from_SimpleSeqRecord()
            DCA_lst_valid.append(DCA)
            SimpleSeqRecord_fate_lst.append(['success', DCA.DemuxConstruct.sample_id, DCA.SimpleSeqRecord.id, 'all_checks_valid'])

        else:
            DCA_lst_invalid.append(DCA)
            SimpleSeqRecord_fate_lst.append(['fail', DCA.DemuxConstruct.sample_id, DCA.SimpleSeqRecord.id, filter])


    # Create one DemuxxedSample for each unique sample_id
    DS_lst = []
    DS_dict = {}
    for sample_id in DC_dict.keys():
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
    barplot = calc_SimpleSeqRecordFates_stats(SimpleSeqRecord_fate_lst, outdir)


    # and write FastqFiles for each demuxxed sample
    for DS in DS_lst:
        fastq_file = DS.init_FastqFile_from_Demuxxed_Sample(outdir=outdir)
        fastq_file.write_FastqFile_to_outdir()


# If this is being imported
if __name__=="__main__":
    main()