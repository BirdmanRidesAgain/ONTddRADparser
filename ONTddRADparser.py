#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from argparse import ArgumentParser
from os import path
import subprocess
import sys
from tqdm import tqdm
import multiprocessing as multi
from itertools import product
from src.utils import *
from src.classes import *
#import time

def main():
    # activate your virtual environment
    script_dir = path.dirname(path.abspath(__file__))
    venv_path = path.join(script_dir, '.venv/bin/activate')
    
    if not path.exists(venv_path):
        print(f"Error: Virtual environment activation script not found at {venv_path}", file=sys.stderr)
        sys.exit(1)
    
    command = f'source {venv_path}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error: Failed to activate virtual environment at {venv_path}", file=sys.stderr)
        if result.stderr:
            print(f"Error details: {result.stderr}", file=sys.stderr)
        sys.exit(1)


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

    # Ensure that our alignment percent arguments are a value between 0 and 1
    if not (min(0, 1) < args.fuzzy_aln_percent <= max(0, 1)) or not (min(0, 1) < args.exact_aln_percent <= max(0, 1)):
        raise ValueError("Alignment arguments (-fa, -ea) must be percentages: eg, -fa .9")
    print_args(args)



    ### PARSE IN FILES
    print("Reading in sequences")
    SimpleSeqRecord_lst = parse_seqfile(args.fastq) # uses FastqGeneralIterator to read big FAs cheaply
    # we want to chunk up the seq_record_lst into smaller items.

    sample_id_dict = make_sample_id_dict(args.demux, args.fuzzy_aln_percent, args.exact_aln_percent, args.buffer)

    ### init aligner to avoid having to recreate it every time we call DemuxAlignment
    aligner=init_aligner()

    # create the inputs to make DemuxConstructAlignments in parallel

    print("Making alignments")
    DCA_lst=[]
    chunk_size = 10000 # this is an arbitrary number
    if len(SimpleSeqRecord_lst) > chunk_size:
        print(f"\tLarge input. Running a burnin of {int(chunk_size/10)} replicates to optimize alignment order.")
        sample_id_dict=optimize_sample_id_dict_order(SimpleSeqRecord_lst, int(chunk_size/10), sample_id_dict, aligner)

    input_lst = list(product(SimpleSeqRecord_lst, [sample_id_dict], [aligner]))
    DCA_lst_valid = []

    # the program is going to spend 95% of its time right in this block.
    input_lst_of_lsts = chunk_input_lst(input_lst, chunk_size)
    pool = multi.Pool(processes = args.threads)
    for tranche in tqdm(input_lst_of_lsts):
        DCA_sublst = pool.map(make_DCA, tranche)
        DCA_lst.extend(DCA_sublst)
    pool.close()
    pool.join()


    # now, we get sumstats for all the DCAs we've made
    print("Checking alignment validity")

    sample_id_dict, invalid_dict = split_DCA_lst(sample_id_dict, DCA_lst)

    # Begin writing plots and FQs to outdir
    # Create one DemuxxedSample for each unique sample_id

    for sample_id, sample_id_info in sample_id_dict.items():
        DS = DemuxxedSample(sample_id, sample_id_info[1])
        fastq_file = DS.init_FastqFile_from_Demuxxed_Sample(outdir=f'{args.prefix}/reads')
        fastq_file.write_FastqFile_to_outdir() # One of these classes might be redundant


    # now we want to plot the results from the 'success' dict.
    plot1, df1 = plot_number_of_SimpleSeqRecords_per_sample_id(sample_id_dict)
    plot1.savefig(f'{args.prefix}/{args.prefix}_demult_success.png', dpi=300)
    df1.to_csv(f'{args.prefix}/{args.prefix}_seqs_per_sampleid_stats.tsv', sep = "\t")

    plot2, df2 = plot_reasons_for_SimpleSeqRecord_invalidity(invalid_dict)
    plot2.savefig(f'{args.prefix}/{args.prefix}_demult_failure.png', dpi=300)
    df2.to_csv(f'{args.prefix}/{args.prefix}_failed_seqs_stats.tsv', sep = "\t")

    # deactivate the venv
    result = subprocess.run('deactivate')

# If this is being imported
if __name__=="__main__":
    main()