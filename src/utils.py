__all__ = ["print_args", "parse_seqfile", "make_outdir", "calc_SimpleSeqRecordFates_stats", "plot_SeqRecordFates", "get_DC_dict", "chunk_input_lst"]

import gzip
#from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import numpy as np
from mimetypes import guess_type
from functools import partial
import seaborn as sns
import matplotlib.pyplot as plt

import os
from collections import defaultdict
import pandas as pd
from src.classes import *


def print_args(args):
    print("User-defined arguments:")
    for key, value in vars(args).items():
        print(f"\t{key}: {value}")

def parse_seqfile(filepath: str):
    '''
    Takes in a path to a sequence file (compressed or uncompressed).
    Returns a list of all records in the file as SimpleSeqRecord objects.
    '''
    encoding = guess_type(filepath)[1]
    _open = partial(gzip.open, mode = 'rt') if encoding == 'gzip' else open
    seqfile_lst = []
    with _open(filepath) as f:
        for id, seq, qual in FastqGeneralIterator(f):
            seqfile_lst.append(SimpleSeqRecord(id, seq, qual))
    return seqfile_lst

def parse_demux_file(filepath: str):
    '''
    Takes in a 5-column TSV file expected to contain a header and checks that columns 2-5 contain only valid DNA nucleotides (A, T, C, G).
    
    A minimal example looks like this:
    sample_id	index_full	index	barcode_full	barcode
    R10N00251	CAA...AATT	CGTGAT	AAT...GCA	ACACCT
    R20N00088	CAA...ATTT	CGTGAT	AAT...GCA	ACAGCA

    Raises a ValueError if either the barcode or index are not 6 or 9 nucleotides long.
    Returns a data frame.
    '''
    df = pd.read_csv(filepath, sep='\t', header='infer')
    header=['sample_id','index_full','index','barcode_full','barcode']

    if list(df.columns) != header:
        raise ValueError(
            "Invalid demux header. "
            f"Expected columns: {header}. "
            f"Found: {list(df.columns)}"
        )
    
    if df.shape[1] != 5:
        raise ValueError(f"Expected 5 columns, found {df.shape[1]}")
    
    # Check columns 3 and 5 (index 2 and 4) for string length 6 or 9
    for col in [df.columns[2], df.columns[4]]:
        invalid_rows = df[~df[col].apply(lambda x: isinstance(x, str) and len(x) in (6, 9))]
        if not invalid_rows.empty:
            raise ValueError(
            f"Barcodes and indices '{col}' must be 6 or 9 nucleotides long. "
            f"Invalid rows:\n{invalid_rows[[col]].to_string(index=True)}"
            )
    return df

def convert_demux_df_to_DemuxConstruct_dict(df: pd.DataFrame, fuzzy_aln_percent: float, exact_aln_percent: float, buffer: int):
    '''
    Ingests a data frame and converts it to a dictionary object organized by sample_id.
    The values are a list of DCs, allowing multiple DCs per sample_id.
    '''
    # Build mapping from sample_id to associated barcode/index values
    DC_dict = defaultdict(list)
    # Group the DataFrame by sample_id, and for each group, iterate over all rows (values)
    for sample_id, group_df in df.groupby('sample_id'):
        sample_id_lst = []
        for _, row in group_df.iterrows():
            index_full_CE = ConstructElement(row['index_full'], 'long', fuzzy_aln_percent, buffer)
            index_CE = ConstructElement(row['index'], 'short', exact_aln_percent)
            barcode_full_CE = ConstructElement(row['barcode_full'], 'long', fuzzy_aln_percent, buffer)
            barcode_CE = ConstructElement(row['barcode'], 'short', exact_aln_percent)
            DC = DemuxConstruct(sample_id, index_full_CE, index_CE, barcode_full_CE, barcode_CE)
            sample_id_lst.append(DC)        
        # We append,
        DC_dict[row['sample_id']] = sample_id_lst

    return DC_dict

def get_DC_dict(filepath: str, fuzzy_aln_percent: float, exact_aln_percent: float, buffer: int):
    '''
    Wraps around `parse_demux_file` and `convert_demux_df_to_DemuxConstruct_lst1
    '''
    demux_df = parse_demux_file(filepath)
    DC_dict = convert_demux_df_to_DemuxConstruct_dict(demux_df, fuzzy_aln_percent, exact_aln_percent, buffer)
    return DC_dict

def make_outdir(prefix: str):
    '''
    Script makes an output directory in your working directory from a prefix.
    It then returns the path of that directory to the main script.
    '''
    os.makedirs(prefix, exist_ok=True)
    return(prefix)

def calc_SimpleSeqRecordFates_stats(fate_lst: list, outdir: str):
    '''
    Converts out sequence record information to a data frame, and then plots it.
    Also saves the relevant data frame as a CSV.
    '''
    fate_df = pd.DataFrame(fate_lst, columns=[ 'outcome', 'sample_id', 'seq_id', 'filter'])
    # you need to add the amount of failure in here too.
    # get your count of inds, then rbind in the length of the failures to a new column, and then plot that
    # what you need to do for 'success_df' is get a count of each 

    success_df = fate_df[fate_df['outcome'] == 'success'][['sample_id', 'seq_id']]
    fail_df = fate_df[fate_df['outcome'] == 'fail'][['sample_id', 'seq_id', 'filter']]
    failed_seqs_lst=['NA',len(fail_df['seq_id'].unique())]
    outcome_seqs_lst=list(success_df.groupby(by=['sample_id']).count().itertuples(index=True, name=None))
    outcome_seqs_lst.append(failed_seqs_lst)
    outcome_seqs_df = pd.DataFrame(outcome_seqs_lst, columns=['sample_id', 'count'])

    barplot = plot_SeqRecordFates(outcome_seqs_df)

    outcome_seqs_df.to_csv(f'{outdir}/ONTddRADparser_all_filter_stats.tsv', sep = "\t")
    barplot.savefig(f'{outdir}/ONTddRADparser_demult_success.png')

    return(barplot)

def plot_SeqRecordFates(df: pd.DataFrame):
    '''
    Plots the dataframe from calc_SeqRecordFates_stats using seaborn.
    '''
    fig = plt.figure(figsize=(14,8), layout='constrained')
    ax = sns.barplot(
        df,
        x='sample_id', y='count'
    )
    ax.set_xlabel('Sample IDs (bins)')
    ax.set_ylabel('Number of reads')
    ax.set_title('Successfully demultiplexed reads per sample', loc='left')

    fig.add_axes(ax)
    return fig

def chunk_input_lst(lst: list, n_rds: int = 10000):
    '''
    Splits the incoming input_lst into groups of n iterations.
    We define 'n' arbitrarily as 10000 - I am unsure what the optimal value is.

    We are doing this to prevent 
    '''