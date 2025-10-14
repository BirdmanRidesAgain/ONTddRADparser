__all__ = ["enzyme_lst", "print_args", "parse_seqfile", "write_seqfile", "make_outdir", "initialize_df","parse_ONT_demux_file"]

import gzip
from Bio import Seq
from Bio import SeqIO
from Bio import Restriction
import numpy as np
from mimetypes import guess_type
from functools import partial
import subprocess
import shutil
import io
import os
import pandas as pd

enzyme_lst=list(Restriction.__dict__)

def print_args(args):
    print("User-defined arguments:")
    for key, value in vars(args).items():
        print(f"\t{key}: {value}")

def parse_seqfile(seqfile, format):
    '''
    Takes in a path to a sequence file (compressed or uncompressed), and a valid string in SeqIO.parse() denoting the file format.
    Returns a list of all records in the file.
    Use case of function generally assumed to be to read fasta/fastq files.
    '''

    encoding = guess_type(seqfile)[1]
    _open = partial(gzip.open, mode = 'rt') if encoding == 'gzip' else open
    with _open(seqfile) as f:
        seqfile_lst = list(SeqIO.parse(f, format))
    return seqfile_lst

def parse_ONT_demux_file(filepath):
    '''
    Reads a 5-column TSV file expected to contain a header and checks that columns 2-5 contain only valid DNA nucleotides (A, T, C, G).
    
    A minimal example looks like this:
    individual	index_full	index	barcode_full	barcode
    R10N00251	CAA...AATT	CGTGAT	AAT...GCA	ACACCT
    R20N00088	CAA...ATTT	CGTGAT	AAT...GCA	ACAGCA

    Raises a ValueError if either the barcode or index are not 6 or 9 nucleotides long.
    Returns a data frame.
    '''
    df = pd.read_csv(filepath, sep='\t', header='infer')
    if df.shape[1] != 5:
        raise ValueError(f"Expected 6 columns, found {df.shape[1]}")
    
    # Check columns 3 and 5 (index 2 and 4) for string length 6 or 9
    for col in [df.columns[2], df.columns[4]]:
        invalid_rows = df[~df[col].apply(lambda x: isinstance(x, str) and len(x) in (6, 9))]
        if not invalid_rows.empty:
            raise ValueError(
            f"Barcodes and indices '{col}' must be 6 or 9 nucleotides long. "
            f"Invalid rows:\n{invalid_rows[[col]].to_string(index=True)}"
            )
    # Convert columns 2-5 to sequence objects
    for col in df.columns[1:5]:
        df[col] = df[col].apply(lambda x: Seq.Seq(str(x)) if pd.notnull(x) else x)
    
    return df

def write_seqfile(filename, seqs, format):
    '''
    Takes in a path to an output sequence file, a list of sequences, and and a valid string in SeqIO.parse() denoting the file format.
    The function then uses SeqIO.write and subprocess to compress the seqfile.
    There is no returned item.
    It attempts to use pigz if it is installed, and otherwise uses gzip.
    '''
    compressor = "pigz" if shutil.which("pigz") else "gzip"
    
    with subprocess.Popen([compressor, "-c"], stdin=subprocess.PIPE, stdout=open(filename, "wb")) as p:
        # textIO is needed to turn the binary pigz output into a string SeqIO can parse
        with io.TextIOWrapper(p.stdin, encoding="utf-8") as handle:
            SeqIO.write(seqs, handle, format)
            p.stdin.close()
        p.wait()

def initialize_df(num_rows, column_lst):
    '''
    Allocates space for a data frame.
    Values are given the placeholder value of 'np.NaN' (float) by default.
    Takes an index and column names as arguments, returns a data frame.
    '''
    df = pd.DataFrame(np.nan, index=np.arange(num_rows), columns=column_lst)
    return(df)

def make_outdir(prefix):
    '''
    Script makes an output directory in your working directory from a prefix.
    It then returns the path of that directory to the main script.
    '''
    os.makedirs(prefix, exist_ok=True)
    return(prefix)