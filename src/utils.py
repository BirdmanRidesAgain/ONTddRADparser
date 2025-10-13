__all__ = ["enzyme_lst", "print_args", "parse_seqfile", "write_seqfile", "make_outdir", "parse_ONT_demux_file","initialize_data_frame","get_full_index_boundaries"]

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

def get_full_index_boundaries(seq_lst, idx_lst, n_mismatch=0):
    '''
    Takes a list of full indexes and sequences, aligns them and returns a list containing their start/end indices.
    Function searches both the forward and reverse complement, resulting in an output list structured like: 
        [seq_1,idx,(f_start,f_end),(r_start,r_end)]
        [seq_2,idx,(f_start,f_end),(r_start,r_end)]
        [seq_3,idx,(f_start,f_end),(r_start,r_end)]
    If not found, it returns a -1 for both slots.
    '''
    idx_loc=[]
    for i in seq_lst:
        for j in idx_lst:
            idx_loc.append([i, j, align_target(i, j, 'f'), align_target(i, j, 'r')])
    return(idx_loc)

def align_target(seq, idx, orientation):
    '''
    Helper function for `check_seq_for_full_index` and others.
    Takes a Bio.Record object and a Bio.Seq object, and aligns them.
    Returns a tuple of the start and end indices of 'idx' to 'seq'.
    '''
        #    # we get the start of the full idx with 'find', and then add the length for the last
        #    f_start_boundary=j.seq.find(i)
        #    r_start_boundary=j.seq.reverse_complement().find(i)
        #    idx_boundaries_f=[j.seq.find(i),j.seq.find(i)+index_len]
        #    idx_in_seq_indices=[j.seq.find(i),j.seq.reverse_complement().find(i)]
        #    #for j in idx_in_seq_indices:
        #    #    f_start_idx,r_start_idx=j.seq.find(i),j.seq.reverse_complement().find(i)
    return(-1,-1)


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
    Reads a 6-column TSV file expected to contain a header and checks that columns 2-5 contain only valid DNA nucleotides (A, T, C, G).
    
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

def initialize_data_frame(num_rows, column_lst):
    '''
    Allocates space for a data frame.
    Values are given the placeholder value of 'np.NaN' (float) by default.
    Takes an index and column names as arguments, returns a data frame.
    '''
    df = pd.DataFrame(np.nan, index=np.arange(num_rows), columns=column_lst)
    return(df)






def make_outdir(prefix):
    """
    Script makes an output directory in your working directory from a prefix.
    It then returns the path of that directory to the main script.
    """
    os.makedirs(prefix, exist_ok=True)
    return(prefix)

#def parse_barcodes(file):
#    '''
#    DEPRECATED. 
#    X. Velkeneers would prefer that we use a different approach to demux; see "parse_ONT_demux_file".
#    Takes in a path to a Stacks-formatted tsv containing 6-nucleotide barcodes in column1 and sample names in column2.
#    A minimal example is presented here:
#    ACACCT  R10N00251
#    ACAGCA  R20N00088
#    ACCTAC  R20N00078
#
#    Returns a dict of containing the barcodes as keys and the individual names as values.
#    '''
#    result = {}
#    with open(file, 'r', encoding="utf-8") as f:
#        for i, line in enumerate(f, 1):
#            parts = line.rstrip("\n").split("\t")
#            print(parts)
#            if (len(parts) != 2 ):
#                raise ValueError(f"Malformed line {i}: {line.strip()}")
#            valid_DNA_nucelotides={'A','T','C','G'}
#            if (len(parts[1]) != 6  or (not set(parts[1]).issubset(valid_DNA_nucelotides))):
#                    raise ValueError(f"Malformed barcode {i}: {line.strip()}")
#            key, value = parts
#            result[key] = value
#    return result