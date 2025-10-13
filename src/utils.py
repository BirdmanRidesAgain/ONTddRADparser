__all__ = ["enzyme_lst", "print_args", "parse_seqfile", "write_seqfile", "make_outdir", "parse_ONT_demux_file","initialize_data_frame","get_fuzzy_alignments","filter_alignment_by_score"]

import gzip
from Bio import Seq
from Bio import SeqIO
from Bio import Restriction
from Bio import Align
import numpy as np
from mimetypes import guess_type
from functools import partial
import subprocess
import shutil
import io
import os
import pandas as pd

enzyme_lst=list(Restriction.__dict__)

def filter_alignment_by_score(aln, max_aln_score, match_percent):
    '''
    Takes an alignment, a max score and a minimum percent of that max score needed to pass.
    Returns a list - either the alignment indices or the invalid indices [-1,-1].
    '''
    min_aln_score=max_aln_score*match_percent
    if (aln.score < min_aln_score):
        return([-1,-1])
    else:
        # the alignment list is always 2 elements long.
        # the first element contains the indices corresponding to the sequence, which we want
        # taking the min and max gets us the boundary of where the index aligned plus whatever gaps were opened
        seq_aln_boundaries=aln[0].aligned[0].flatten()
        return(min(seq_aln_boundaries), max(seq_aln_boundaries)) # first element of the alignment

def align_target(seq, idx, orientation):
    '''
    Helper function for `check_seq_for_full_index` and others.
    Takes a Bio.Record object and a Bio.Seq object, and aligns them.
    Returns a tuple of the start and end indices of 'idx' to 'seq'.
    '''
    # We penalize opening gaps because our markers should theoretically be one group
    aligner=Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    if (orientation=='f'):
        alignment=aligner.align(seq, idx)
    elif (orientation=='r'):
        alignment=aligner.align(seq.reverse_complement(),idx)
    else:
        raise ValueError(f"Ambiguous alignment orientation.\n\tAccepted values: 'f', 'r'.\n\tActual value: {orientation}")
    return(alignment)

def get_fuzzy_alignments(seq_lst, subseq_lst, percent_match):
    '''
    Takes a list of full indexes and sequences, aligns them and returns a list of the alignments.
    Function searches both the forward and reverse complement, resulting in an output list structured like: 
        [seq1,idx,orientation,(subseq_start, subseq_end)]
        [seq2,idx,orientation,(subseq_start, subseq_end)]
        [seq3,idx,orientation,(subseq_start, subseq_end)]
    '''
    # for every sequence/index combination, align and find the indices of all high-qual alignments

    subseq_loc=[]
    for i in seq_lst:
        for j in subseq_lst:
            max_aln_score=len(j)
            for k in ['f','r']:
                alignment=align_target(i, j, k)
                subseq_boundary_lst=[]
                subseq_boundary_lst=filter_alignment_by_score(alignment, max_aln_score, percent_match)
                subseq_loc.append([i.name,j.__str__,k, subseq_boundary_lst])
    return(subseq_loc)






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