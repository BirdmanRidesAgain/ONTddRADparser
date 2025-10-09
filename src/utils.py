__all__ = ["enzyme_list", "print_args", "parse_seqfile", "write_seqfile", "make_outdir", "parse_barcodes"]

import gzip
from Bio import SeqIO
from Bio import Restriction
from mimetypes import guess_type
from functools import partial
import subprocess
import shutil
import io
import datetime
import os

enzyme_list=list(Restriction.__dict__)

def print_args(args):
    for key, value in vars(args).items():
        print(f"{key}: {value}")

def parse_seqfile(seqfile, format):
    '''
    Takes in a path to a sequence file (compressed or uncompressed), and a valid string in SeqIO.parse() denoting the file format.
    Returns a list of all records in the file.
    Use case of function generally assumed to be to read fasta/fastq files.
    '''

    encoding = guess_type(seqfile)[1]
    _open = partial(gzip.open, mode = 'rt') if encoding == 'gzip' else open
    with _open(seqfile) as f:
        seqfile_list = list(SeqIO.parse(f, format))
    return seqfile_list

def parse_barcodes(barcodefile):
    '''
    Takes in a path to a Stacks-formatted tsv containing 6-nucleotide barcodes in column1 and sample names in column2.
    A minimal example is presented here:
    ACACCT  R10N00251
    ACAGCA  R20N00088
    ACCTAC  R20N00078

    Returns a dict of containing the barcodes as keys and the individual names as values.
    '''
    result = {}
    with open(barcodefile, "r", encoding="utf-8") as f:
        for i, line in enumerate(f, 1):
            parts = line.rstrip("\n").split("\t")
            print(parts)
            if (len(parts) != 2 ):
                raise ValueError(f"Malformed line {i}: {line.strip()}")
            valid_DNA_nucelotides={'A','T','C','G'}
            if (len(parts[1]) != 6  or (not set(parts[1]).issubset(valid_DNA_nucelotides))):
                    raise ValueError(f"Malformed barcode {i}: {line.strip()}")
            key, value = parts
            result[key] = value
    return result

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
            SeqIO.write(seqs, handle, "fastq")
            p.stdin.close()
        p.wait()

def make_outdir(prefix):
    """
    Script makes an output directory in your working directory from a prefix and current timestamp.
    It then returns the path of that directory to the main script.
    """
    timestamp=datetime.datetime.now().strftime("%Y_%b_%d_%H%M%S")
    outdir=f"{prefix}_{timestamp}"
    os.makedirs(outdir, exist_ok=True)
    return(outdir)

def add(int1,int2):
    return(int1+int2)