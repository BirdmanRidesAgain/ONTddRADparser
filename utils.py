__all__ = ["print_args", "parse_seqfile", "write_seqfile"]

import gzip
from Bio import SeqIO
from mimetypes import guess_type
from functools import partial
import subprocess
import shutil
import io

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