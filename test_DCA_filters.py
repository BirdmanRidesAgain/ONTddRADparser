#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test file designed to check that ConstructElementAlignments work as designed.
"""
# Add src directory to path so we can import modules

from poplib import CR
import pytest
import sys
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# import custom classes
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from src.classes import ConstructElementAlignment
from src.utils import *
from src.classes import *


# Create testing dataset
# basic parts for dataset
seq_stem=Seq('CGTGTGCTCTTCCGATCTAATTCTTTGGTTTATTTAAGATTAAACGCATACATGTTTGCAGAATTGCTGCCACAATAAGGTACTTTATCTTTTATACACA')

seq_name='input_seq'
index_full=Seq('CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATT')
barcode_full=Seq('AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACACCTTGCA')
index=Seq('CGTGAT')
barcode=Seq('ACACCT')
input_seq = index_full + seq_stem + barcode_full

aligner = init_aligner()
print(align_target(index_full, index, aligner, 1))




seq, DC, DCA = create_sample_data(seq_name, Seq(input_seq), Seq(index_full), Seq(index), Seq(barcode_full), Seq(barcode))




def valid_seq_passes_filters():
    """Basic test to see if valid seq. passes all filters"""
    assert True