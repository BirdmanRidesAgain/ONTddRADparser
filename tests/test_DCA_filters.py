#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test file designed to check that ConstructElementAlignments work as designed.
"""
# Add src directory to path so we can import modules

from poplib import CR
from re import sub
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
@pytest.fixture
def seq_stem():
    seq_stem=Seq('CGTGTGCTCTTCCGATCTAATTCTTTGGTTTATTTAAGATTAAACGCATACATGTTTGCAGAATTGCTGCCACAATAAGGTACTTTATCTTTTATACACA')
    yield seq_stem

@pytest.fixture
def index_full():
    index_full=Seq('CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATT')
    index_full_rev=Seq('AATTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG')
    yield index_full

@pytest.fixture
def barcode_full():
    barcode_full=Seq('AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACACCTTGCA')
    barcode_full_rev=Seq('TGCAAGGTGTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT')
    yield barcode_full

@pytest.fixture
def index():
    index=Seq('CGTGAT')
    index_rev=Seq('ATCACG')
    yield index

@pytest.fixture
def barcode():
    barcode=Seq('ACACCT')
    barcode_rev=Seq('AGGTGT')
    yield barcode

@pytest.fixture
def aligner():
    aligner = init_aligner()
    yield aligner

@pytest.fixture
def valid_DCA_indexF(seq_stem, index_full, index, barcode_full, barcode, aligner):
    valid_seq_index_F = index_full + seq_stem + barcode_full.reverse_complement()
    testSeqRec = SeqRecord(valid_seq_index_F, id='valid_seq_index_F')
    exact_aln_percent = 1
    fuzzy_aln_percent = 0.9 # to simplify initial testing. We already know that this param works
    exact_match_buffer=0
    long_match_buffer=9

    index_full_CE = ConstructElement(index_full, 'long', fuzzy_aln_percent, long_match_buffer)
    index_CE = ConstructElement(index, 'short', exact_aln_percent, exact_match_buffer)
    barcode_full_CE = ConstructElement(barcode_full, 'long', fuzzy_aln_percent, long_match_buffer)
    barcode_CE = ConstructElement(barcode, 'short', exact_aln_percent, exact_match_buffer)

    index_full_CEA = ConstructElementAlignment(testSeqRec, index_full_CE, aligner)
    index_CEA = ConstructElementAlignment(testSeqRec, index_full_CE, aligner)
    barcode_full_CEA = ConstructElementAlignment(testSeqRec, index_full_CE, aligner)
    barcode_CEA = ConstructElementAlignment(testSeqRec, index_full_CE, aligner)

    index_CEAP = ConstructElementAlignmentPair(index_full_CEA, index_CEA)
    barcode_CEAP = ConstructElementAlignmentPair(barcode_full_CEA, barcode_CEA)    

    DC = DemuxConstruct(
        sample_id='test_sample',index_full=index_full_CE,index=index_CE,barcode_full=barcode_full_CE,barcode=barcode_CE)
    valid_DCA_indexF = DemuxConstructAlignment(testSeqRec, DC, aligner)

    yield valid_DCA_indexF

def test_index_index_full_pair_is_valid(index, index_full, aligner):
    '''
    Each pair of short/long barcodes should be contained within each other.
    '''
    orientationF = align_target(seq=index_full, subseq=index, aligner=aligner, aln_percent=1)
    orientationR = align_target(seq=index_full.reverse_complement(), subseq=index.reverse_complement(), aligner=aligner, aln_percent=1)
    assert orientationF == [24,30] and orientationR == [38,44]

def test_barcode_barcode_full_pair_is_valid(barcode, barcode_full, aligner):
    '''
    Each pair of short/long barcodes should be contained within each other.
    '''
    orientationF = align_target(seq=barcode_full, subseq=barcode, aligner=aligner, aln_percent=1)
    orientationR = align_target(seq=barcode_full.reverse_complement(), subseq=barcode.reverse_complement(), aligner=aligner, aln_percent=1)
    assert orientationF == [58,64] and orientationR == [4,10]

def test_valid_seq_index_orientationF_passes_CEA_checks(valid_DCA_indexF):
    valid_DCA_indexF.check_all_ConstructElementAlignments_validity() 
    assert valid_DCA_indexF.valid

def test_valid_seq_index_orientationF_CEAP_check_finds_orientationF(valid_DCA_indexF):
    valid_DCA_indexF.check_all_ConstructElementAlignments_validity() 
    valid_DCA_indexF.check_all_ConstructElementAlignmentPairs_validity()
    assert valid_DCA_indexF.index_CEAP.orientation == 'F' #and valid_DCA_indexF.barcode_CEAP.orientation == 'R'

def test_valid_seq_index_orientationF_passes_CEA_and_CEAP_checks(valid_DCA_indexF):
    valid_DCA_indexF.check_all_ConstructElementAlignments_validity() 
    valid_DCA_indexF.check_all_ConstructElementAlignmentPairs_validity()
    assert valid_DCA_indexF.valid


def test_valid_seq_index_orientationF_passes_CEA_and_CEAP_checks(valid_DCA_indexF):
    valid_DCA_indexF.check_all_ConstructElementAlignments_validity()
    valid_DCA_indexF.check_all_ConstructElementAlignments_concatamer_validity()
    valid_DCA_indexF.check_all_ConstructElementAlignmentPairs_validity()
    assert valid_DCA_indexF.valid

def test_valid_seq_index_orientationF_passes_all_checks(valid_DCA_indexF):
    valid_DCA_indexF.align_all_ConstructElements()
    valid_DCA_indexF.check_all_ConstructElementAlignments_validity() 
    valid_DCA_indexF.check_all_ConstructElementAlignmentPairs_validity()
    valid_DCA_indexF.check_DemuxConstructAlignment_validity()
    assert valid_DCA_indexF.valid

