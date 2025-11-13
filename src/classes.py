__all__ = ["Boundary","Construct","DemuxAlignment"]

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

class Boundary:
    '''
    Simple class representing the zero-indexed start/end boundaries of two sequences.
    Alignments are made with BioPython's Align module.
    '''
    def __init__(self, type, start_idx, end_idx):
        self.type = type
        self.start_idx = start_idx
        self.end_idx = end_idx

class Construct:
    '''
    Represents the construct sequence data we add to our sequences to demux them.
    Consists of a sample ID and four Seq objects: two 6-9bp ones which must be exact matches, and two longer ones where error is allowed.
    Canonically, should be read in from your demux file.
    '''
    def __init__(self, sample_id, index_full, index, barcode_full, barcode):
        self.sample_id = Seq(sample_id)
        self.index_full = Seq(index_full)
        self.index = Seq(index)
        self.barcode_full = Seq(barcode_full)
        self.barcode = Seq(barcode)

    def __str__(self):
        str = f"""Construct object information is:
        Sample id: {self.sample_id}
        Full index: {self.index_full}
        Short index: {self.index}
        Full barcode: {self.barcode_full}
        Short barcode: {self.barcode}"""
        return(str)

class DemuxAlignment:
    '''
    Represents the alignments between a fastq object and its demux Construct.
    Contains four alignments, and a 'valid' tag.
    '''
    uninitialized_boundary = Boundary('undefined',np.nan,np.nan)
    index_full_boundaries=uninitialized_boundary
    index_boundaries=uninitialized_boundary
    barcode_full_boundaries=uninitialized_boundary
    barcode_boundaries=uninitialized_boundary

    valid=True

    def __init__(self, SeqRecord, construct):
        self.SeqRecord = SeqRecord
        self.Construct = construct

    def __str__(self):
        str = f"""Alignment between SeqRecord {self.SeqRecord.id} + Construct {self.Construct.sample_id}.
        Construct element\tstart_idx\tend_idx\t
        Full index\t{self.index_full_boundaries.start_idx}\t{self.index_full_boundaries.start_idx}
        Short index\t{self.index_boundaries.start_idx}\t{self.index_boundaries.end_idx}
        Full barcode\t{self.barcode_full_boundaries.start_idx}\t{self.barcode_full_boundaries.end_idx}
        Short barcode\t{self.barcode_boundaries.start_idx}\t{self.barcode_boundaries.end_idx}
    
        Validity of alignment={self.valid}
        """
        
        return(str)

def main():
    # tests the class
    id='R10N00251'	
    index_full='CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATT'
    index='CGTGAT'
    barcode_full='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACACCTTGCA'
    barcode='ACACCT'

    construct=Construct(id, index_full, index, barcode_full, barcode)
    print(construct)

    seqrec=SeqRecord("ATGC",id='test_sample')
    alignment=DemuxAlignment(seqrec, construct)
    alignment.barcode_boundaries=Boundary('barcode',0,-1)
    print(alignment)


if __name__=="__main__":
    main()