__all__ = ["Boundary","ConstructElement","DemuxConstruct","DemuxAlignment"]

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align
import numpy as np

class Boundary:
    '''
    Simple class representing the zero-indexed start/end boundaries of two sequences.
    Alignments are made with BioPython's Align module.
    '''
    def __init__(self, type, orientation, start_idx, end_idx):
        self.type = type
        self.orientation = orientation
        self.start_idx = start_idx
        self.end_idx = end_idx
    def __str__(self):
        str=f"""
        Type {self.type}
        Orientation {self.orientation}
        Alignment start {self.start_idx}
        Alignment end {self.end_idx}
        """
        return(str)

class ConstructElement:
    '''
    Individual elements of a larger DNA construct, implemented as a Seq object with metadata.
    In the context of ONTddRADparse, they correspond to the index, barcode, index_full and barcode_full elements.
    We include the 'type' of construct, (index, barcode, etc.), the alignment specificity (`aln_percent`), and the raw sequence.
    '''
    def __init__(self, type, aln_percent, seq):
        self.type = type
        self.aln_percent = aln_percent
        self.seq = Seq(seq)

    def __str__(self):
        str = f"""ConstructElement type {self.type}
        ConstructElement alignment percentage {self.aln_percent}
        ConstructElement sequence {self.seq}
        """
        return(str)

class DemuxConstruct:
    '''
    Represents the construct sequence data we add to our sequences to demux them.
    Consists of a sample ID and four ConstructElement objects: two 6-9bp ones which must be exact matches, and two longer ones where error is allowed.
    Canonically, should be read in from your demux file.
    '''
    def __init__(self, sample_id, index_full, index, barcode_full, barcode, fuzzy_aln_percent, exact_aln_percent):
        self.sample_id = sample_id
        self.index_full = ConstructElement('index_full', fuzzy_aln_percent, index_full)
        self.index = ConstructElement('index', exact_aln_percent, index)
        self.barcode_full = ConstructElement('barcode_full', fuzzy_aln_percent, barcode_full)
        self.barcode = ConstructElement('barcode', exact_aln_percent, barcode)

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
    Represents the alignments between a SecRecord object and a DemuxConstruct.
    Contains alignments of four ConstructElements (forward and reverse, so 8 in total), and a 'valid' tag.
    Also contains
    '''

    # We add two boundaries for forward and reverse orientations
    uninitialized_boundary = Boundary('undefined','undefined', np.nan, np.nan)
    index_full_boundaries = [uninitialized_boundary, uninitialized_boundary]
    index_boundaries = [uninitialized_boundary, uninitialized_boundary]
    barcode_full_boundaries = [uninitialized_boundary, uninitialized_boundary]
    barcode_boundaries = [uninitialized_boundary, uninitialized_boundary]

    # we set validity true until proven otherwise
    valid=True

    def __init__(self, SeqRecord, construct, aligner):
        self.SeqRecord = SeqRecord
        self.DemuxConstruct = construct
        self.aligner = aligner

    def __str__(self):
        str = f"""Alignment between SeqRecord {self.SeqRecord.id} + Construct {self.DemuxConstruct.sample_id}.
        Construct element\tstart_idx\tend_idx\t
        Full index F\t{self.index_full_boundaries[0].start_idx}\t{self.index_full_boundaries[0].end_idx}
        Full index R\t{self.index_full_boundaries[1].start_idx}\t{self.index_full_boundaries[1].end_idx}
        Short index F\t{self.index_boundaries[0].start_idx}\t{self.index_boundaries[0].end_idx}
        Short index R\t{self.index_boundaries[1].start_idx}\t{self.index_boundaries[1].end_idx}
        Full barcode F\t{self.barcode_full_boundaries[0].start_idx}\t{self.barcode_full_boundaries[0].end_idx}
        Full barcode R\t{self.barcode_full_boundaries[1].start_idx}\t{self.barcode_full_boundaries[1].end_idx}
        Short barcode F\t{self.barcode_boundaries[0].start_idx}\t{self.barcode_boundaries[0].end_idx}
        Short barcode R\t{self.barcode_boundaries[1].start_idx}\t{self.barcode_boundaries[1].end_idx}
    
        Validity of alignment={self.valid}
        """
        return(str)
    
    def align_index_full(self):
        aln_f = align_target(self.SeqRecord.seq, self.DemuxConstruct.index_full.seq, self.aligner, 'f', self.DemuxConstruct.index_full.aln_percent)
        aln_r = align_target(self.SeqRecord.seq, self.DemuxConstruct.index_full.seq, self.aligner, 'r', self.DemuxConstruct.index_full.aln_percent)
        f_boundary = Boundary(self.DemuxConstruct.index_full.type, 'f', aln_f[0], aln_f[1])
        r_boundary = Boundary(self.DemuxConstruct.index_full.type, 'r', aln_r[0], aln_r[1])
        self.index_full_boundaries=[f_boundary,r_boundary]

    def align_index(self):
        aln_f = align_target(self.SeqRecord.seq, self.DemuxConstruct.index.seq, self.aligner, 'f', self.DemuxConstruct.index.aln_percent)
        aln_r = align_target(self.SeqRecord.seq, self.DemuxConstruct.index.seq, self.aligner, 'r', self.DemuxConstruct.index.aln_percent)
        f_boundary = Boundary(self.DemuxConstruct.index.type, 'f', aln_f[0], aln_f[1])
        r_boundary = Boundary(self.DemuxConstruct.index.type, 'r', aln_r[0], aln_r[1])
        self.index_boundaries=[f_boundary,r_boundary]

    def align_barcode_full(self):
        aln_f = align_target(self.SeqRecord.seq, self.DemuxConstruct.barcode_full.seq, self.aligner, 'f', self.DemuxConstruct.barcode_full.aln_percent)
        aln_r = align_target(self.SeqRecord.seq, self.DemuxConstruct.barcode_full.seq, self.aligner, 'r', self.DemuxConstruct.barcode_full.aln_percent)
        f_boundary = Boundary(self.DemuxConstruct.barcode_full.type, 'f', aln_f[0], aln_f[1])
        r_boundary = Boundary(self.DemuxConstruct.barcode_full.type, 'r', aln_r[0], aln_r[1])
        self.barcode_full_boundaries=[f_boundary,r_boundary]


    def align_barcode(self):
        aln_f = align_target(self.SeqRecord.seq, self.DemuxConstruct.barcode.seq, self.aligner, 'f', self.DemuxConstruct.barcode.aln_percent)
        aln_r = align_target(self.SeqRecord.seq, self.DemuxConstruct.barcode.seq, self.aligner, 'r', self.DemuxConstruct.barcode.aln_percent)
        f_boundary = Boundary(self.DemuxConstruct.barcode.type, 'f', aln_f[0], aln_f[1])
        r_boundary = Boundary(self.DemuxConstruct.barcode.type, 'r', aln_r[0], aln_r[1])
        self.barcode_boundaries=[f_boundary,r_boundary]

def align_target(seq, subseq, aligner, orientation, aln_percent):
    '''
    Helper function for `check_seq_for_full_index` and others.
    Takes a Bio.Record object and a Bio.Seq object, and aligns them.
    Returns a tuple of the start and end indices of 'subseq' to 'seq'.
    Returns a list of two indices which can be formatted to a Boundary object.
    
    '''
    aligner.target_end_gap_score = aligner.query_end_gap_score = 0.0
    min_aln_score=len(subseq)*aln_percent

    if (orientation=='f'):
        aln=aligner.align(seq, subseq)
    elif (orientation=='r'):
        aln=aligner.align(seq.reverse_complement(),subseq)
    else:
        raise ValueError(f"Ambiguous alignment orientation.\n\tAccepted values: 'f', 'r'.\n\tActual value: {orientation}")
    
    if (aln.score < min_aln_score):
        return([np.nan, np.nan])
    else:
        aln_boundaries=(aln[0].aligned[0].flatten())
        return([min(aln_boundaries), max(aln_boundaries)])

def init_aligner(open_gap_score=-.5, extend_gap_score=-.1):
    # Generic aligner we'll initialize once
    # We penalize opening gaps because our markers should theoretically be one group
    aligner=Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    return(aligner)





def main():
    # Define sample variables
    id='R10N00251'	
    index_full='CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATT'
    index='CGTGAT'
    barcode_full='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACACCTTGCA'
    barcode='ACACCT'
    seqrec=SeqRecord(Seq("CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATTAATGATACGGCGACGTGATCCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACACCTTGCA"), id='test_sample')

    fuzzy_aln=.9
    exact_aln=1

    # construct relevant class examples
    construct=DemuxConstruct(id, index_full, index, barcode_full, barcode, fuzzy_aln, exact_aln)
    aligner = init_aligner()
    alignment=DemuxAlignment(seqrec, construct, aligner)
    alignment.barcode_boundaries=[Boundary('barcode','f' , 0, -1),Boundary('barcode','r' , 0, -1)]
    print(alignment)

    # testing method for setting attributes
    alignment.align_index_full()
    print(alignment)
    alignment.align_index()
    alignment.align_barcode_full()
    alignment.align_barcode()
    print(alignment)



if __name__=="__main__":
    main()