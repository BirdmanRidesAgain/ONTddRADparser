__all__ = ["Boundary", "ConstructElement", "DemuxConstruct", "DemuxConstructAlignment", "DemuxxedSample", "FastqFile", "init_aligner"]

import gzip
import subprocess
from io import TextIOWrapper
from shutil import which

import numpy as np
from Bio import Align
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Boundary:
    '''
    Simple class representing the zero-indexed start/end boundaries of two sequences.
    Alignments are made with BioPython's Align module.
    '''

    def __init__(self, type='undefined', orientation='undefined', start_idx=np.nan, end_idx=np.nan):
        self.type = type
        self.orientation = orientation
        self.start_idx = start_idx
        self.end_idx = end_idx

    def __str__(self):
        str=f"""
        Type\tOrientation\tStart_idx\tEnd_index
        {self.type}\t{self.orientation}\t{self.start_idx}\t{self.end_idx}
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
        str = f"""
        type\taln_percent\tseq
        {self.type}\t{self.aln_percent}\t{self.seq}
        """
        return(str)

class ConstructElementAlignment:
    '''
    Represents a ConstructElement aligned to a SeqRecord.
    Generates/houses boundary objects and stores the validity of their alignment.
    A ConstructElementAlignment is valid if XORing RBoundary and RBoundary is true.
    '''
    FBoundary = Boundary()
    RBoundary = Boundary()

    valid = False

    def __init__(self, SeqRecord, ConstructElement, aligner):
        self.SeqRecord = SeqRecord
        self.ConstructElement = ConstructElement
        self.aligner = aligner

    def __str__(self):
        str=f'''
        SeqRecord\tConstructElement_type\taln_valid\tFBoundary_start_idx\tFBoundary_end_idx\tRBoundary_start_idx\tRBoundary_end_idx
        {self.SeqRecord.id}\t{self.ConstructElement.type}\t{self.valid}\t{self.FBoundary.start_idx}\t{self.FBoundary.end_idx}\t{self.RBoundary.start_idx}\t{self.RBoundary.end_idx}
        '''
        return(str)

    def align_ConstructElement(self, orientation):
        '''Aligns the subset to the seq in either the 5->3 ('f') or 3->5 ('r') direction.'''
        aln = align_target(self.SeqRecord.seq, self.ConstructElement.seq, self.aligner, orientation, self.ConstructElement.aln_percent)
        if (orientation == 'f'):
            self.FBoundary=Boundary(self.ConstructElement.type, orientation, aln[0], aln[1])
        elif (orientation == 'r'):
            self.RBoundary=Boundary(self.ConstructElement.type, orientation, aln[0], aln[1])
        else:
            raise ValueError(f"Ambiguous alignment orientation.\n\tAccepted values: 'f', 'r'.\n\tActual value: {orientation}")
        self.check_ConstructElementAlignment_validity()

    def check_ConstructElementAlignment_validity(self):
        # explicitly spelling out conditions for idiot future self
        no_construct_aligned = (np.isnan(self.FBoundary.start_idx)) & (np.isnan(self.RBoundary.start_idx == np.nan))
        construct_aligned_f_and_r = (not np.isnan(self.FBoundary.start_idx)) & (not np.isnan(self.RBoundary.start_idx))
        if (no_construct_aligned | construct_aligned_f_and_r):
            self.valid = False
        else:
            self.valid = True

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

class DemuxConstructAlignment:
    '''
    Represents the alignments between a SecRecord object and a DemuxConstruct.
    Contains alignments of four ConstructElements (forward and reverse, so 8 in total), and a 'valid' tag.
    Also contains
    '''

    # we set validity true until proven otherwise
    valid=False

    def __init__(self, SeqRecord, DemuxConstruct, aligner):
        self.SeqRecord = SeqRecord
        self.DemuxConstruct = DemuxConstruct
        self.index_full_aln=ConstructElementAlignment(SeqRecord, DemuxConstruct.index_full, aligner)
        self.index_aln=ConstructElementAlignment(SeqRecord, DemuxConstruct.index, aligner)
        self.barcode_full_aln=ConstructElementAlignment(SeqRecord, DemuxConstruct.barcode_full, aligner)
        self.barcode_aln=ConstructElementAlignment(SeqRecord, DemuxConstruct.barcode, aligner)

    def __str__(self):

        str = f"""
        SeqRecord\tDemuxConstruct_sample_id\tDemuxConstructAlignment_valid\tindex_full_valid\tindex_valid\tbarcode_full_valid\tbarcode_valid
        {self.SeqRecord.id}\t{self.DemuxConstruct.sample_id}\t{self.valid}\t{self.index_full_aln.valid}\t{self.index_aln.valid}\t{self.barcode_full_aln.valid}\t{self.barcode_aln.valid}
        """
        return(str)

    def check_DemuxConstructAlignment_validity(self):
        
        # check 1 - each ConstructElement must be valid as an element
        ConstructElement_list = [self.index_full_aln.valid, self.index_aln.valid, self.barcode_full_aln.valid, self.barcode_aln.valid]
        if not all(ConstructElement_list):
            self.valid = False
        else:
            self.valid = True

        # check 2 - each ConstructElement must 

    def align_all_ConstructElements(self):
        construct_element_alignments = [self.index_full_aln, self.index_aln, self.barcode_full_aln, self.barcode_aln]
        orientations = ['f', 'r']
        for alignment in construct_element_alignments:
            for orientation in orientations:
                alignment.align_ConstructElement(orientation)
        self.check_DemuxConstructAlignment_validity()

class DemuxxedSample:
    '''
    Represents each individual sample with that has sequence data associated with it.
    Aggregates SeqRecords from DemuxConstructAlignments that share the same sample_id.
    Uses that information to create a FastqFile by calling that class's method.
    '''
    def __init__(self, sample_id):
        self.sample_id = sample_id
        self.SeqRecord_lst = []
    
    def gather_SeqRecords_from_DemuxConstructAlignment(self, DemuxConstructAlignment):
        sample_ids_match = (self.sample_id == DemuxConstructAlignment.DemuxConstruct.sample_id)
        demux_construct_alignment_is_valid = (DemuxConstructAlignment.valid)
        if (sample_ids_match & demux_construct_alignment_is_valid):
            self.SeqRecord_lst.append(DemuxConstructAlignment.SeqRecord)

    def init_FastqFile_from_Demuxxed_Sample(self, outdir='.'):
        f = FastqFile(filename = self.sample_id, outdir=outdir, SeqRecord_lst = self.SeqRecord_lst)
        return(f)
    
class FastqFile:
    '''
    Represents a name and a set of associated sequences.
    Written to the working directory by default, optionally takes an output directory.
    '''
    format = 'fastq'

    def __init__(self, filename, outdir='.', SeqRecord_lst = []):
        if (not filename.endswith('.fq.gz')):
            filename = filename + '.fq.gz'
        self.filename = filename
        if (not outdir.endswith('/')):
            outdir = outdir + '/'
        self.outdir = outdir
        self.filepath = f'{self.outdir}{self.filename}'
        self.SeqRecord_lst = SeqRecord_lst

    def append_SeqRecord_to_FastqFile(self, SeqRecord):
        '''
        Appends SeqRecord objects to SeqRecord_arr attribute.
        Thin wrapper around a list.
        '''
        self.SeqRecord_lst.append(SeqRecord)

    def write_FastqFile_to_outdir(self):
        '''
        Thin wrapper around SeqIO.write.
        Uses `subprocess` to compress files - looks for `pigz` by default; otherwise uses `gzip`.
        '''
        pigz_path = which("pigz")
        if pigz_path:
            with open(self.filepath, "wb") as outfile, \
                    subprocess.Popen([pigz_path, "-c"], stdin=subprocess.PIPE, stdout=outfile) as proc:
                with TextIOWrapper(proc.stdin, encoding="utf-8") as handle:
                    SeqIO.write(self.SeqRecord_lst, handle, self.format)
                    handle.flush()
                proc.stdin.close()
                proc.wait()
                if proc.returncode != 0:
                    raise RuntimeError("pigz failed while writing FASTQ output.")
        else:
            with gzip.open(self.filepath, "wt") as handle:
                SeqIO.write(self.SeqRecord_lst, handle, self.format)









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
    # kept for testing purposes.

    # Define sample variables
    DEFINE_SAMPLES=True
    if DEFINE_SAMPLES:
        id='R10N00251'	
        index_full='CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATT'
        index='CGTGAT'
        barcode_full='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACACCTTGCA'
        barcode='ACACCT'
        seqrec1 = SeqRecord(Seq("CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATTAATGATACGGCGACGTGATCCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACACCTTGCA"), id='test_seq')
        seqrec1.letter_annotations["phred_quality"] = [40] * len(seqrec1.seq) # added to make it a valid fastq
        seqrec2 = SeqRecord(Seq("CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATTAATGATACGGCGACGTGATCCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACACCTTGCA"), id='test_seq')
        seqrec2.letter_annotations["phred_quality"] = [40] * len(seqrec2.seq) # added to make it a valid fastq

    testfile=FastqFile(id)
    testfile.append_SeqRecord_to_FastqFile(seqrec1)
    testfile.append_SeqRecord_to_FastqFile(seqrec2)

    testfile.write_FastqFile_to_outdir()

    fuzzy_aln=.9
    exact_aln=1

    # construct relevant class examples
    construct=DemuxConstruct(id, index_full, index, barcode_full, barcode, fuzzy_aln, exact_aln)
    aligner = init_aligner()
    alignment=DemuxConstructAlignment(seqrec1, construct, aligner)
    alignment.barcode_boundaries=[Boundary('barcode','f' , 0, -1),Boundary('barcode','r' , 0, -1)]

    # testing method for setting attributes
    #print(alignment)


    barcode_element=ConstructElement('barcode', exact_aln, barcode)
    elementalignment=ConstructElementAlignment(seqrec1, barcode_element, aligner)
    elementalignment.align_ConstructElement('f')
    elementalignment.align_ConstructElement('r')
    #print(elementalignment)

if __name__=="__main__":
    main()