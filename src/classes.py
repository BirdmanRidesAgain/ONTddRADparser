#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

from altair import Orient

__all__ = ["Boundary", "ConstructElement", "ConstructElementAlignmentPair", "DemuxConstruct", "DemuxConstructAlignment", "DemuxxedSample", "FastqFile", "init_aligner", "align_target","create_sample_data"]

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
    def __init__(self, length: str ='', start_idx=np.nan, end_idx=np.nan, buffer: int = 0, valid: bool = False):
        self.length = length
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.buffer = buffer
        self.valid = valid


    def __str__(self):
        str=f'''
        Length\tOrientation\tStart_idx\tEnd_index\tBuffer
        {self.length}\t{self.valid}\t{self.start_idx}\t{self.end_idx}\t{self.buffer}
        '''
        return(str)

    def get_span(self):
        '''
        Calculates the span of the boundary, taking the buffer into account when present.
        If the buffer puts the span past the start of the seq (ie, if it's less than 0), defaults to 0.
        '''
        span_start_idx = self.start_idx - self.buffer
        if (span_start_idx < 0):
            span_start_idx=0
        span_end_idx = self.end_idx + self.buffer
        span_list=[span_start_idx, span_end_idx]
        return span_list
    
    def set_Boundary(self, length: str, start_idx: int, end_idx: int, buffer: int):
        '''
        Fills in a boundary and checks its validity.
        '''
        self.length=length
        self.start_idx=start_idx
        self.end_idx=end_idx
        self.buffer=buffer
        self.check_Boundary_validity()

    def check_Boundary_validity(self):
        '''Boundaries are valid if both the start and end indices are valid.'''
        if np.isnan(self.start_idx) or np.isnan(self.end_idx):
            self.valid = False
        else:
            self.valid = True

class ConstructElement:
    '''
    Individual elements of a larger DNA construct, implemented as a Seq object with metadata.
    In the context of ONTddRADparse, they correspond to the index, barcode, index_full and barcode_full elements.
    We include the 'length' of construct, (long, short, etc.), the alignment specificity (`aln_percent`), and the raw sequence.
    '''
    def __init__(self, seq: Seq, length: str, aln_percent: float, buffer: int=0):
        self.seq = seq
        self.length = length
        self.aln_percent = aln_percent
        self.buffer = buffer

    def __str__(self):
        str = f'''
        seq\tlength\taln_percent\tbuffer
        {self.seq}\t{self.length}\t{self.aln_percent}\t{self.buffer}
        '''
        return(str)

class ConstructElementAlignment:
    '''
    Represents a ConstructElement aligned to a SeqRecord.
    Generates/houses boundary objects and stores the validity of their alignment.
    A ConstructElementAlignment is valid if XORing RBoundary and RBoundary is true.
    '''

    def __init__(self, SeqRecord: 'SeqRecord', ConstructElement: 'ConstructElement', aligner):
        self.valid = False
        self.SeqRecord = SeqRecord
        self.ConstructElement = ConstructElement
        self.aligner = aligner
        self.orientation = []
        self.FBoundary = Boundary()
        self.RBoundary = Boundary()

    def __str__(self):
        str=f'''
        SeqRecord\tConstructElement_length\Orientation\taln_valid\tFBoundary_start_idx\tFBoundary_end_idx\tRBoundary_start_idx\tRBoundary_end_idx
        {self.SeqRecord.id}\t{self.ConstructElement.length}\t{self.ConstructElement.orientation}\t{self.valid}\t{self.FBoundary.start_idx}\t{self.FBoundary.end_idx}\t{self.RBoundary.start_idx}\t{self.RBoundary.end_idx}
        '''
        return(str)

    def align_ConstructElement(self):
        '''
        Aligns the subset to the seq in the 5->3 ('F') and 3->5 ('R') direction and fills in boundaries.
        Also updates the "orientation" tag of the ConstructElementAlignment.

        '''
        self.set_ConstructElement_FBoundary()
        self.set_ConstructElement_RBoundary()
        
        # We allow both to be appended to account for short sequences possibly appearing more than once by chance
        if self.FBoundary.valid: # aligned in F
            self.orientation.append('F')
        if self.RBoundary.valid: # aligned in R
            self.orientation.append('R')

    def set_ConstructElement_FBoundary(self):     
        '''Align in the forward direction.'''
        aln = align_target(self.SeqRecord.seq, self.ConstructElement.seq, self.aligner, self.ConstructElement.aln_percent)
        self.FBoundary.set_Boundary(self.ConstructElement.length, aln[0], aln[1], self.ConstructElement.buffer)

    def set_ConstructElement_RBoundary(self):
        '''Align in the reverse direction.'''
        aln = align_target(self.SeqRecord.seq.reverse_complement(), self.ConstructElement.seq, self.aligner, self.ConstructElement.aln_percent)
        self.RBoundary.set_Boundary(self.ConstructElement.length, aln[0], aln[1], self.ConstructElement.buffer)

    def check_ConstructElementAlignment_validity(self):
        '''
        Checks the orientation of each boundary, as well as the length of the construct to determine validity.
        Valid either when XORed or when a short element is aligned in both directions.
        '''
        # Base rule: valid when exactly one of aligned_in_F_orientation / aligned_in_R_orientation is True (logical XOR)
        F_orientation = bool(('F' in self.orientation) and ('R' not in self.orientation))
        R_orientation = bool(('F' not in self.orientation) and ('R' in self.orientation))
        xor_valid = bool(F_orientation ^ R_orientation)

        # Additional rule: for short constructs, allow both boundaries aligned
        short_FR_orientation_valid = bool((self.ConstructElement.length == 'short') and ('F' in self.orientation) and ('R' in self.orientation))

        self.valid = xor_valid or short_FR_orientation_valid
        return True

class ConstructElementAlignmentPair:
    '''
    Two complementary `ConstructElements` - one long, and one short. In the context of ONTddRADparser, they consist of:
        - `index` and `index_full`
        - `barcode` and `barcode_full`
    The short one is around 6-9 bp, and is allowed to exist more than once in the SeqRecord.
    The long one should only be found once - otherwise, it is a concatamer.
    '''
    def __init__(self, long_CEA: 'ConstructElementAlignment', short_CEA: 'ConstructElementAlignment'):
        if long_CEA.SeqRecord.id != short_CEA.SeqRecord.id:
            ValueError("ConstructElementAlignments must be based on the same sequence.")
        self.valid = False
        self.long_CEA = long_CEA
        self.short_CEA = short_CEA
        self.orientation = [] # can be ['F','R','invalid']
        self.short_CEA_in_long_CEA = '' # can be T or F
    
    def __str__(self):
        str = f'''
        overall_valid\torientation\tshort_CE_inside_long_CE
        {self.valid}\t{self.orientation}\t{self.short_CEA_in_long_CEA}
        '''
        return(str)   

    def get_ConstructElementAlignmentPair_orientation(self):
        '''
        Determines whether the pair is `F`, `R`, or `invalid`, and updates `self.orientation`.
        '''
        for i in ['F','R']:
            if (i in self.long_CEA.orientation):
                if (i not in self.short_CEA.orientation):
                    self.orientation = [] # this may not be necessary. But it'll return false in any checks
                    return False
                self.orientation.append(i)
                return True

        # It is valid to have more than one orientation in a CEA. However, a CEAP should never have more than one alignment
        if len(self.orientation) > 1:
            raise ValueError("Something has gone wrong here.")

    def check_short_ConstructElementAlignment_in_long_ConstructElementAlignment(self):
        '''
        Performs position checking to ensure that short_CEA is inside long_CEA.
        '''
        if not self.orientation: # if you didn't have an orientation before, get one first.
            self.valid = self.get_ConstructElementAlignmentPair_orientation()

        # now we need to check that the short element is found inside of the long element
        if ('F' in self.orientation):
            short_span = self.short_CEA.FBoundary.get_span()
            long_span = self.long_CEA.FBoundary.get_span()
        elif ('R' in self.orientation): # now we need to check that the short element is found inside of the long element
            short_span=self.short_CEA.RBoundary.get_span()
            long_span=self.long_CEA.RBoundary.get_span()
        else:
            raise ValueError("Something has gone wrong here.")
            
        if ((short_span[0] < long_span[0]) or (short_span[1] > long_span[1])):
            self.short_CEA_in_long_CEA = False
            return False
        self.short_CEA_in_long_CEA = True
        return True

    def check_ConstructElementAlignmentPair_validity(self):
        '''
        Checks if the CEPA is valid, based on self.orientation and self.short_CEA_in_long_CEA, and updates the self.valid tag.
        '''
        if self.get_ConstructElementAlignmentPair_orientation():
            if self.check_short_ConstructElementAlignment_in_long_ConstructElementAlignment():
                self.valid = True
                return True
        self.valid = False
        return False


class DemuxConstruct:
    '''
    Represents the construct sequence data we add to our sequences to demux them.
    Consists of a sample ID and four ConstructElement objects: two 6-9bp ones which must be exact matches, and two longer ones where error is allowed.
    Canonically, should be read in from your demux file.
    '''
    def __init__(self, sample_id, index_full: 'ConstructElement', index:'ConstructElement', barcode_full:'ConstructElement', barcode: 'ConstructElement'):
        self.sample_id = sample_id
        self.index_full = index_full
        self.index = index
        self.barcode_full = barcode_full
        self.barcode = barcode

    def __str__(self):
        str = f'''
        SeqRecord\t{self.sample_id}
        index_full\t{self.index_full}
        index\t\t{self.index}
        barcode_full\t{self.barcode_full}
        barcode\t{self.barcode}'''
        return(str)

class DemuxConstructAlignment:
    '''
    Represents the alignments between a SecRecord object and a DemuxConstruct.
    Contains the SeqRecord and DemuxConstruct.
    Also contains two ConstructElementPairAlignments, (idx and barcode).
    Finally, contains a 'valid' and 'reason' tag.

    '''

    def __init__(self, SeqRecord: 'SeqRecord', DemuxConstruct: 'DemuxConstruct', aligner):
        self.valid = False     # we initialize validity as false until proven otherwise
        self.SeqRecord = SeqRecord
        self.DemuxConstruct = DemuxConstruct

        # Ugh. I really don't like that these are accessible from the top level, but fine; whatever.
        self.index_full_CEA=ConstructElementAlignment(SeqRecord, DemuxConstruct.index_full, aligner)
        self.index_CEA=ConstructElementAlignment(SeqRecord, DemuxConstruct.index, aligner)
        self.barcode_full_CEA=ConstructElementAlignment(SeqRecord, DemuxConstruct.barcode_full, aligner)
        self.barcode_CEA=ConstructElementAlignment(SeqRecord, DemuxConstruct.barcode, aligner)

        self.align_all_ConstructElements()

        self.index_CEAP = ConstructElementAlignmentPair(long_CEA=self.index_full_CEA, short_CEA=self.index_CEA)
        self.barcode_CEAP = ConstructElementAlignmentPair(long_CEA=self.barcode_full_CEA, short_CEA=self.barcode_CEA)


    def __str__(self):

        str = f'''
        SeqRecord\t{self.SeqRecord.id}
        index_CEAP\t{self.index_CEAP}
        barcode_CEAP\t{self.barcode_CEAP}
        '''
        return(str)

    def align_all_ConstructElements(self):
        CEAs = [self.index_full_CEA, self.index_CEA, self.barcode_full_CEA, self.barcode_CEA]
        for CEA in CEAs:
            CEA.align_ConstructElement()

    def check_all_ConstructElementAlignments_validity(self):
        '''
        Wrapper around ConstructElementAlignment.check_ConstructElementAlignment_validity().
        Ports the logic into the main script so we can weed out SeqQbjects that fail the first line of validation.
        '''
        CEAs = [self.index_full_CEA, self.index_CEA, self.barcode_full_CEA, self.barcode_CEA]
        for CEA in CEAs:
            CEA.check_ConstructElementAlignment_validity()
            if not CEA.valid:
                self.valid = False
                return False
            self.valid = True

    def check_all_ConstructElementAlignments_concatamer_validity(self):
        '''
        Checks validated CEs for concatamers
        '''
        long_CEAs = [self.index_full_CEA, self.barcode_full_CEA]
        for long_CEA in long_CEAs:
            long_CEA.check_ConstructElementAlignment_concatamer_validity()
            if not long_CEA.valid:
                self.valid = False
                return False
            self.valid = True        

    def check_all_ConstructElementAlignmentPairs_validity(self):
        '''
        Wrapper around ConstructElementPairAlignment.check_ConstructElementPairAlignment_validity().
        Assumes that all ConstructElements are valid.
        '''
        CEAPs = [self.index_CEAP, self.barcode_CEAP]
        for CEAP in CEAPs:
            CEAP.check_ConstructElementAlignmentPair_validity()
            if not CEAP.valid:
                self.valid = False
                return False
            self.valid = True

    def check_DemuxConstructAlignment_validity(self):
        '''
        Checks if there are any violations in the interaction of the two CEAPs and updates 'self.valid' as necessary.
        Barcode and index must have alternating orientations.
        '''
        for i in [self.index_CEAP, self.barcode_CEAP]:
            if not i.orientation:
                i.get_ConstructElementAlignmentPair_orientation()

        if self.index_CEAP.orientation == self.barcode_CEAP.orientation:
            self.valid = False
            return False
        self.valid = True
        return True

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
        DCA_valid = (DemuxConstructAlignment.valid)
        if (sample_ids_match & DCA_valid):
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

# functions
def align_target(seq: Seq, subseq: Seq, aligner: Align.PairwiseAligner, aln_percent: float = 1):
    '''
    Helper function for `check_seq_for_full_index` and others.
    Takes a Bio.Record object and a Bio.Seq object, and aligns them.
    Returns a tuple of the start and end indices of 'subseq' to 'seq'.
    Returns a list of two indices which can be formatted to a Boundary object.
    '''
    aligner.target_end_gap_score = aligner.query_end_gap_score = 0.0
    min_aln_score=len(subseq)*aln_percent

    aln=aligner.align(seq, subseq)

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

def create_sample_data(seq_name, input_seq: Seq, index_full: Seq, index: Seq, barcode: Seq, barcode_full: Seq):
    '''
    Creates a full contingent of data structures from a single set of complete inputs.
    Inputs should be given as strings.
    '''
    testSeqRec = SeqRecord(input_seq, id = seq_name)
    testSeqRec.letter_annotations["phred_quality"] = [40] * len(testSeqRec.seq) # added to make it a valid fastq

    index_full=Seq(index_full)
    index=Seq(index)
    barcode_full=Seq(barcode_full)
    barcode=Seq(barcode)

    aligner=init_aligner()
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

    # Create DemuxConstruct from sample values
    test_DC = DemuxConstruct(
        sample_id='test_sample',
        index_full=index_full_CE,
        index=index_CE,
        barcode_full=barcode_full_CE,
        barcode=barcode_CE
    )
        # Create DemuxConstructAlignment from test_DC
    test_DCA = DemuxConstructAlignment(testSeqRec, test_DC, aligner)

    return [testSeqRec, test_DC, test_DCA]

def main():

    seq_name='test_seq'
    input_seq='TTTTCTGTCCTGTACTTCGTTCAGTTAGGTATTGTTCAAGCAGAAGACAGAATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATTCTTTGGTTTATTTAAGATTAAACGCATACATGTTTGCAGAATTGCTGCCACAATAAGGTACTTTATCTTTTATACACAGGGAGCATTTACATTTTATACATATGAAAATATACCTCTAAGGACTTTTTTTTTTGCAACAAAACACAGCAGTTACCGACGCCTGATTCCCAGCTGGGGTAAGTCAGCTTGCAGACACTGTAGGAGCTGTGATGGTTGTAGCAGCTGAGATCTAGATGTACAGTCATGTCGAGTTTCCTTAAATATATGACACAAATGTACTACTCTCACACTCAGAGTGGCAACTTAGCACAACTATCTGCCAGCGCATGAGCATCTCTCAGTCCCAGAGAAAACCTTTATCCCTTACCTACACAGGTAATCTTTGAAACACTGACAGCACAACAACTAAACAAAATCTTGTATCATAAAGGCATAGAATTTAGGTTCTTTTTGATGAGAAAATCATCAAATACAGCTCTGAGCCAAAGCCCAGTGATGTTCCAGCCCTTTCTGCTGCCACCACCAGGCATGTCCCATGAAGGCCTGCAGAAGCTAGATAGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAGCAATACGT'
    index_full='CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATT'
    barcode_full='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACACCTTGCA'
    index='CGTGAT'
    barcode='ACACCT'

    seq, DC, DCA=create_sample_data(seq_name, input_seq, index_full, index, barcode_full, barcode)

    print(DC)

    print(DCA)

if __name__=="__main__":
    main()