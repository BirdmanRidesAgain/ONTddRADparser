#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
from tracemalloc import start

__all__ = ["Boundary", "ConstructElement", "ConstructElementAlignmentPair", "DemuxConstruct", "DemuxConstructAlignment", "DemuxxedSample", "FastqFile", "init_aligner", "align_target"]

import gzip
import subprocess
from io import TextIOWrapper
from shutil import which
from itertools import chain 

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
    def __init__(self, seq_len: int, start_idx=np.nan, end_idx=np.nan, buffer: int = 0, valid: bool = False):
        self.seq_len = seq_len
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.buffer = buffer
        self.span = [np.nan, np.nan]
        self.valid = valid

    def __str__(self):
        str=f'''
        Full_seq_len\tAln_start_idx\tAln_end_index\tBuffer\tValidity
        {self.seq_len}\t{self.start_idx}\t{self.end_idx}\t{self.buffer}\t{self.valid}
        '''
        return(str)

    def get_Boundary_span(self):
        '''
        Calculates the span of the boundary, taking the buffer into account when present.
        If the buffer puts the span past the start of the seq (ie, if it's less than 0), defaults to 0.
        If the buffer puts it over the length of the seq, defaults to the seq length.

        If the start or end idx are np.nan, then just leave this unset.
        '''
        if not self.valid:
            raise ValueError("Cannot get span of an invalid alignment.")
        span_start_idx = self.start_idx - self.buffer
        if (span_start_idx < 0):
            span_start_idx=np.int64(0)
        span_end_idx = self.end_idx + self.buffer
        if (span_end_idx > self.seq_len):
            span_end_idx = self.seq_len
        self.span=[span_start_idx, span_end_idx]
        
    def check_Boundary_validity(self):
        '''Boundaries are valid if both the start and end indices are valid.'''
        if np.isnan(self.start_idx) or np.isnan(self.end_idx):
            self.valid = False
        else:
            self.valid = True
    
    def set_Boundary(self, start_idx: int, end_idx: int, buffer: int):
        '''
        Fills in a boundary and checks its validity.
        '''
        self.start_idx=start_idx
        self.end_idx=end_idx
        self.buffer=buffer
        self.check_Boundary_validity()
        if self.valid:
            self.get_Boundary_span()

class ConstructElement:
    '''
    Individual elements of a larger DNA construct, implemented as a Seq object with metadata.
    In the context of ONTddRADparse, they correspond to the index, barcode, index_full and barcode_full elements.
    We include the 'length' of construct, (long, short, etc.), the alignment specificity (`aln_percent`), and the raw sequence.
    '''
    def __init__(self, seq: Seq, CE_type: str, aln_percent: float, buffer: int=0):
        self.seq = seq
        self.CE_type = CE_type
        self.aln_percent = aln_percent
        self.buffer = buffer

    def __str__(self):
        str = f'''
        seq\tCE_ype\taln_percent\tbuffer
        {self.seq}\t{self.CE_type}\t{self.aln_percent}\t{self.buffer}
        '''
        return(str)

class ConstructElementAlignment:
    '''
    Represents a ConstructElement aligned to a SeqRecord.
    Generates/houses boundary objects and stores the validity of their alignment.
    A ConstructElementAlignment is valid if XORing RBoundary and RBoundary is true.
    '''

    def __init__(self, SeqRecord: 'SeqRecord', ConstructElement: 'ConstructElement', aligner):
        self.SeqRecord = SeqRecord
        self.ConstructElement = ConstructElement
        self.aligner = aligner
        self.orientation = []
        seq_len=len(self.SeqRecord.seq)
        self.FBoundary = Boundary(seq_len)
        self.RBoundary = Boundary(seq_len)
        self.valid = False

    def __str__(self):
        str=f'''
        SeqRecord\tConstructElement_type\Orientation\taln_valid\tFBoundary_start_idx\tFBoundary_end_idx\tRBoundary_start_idx\tRBoundary_end_idx
        {self.SeqRecord.id}\t{self.ConstructElement.CE_type}\t{self.ConstructElement.orientation}\t{self.valid}\t{self.FBoundary.start_idx}\t{self.FBoundary.end_idx}\t{self.RBoundary.start_idx}\t{self.RBoundary.end_idx}
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
        # We will only use the alignment method if we have to. Exact matches can get string methods.
        if self.ConstructElement.CE_type == 'long':
            aln = align_target(self.SeqRecord.seq, self.ConstructElement.seq, self.aligner, self.ConstructElement.aln_percent)

        elif self.ConstructElement.CE_type == 'short':
            start_idx = self.SeqRecord.seq.lower().find(self.ConstructElement.seq.lower())
            if start_idx != -1:
                aln = [start_idx, start_idx+(len(self.ConstructElement.seq.lower())-1)]
            else:
                aln = [np.nan, np.nan]
        else:
            raise ValueError("Your CE_type should be either 'short' or 'long'.")

        # set the boundary    
        self.FBoundary.set_Boundary(aln[0], aln[1], self.ConstructElement.buffer)

    def set_ConstructElement_RBoundary(self):
        '''Align in the reverse direction.'''
        if self.ConstructElement.CE_type == 'long':
            aln = align_target(self.SeqRecord.seq.reverse_complement(), self.ConstructElement.seq, self.aligner, self.ConstructElement.aln_percent)
        elif self.ConstructElement.CE_type == 'short':
            start_idx = self.SeqRecord.seq.reverse_complement().lower().find(self.ConstructElement.seq.lower())
            if start_idx != -1:
                aln = [start_idx, start_idx+(len(self.ConstructElement.seq.lower())-1)]
            else:
                aln = [np.nan, np.nan]
        else:
            raise ValueError("Your CE_type should be either 'short' or 'long'.")
        # set the boundary   
        self.RBoundary.set_Boundary(aln[0], aln[1], self.ConstructElement.buffer)

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
        short_FR_orientation_valid = bool((self.ConstructElement.CE_type == 'short') and ('F' in self.orientation) and ('R' in self.orientation))

        self.valid = xor_valid or short_FR_orientation_valid
        return True

    def check_ConstructElementAlignment_concatamer_validity(self):
        '''
        Checks for concatamers in binding in cis to each other.
        Concatamers resulting from the same seq binding in trans are already removed by checkConstructElementAlignment_validity.
        '''
        self.aligner.target_end_gap_score = self.aligner.query_end_gap_score = 0.0
        min_aln_score = len(self.ConstructElement.seq)* self.ConstructElement.aln_percent

        if ('F' in self.orientation):
            aln=self.aligner.align(self.SeqRecord.seq, self.ConstructElement.seq)
        elif ('R' in self.orientation):
            aln=self.aligner.align(self.SeqRecord.seq.reverse_complement(), self.ConstructElement.seq)
        
        if len(aln)>1:
            if aln.score >= min_aln_score:
                self.valid = False
                return False
        self.valid=True
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
        Determines whether the pair is `F`, `R`, or [], and updates `self.orientation`.
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
            short_span = self.short_CEA.FBoundary.span
            long_span = self.long_CEA.FBoundary.span
        elif ('R' in self.orientation): # now we need to check that the short element is found inside of the long element
            short_span=self.short_CEA.RBoundary.span
            long_span=self.long_CEA.RBoundary.span
        else:
            # failing here means you have an invalid alignment.
            self.short_CEA_in_long_CEA = False
            return False
            
        if ((short_span[0] < long_span[0]) or (short_span[1] > long_span[1])):
            self.short_CEA_in_long_CEA = False
            return False
        self.short_CEA_in_long_CEA = True
        return True

    def check_ConstructElementAlignmentPair_validity(self):
        '''
        Checks if the CEPA is valid, based on self.orientation and self.short_CEA_in_long_CEA, and updates the self.valid tag.
        '''
        if self.orientation:
            if self.short_CEA_in_long_CEA:
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
        self.orientation = []

        # construct and validate CEAs
        index_full_CEA=ConstructElementAlignment(SeqRecord, DemuxConstruct.index_full, aligner)
        index_CEA=ConstructElementAlignment(SeqRecord, DemuxConstruct.index, aligner)
        barcode_full_CEA=ConstructElementAlignment(SeqRecord, DemuxConstruct.barcode_full, aligner)
        barcode_CEA=ConstructElementAlignment(SeqRecord, DemuxConstruct.barcode, aligner)

        self.align_all_ConstructElements([index_full_CEA, index_CEA, barcode_full_CEA, barcode_CEA])

        self.index_CEAP = ConstructElementAlignmentPair(long_CEA=index_full_CEA, short_CEA=index_CEA)
        self.barcode_CEAP = ConstructElementAlignmentPair(long_CEA=barcode_full_CEA, short_CEA=barcode_CEA)
        # now get their orientations
        self.index_CEAP.get_ConstructElementAlignmentPair_orientation()
        self.index_CEAP.check_short_ConstructElementAlignment_in_long_ConstructElementAlignment()
        self.barcode_CEAP.get_ConstructElementAlignmentPair_orientation()
        self.barcode_CEAP.check_short_ConstructElementAlignment_in_long_ConstructElementAlignment()

    def __str__(self):

        str = f'''
        SeqRecord\t{self.SeqRecord.id}
        index_CEAP\t{self.index_CEAP}
        barcode_CEAP\t{self.barcode_CEAP}
        '''
        return(str)

    def align_all_ConstructElements(self, CEA_lst):
        for CEA in CEA_lst:
            CEA.align_ConstructElement()
            CEA.check_ConstructElementAlignment_validity()
            if not CEA.valid:
                break

    def check_all_ConstructElementAlignments_validity(self):
        '''
        Wrapper around ConstructElementAlignment.check_ConstructElementAlignment_validity().
        Ports the logic into the main script so we can weed out SeqQbjects that fail the first line of validation.
        '''
        CEA_lst = [self.index_CEAP.long_CEA, self.index_CEAP.short_CEA, self.barcode_CEAP.long_CEA, self.barcode_CEAP.short_CEA]
        for CEA in CEA_lst:
            CEA.check_ConstructElementAlignment_validity()
            if not CEA.valid:
                self.valid = False
                return False
            self.valid = True

    def check_all_ConstructElementAlignments_concatamer_validity(self):
        '''
        Checks validated CEAs for concatamers.
        '''
        long_CEAs = [self.index_CEAP.long_CEA, self.barcode_CEAP.long_CEA]
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
        
        self.orientation.append(self.index_CEAP.orientation)
        self.orientation.append(self.barcode_CEAP.orientation)
        self.orientation = list(chain(*self.orientation))
        self.valid = True
        return True

    def trim_ConstructElements_from_SeqRecord(self):
        '''
        Use the coordinates from the long CEAs to remove them from the SeqRecord.
        '''
        # psuedocode
        # get orientation of index_CEAP and barcode_CEAP.
        # if forward, split the SeqRecord into two separate things and then add them together.
            # include edge case where the long CEA is on the edge of a sequence
        # edit the seqrecord, leave everything else unchanged.

        if self.orientation == ['F','R']:
            FSpan = self.index_CEAP.long_CEA.FBoundary.span
            RSpan = self.barcode_CEAP.long_CEA.RBoundary.span

        elif self.orientation == ['R','F']:
            FSpan = self.barcode_CEAP.long_CEA.FBoundary.span
            RSpan = self.index_CEAP.long_CEA.RBoundary.span

        else:
            ValueError("Orientation is invalid for DCA")

        # flip RSpan around to compensate for it being taken from the reverse complement

        seq_len=np.int64(len(self.SeqRecord.seq))
        RSpan[0]=seq_len - RSpan[0]
        RSpan[1]=seq_len - RSpan[1]

        # swap the values with a tempvar
        temp_val=RSpan[0]
        RSpan[0] = RSpan[1]
        RSpan[1] = temp_val
        #print(f'FSpan:{FSpan}\tRSpan:{RSpan}')

        # now we use the cut sites we derived to actually cut the seq
        first_segment=self.SeqRecord[0:FSpan[0]]
        second_segment=self.SeqRecord[FSpan[1]:RSpan[0]]
        third_segment=self.SeqRecord[RSpan[1]:seq_len]

        trim_SeqRecord = first_segment + second_segment + third_segment
        self.SeqRecord = trim_SeqRecord
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
def init_aligner(open_gap_score=-.5, extend_gap_score=-.1):
    # Generic aligner we'll initialize once
    # We penalize opening gaps because our markers should theoretically be one group
    aligner=Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    return(aligner)

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

def main():

    seq_name='test_seq'
    input_seq='TTTTCTGTCCTGTACTTCGTTCAGTTAGGTATTGTTCAAGCAGAAGACAGAATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATTCTTTGGTTTATTTAAGATTAAACGCATACATGTTTGCAGAATTGCTGCCACAATAAGGTACTTTATCTTTTATACACAGGGAGCATTTACATTTTATACATATGAAAATATACCTCTAAGGACTTTTTTTTTTGCAACAAAACACAGCAGTTACCGACGCCTGATTCCCAGCTGGGGTAAGTCAGCTTGCAGACACTGTAGGAGCTGTGATGGTTGTAGCAGCTGAGATCTAGATGTACAGTCATGTCGAGTTTCCTTAAATATATGACACAAATGTACTACTCTCACACTCAGAGTGGCAACTTAGCACAACTATCTGCCAGCGCATGAGCATCTCTCAGTCCCAGAGAAAACCTTTATCCCTTACCTACACAGGTAATCTTTGAAACACTGACAGCACAACAACTAAACAAAATCTTGTATCATAAAGGCATAGAATTTAGGTTCTTTTTGATGAGAAAATCATCAAATACAGCTCTGAGCCAAAGCCCAGTGATGTTCCAGCCCTTTCTGCTGCCACCACCAGGCATGTCCCATGAAGGCCTGCAGAAGCTAGATAGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAGCAATACGT'
    index_full='CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATT'
    barcode_full='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACACCTTGCA'
    index='CGTGAT'
    barcode='ACACCT'

    #seq, DC, DCA=create_sample_data(seq_name, input_seq, index_full, index, barcode_full, barcode)

    #print(DC)

    #print(DCA)

if __name__=="__main__":
    main()