from __future__ import annotations

__all__ = ["Boundary", "ConstructElement", "ConstructElementAlignmentPair", "DemuxConstruct", "DemuxConstructAlignment", "DemuxxedSample", "FastqFile", "init_aligner"]

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
    def __init__(self, type='undefined', orientation='undefined', start_idx=np.nan, end_idx=np.nan, buffer=0):
        self.type = type
        self.orientation = orientation
        self.start_idx = start_idx
        self.end_idx = end_idx
        self.buffer = buffer

    def __str__(self):
        str=f'''
        Type\tOrientation\tStart_idx\tEnd_index\tBuffer
        {self.type}\t{self.orientation}\t{self.start_idx}\t{self.end_idx}\t{self.buffer}
        '''
        return(str)

    def get_span(self):
        '''Calculates the span of the boundary, taking the buffer into account when present.
        If the buffer puts the span past the start of the seq (ie, if it's less than 0), defaults to 0.'''
        span_start_idx=self.start_idx - self.buffer
        if (span_start_idx < 0):
            span_start_idx=0
        span_end_idx = self.end_idx + self.buffer
        span_list=[span_start_idx, span_end_idx]
        return span_list

class ConstructElementAlignmentPair:
    '''
    Two complementary `ConstructElements` - one long, and one short. In the context of ONTddRADparser, they consist of:
        - `index` and `index_full`
        - `barcode` and `barcode_full`
    The short one is around 6-9 bp, and is allowed to exist more than once in the SeqRecord.
    The long one should only be found once - otherwise, it is a concatamer.
    '''


    def __init__(self, CEA_long: 'ConstructElementAlignment', CEA_short: 'ConstructElementAlignment'):
        self.valid = False
        self.CEA_long = CEA_long
        self.CEA_short = CEA_short
        self.orientation = '' # can be ['F','R','invalid']
        self.short_in_long = '' # can be T or F
    
    def __str__(self):
        str = f'''
        CEA_short\tCEA_long
        {self.CEA_short}\t{self.CEA_long}
        '''
        return(str)   

    def check_ConstructElementAlignmentPair_validity(self):
        '''
        Checks if the CEPA is valid, based on self.orientation and self.short_in_long, and updates the self.valid tag.
        '''
        if (self.get_paired_orientation()):
            if (self.check_short_in_long()):
                self.valid = True
                return True
        self.valid = False
        return False

    def check_short_in_long(self):
        '''
        Performs position checking to ensure that CEA_short is inside CEA_long.
        '''
        if self.orientation == '':
            self.valid = self.get_paired_orientation()

        # now we need to check that the short element is found inside of the long element
        if self.orientation == 'F':
            short_span=self.CEA_short.FBoundary.get_span()
            long_span=self.CEA_long.FBoundary.get_span()
            if ((short_span[0] < long_span[0]) or (short_span[1] > long_span[1])):
                self.short_in_long = False
                return False

        else: # now we need to check that the short element is found inside of the long element
            short_span=self.CEA_short.RBoundary.get_span()
            long_span=self.CEA_long.RBoundary.get_span()
            if ((short_span[0] < long_span[0]) or (short_span[1] > long_span[1])):
                self.short_in_long = False
                return False
        self.short_in_long = True
        return True

    def get_paired_orientation(self):
        '''
        Determines whether the pair is `F`, `R`, or `invalid`, and updates `self.orientation`.
        '''
        CEA_long_orientation_F = (not np.isnan(self.CEA_long.FBoundary.start_idx) and not np.isnan(self.CEA_long.FBoundary.end_idx))
        CEA_long_orientation_R = (not np.isnan(self.CEA_long.RBoundary.start_idx) and not np.isnan(self.CEA_long.RBoundary.end_idx))

        # We expect exactly one aligned boundary in the long element
        if CEA_long_orientation_F == CEA_long_orientation_F:
            # Either none or both aligned → invalid pairing
            self.orientation = 'invalid'
            return False

        # Make checks for alignment and orientation of long+short elements        
        if CEA_long_orientation_F:
            # Long uses FBoundary; require short to also use FBoundary with same orientation
            CEA_short_orientation_F = (not np.isnan(self.CEA_short.FBoundary.start_idx)and not np.isnan(self.CEA_short.FBoundary.end_idx))  
            if not CEA_short_orientation_F:
                self.orientation = 'invalid'
                return False

        # do the same thing as above, except on the RBoundary
        else:
            CEA_short_orientation_R = (not np.isnan(self.CEA_short.RBoundary.start_idx) and not np.isnan(self.CEA_short.RBoundary.end_idx))
            if not CEA_short_orientation_R:
                self.orientation = 'invalid'
                return False

        # Assuming all other checks are passed, return the orientation of both of the values:
        if (CEA_long_orientation_F):
            self.orientation = 'F'
            return True
        else:
            self.orientation = 'R'
            return True

class ConstructElement:
    '''
    Individual elements of a larger DNA construct, implemented as a Seq object with metadata.
    In the context of ONTddRADparse, they correspond to the index, barcode, index_full and barcode_full elements.
    We include the 'type' of construct, (long, short, etc.), the alignment specificity (`aln_percent`), and the raw sequence.
    '''
    def __init__(self, seq, type, aln_percent: 'float', buffer):
        self.seq = Seq(seq)
        self.type = type
        self.aln_percent = aln_percent
        self.buffer = buffer

    def __str__(self):
        str = f'''
        seq\ttype\taln_percent\tbuffer
        {self.type}\t{self.aln_percent}\t{self.seq}\t{self.buffer}
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
        self.FBoundary = Boundary()
        self.RBoundary = Boundary()

    def __str__(self):
        str=f'''
        SeqRecord\tConstructElement_type\taln_valid\tFBoundary_start_idx\tFBoundary_end_idx\tRBoundary_start_idx\tRBoundary_end_idx
        {self.SeqRecord.id}\t{self.ConstructElement.type}\t{self.valid}\t{self.FBoundary.start_idx}\t{self.FBoundary.end_idx}\t{self.RBoundary.start_idx}\t{self.RBoundary.end_idx}
        '''
        return(str)

    def align_ConstructElement(self, orientation):
        '''Aligns the subset to the seq in either the 5->3 ('f') or 3->5 ('r') direction.'''
        aln = align_target(self.SeqRecord.seq, self.ConstructElement.seq, self.aligner, orientation, self.ConstructElement.aln_percent)
        if (orientation == 'f' or orientation == 'F'):
            self.FBoundary=Boundary(self.ConstructElement.type, orientation, aln[0], aln[1], self.ConstructElement.buffer)
        elif (orientation == 'r' or orientation == 'R'):    
            self.RBoundary=Boundary(self.ConstructElement.type, orientation, aln[0], aln[1], self.ConstructElement.buffer)
        else:
            raise ValueError(f"Ambiguous alignment orientation.\n\tAccepted values: 'f,F', 'r,R'.\n\tActual value: {orientation}")

    def check_ConstructElementAlignment_validity(self):
        '''
        A ConstructElementAlignment is valid if:
            - exactly ONE of the boundaries (FBoundary or RBoundary) is aligned, OR
            - BOTH boundaries are aligned AND the ConstructElement type is 'short'.

        A boundary is considered "aligned" when:
            - that boundary has numeric start_idx and end_idx
            - the other boundary has np.nan for both indices
        '''
        # Boundary is considered "aligned" only if BOTH indices are numeric
        aligned_in_F_orientation = (not np.isnan(self.FBoundary.start_idx)) and (not np.isnan(self.FBoundary.end_idx))
        aligned_in_R_orientation = (not np.isnan(self.RBoundary.start_idx)) and (not np.isnan(self.RBoundary.end_idx))

        # Base rule: valid when exactly one of aligned_in_F_orientation / aligned_in_R_orientation is True (logical XOR)
        xor_valid = bool(aligned_in_F_orientation ^ aligned_in_R_orientation)

        # Additional rule: for short constructs, allow both boundaries aligned
        short_both_aligned_valid = (self.ConstructElement.type == 'short') and aligned_in_F_orientation and aligned_in_R_orientation

        self.valid = xor_valid or short_both_aligned_valid    

class DemuxConstruct:
    '''
    Represents the construct sequence data we add to our sequences to demux them.
    Consists of a sample ID and four ConstructElement objects: two 6-9bp ones which must be exact matches, and two longer ones where error is allowed.
    Canonically, should be read in from your demux file.
    '''
    def __init__(self, sample_id, index_full: 'ConstructElement', index:'ConstructElement', barcode_full:'ConstructElement', barcode: 'ConstructElement', fuzzy_aln_percent, exact_aln_percent, buffer):
        exact_match_buffer=0
        self.sample_id = sample_id
        self.index_full = ConstructElement(index_full, 'long', fuzzy_aln_percent, buffer)
        self.index = ConstructElement(index, 'short', exact_aln_percent, exact_match_buffer)
        self.barcode_full = ConstructElement(barcode_full, 'long', fuzzy_aln_percent, buffer)
        self.barcode = ConstructElement(barcode, 'short', exact_aln_percent, exact_match_buffer)

    def __str__(self):
        str = f'''Construct object information is:
        Sample id: {self.sample_id}
        Full index: {self.index_full}
        Short index: {self.index}
        Full barcode: {self.barcode_full}
        Short barcode: {self.barcode}'''
        return(str)

class DemuxConstructAlignment:
    '''
    Represents the alignments between a SecRecord object and a DemuxConstruct.
    Contains the SeqRecord and DemuxConstruct.
    Also contains two ConstructElementPairAlignments, (idx and barcode).
    Finally, contains a 'valid' and 'reason' tag.

    '''
    invalidity_reason_lst = [] # list appended to whenever a check fails

    def __init__(self, SeqRecord: 'SeqRecord', DemuxConstruct: 'DemuxConstruct', aligner):
        self.valid = False     # we initialize validity as false until proven otherwise
        self.SeqRecord = SeqRecord
        self.DemuxConstruct = DemuxConstruct

        index_full_aln=ConstructElementAlignment(SeqRecord, DemuxConstruct.index_full, aligner)
        index_aln=ConstructElementAlignment(SeqRecord, DemuxConstruct.index, aligner)
        barcode_full_aln=ConstructElementAlignment(SeqRecord, DemuxConstruct.barcode_full, aligner)
        barcode_aln=ConstructElementAlignment(SeqRecord, DemuxConstruct.barcode, aligner)

        self.index_pair_aln = ConstructElementAlignmentPair(index_aln, index_full_aln)
        self.barcode_pair_aln = ConstructElementAlignmentPair(barcode_aln, barcode_full_aln)

    def __str__(self):

        str = f'''
        SeqRecord\tDemuxConstruct_sample_id\tDemuxConstructAlignment_valid\tindex_full_valid\tindex_valid\tbarcode_full_valid\tbarcode_valid
        {self.SeqRecord.id}\t{self.DemuxConstruct.sample_id}\t{self.valid}\t{self.index_full_aln.valid}\t{self.index_aln.valid}\t{self.barcode_full_aln.valid}\t{self.barcode_aln.valid}
        '''
        return(str)

    def align_all_ConstructElements(self):
        CEAs = [self.index_pair_aln.CEA_long, self.index_pair_aln.CEA_short, self.barcode_pair_aln.CEA_long, self.index_pair_aln.CEA_short]
        orientations = ['F', 'R']
        for CEA in CEAs:
            for orientation in orientations:
                CEA.align_ConstructElement(orientation)

    def check_all_ConstructElementAlignments_validity(self):
        '''
        Wrapper around ConstructElementAlignment.check_ConstructElementAlignment_validity().
        Ports the logic into the main script so we can weed out SeqQbjects that fail the first line of validation.
        '''
        CEAs = [self.index_pair_aln.CEA_long, self.index_pair_aln.CEA_short, self.barcode_pair_aln.CEA_long, self.index_pair_aln.CEA_short]
        for CEA in CEAs:
            CEA.check_ConstructElementAlignment_validity()
            if not CEA.valid:
                self.valid = False
                return False
            self.valid = True

    def check_all_ConstructElementAlignmentPairs_validity(self):
        '''
        Wrapper around ConstructElementPairAlignment.check_ConstructElementPairAlignment_validity().
        Assumes that all ConstructElements are valid.
        '''
        CEAPs = [self.index_pair_aln, self.barcode_pair_aln]
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
        if self.index_pair_aln.orientation == '':
            self.index_pair_aln.get_paired_orientation()
        if self.barcode_pair_aln.orientation == '':
            self.barcode_pair_aln.get_paired_orientation()

        if self.index_pair_aln.orientation == self.barcode_pair_aln.orientation:
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
def align_target(seq, subseq, aligner, orientation, aln_percent):
    '''
    Helper function for `check_seq_for_full_index` and others.
    Takes a Bio.Record object and a Bio.Seq object, and aligns them.
    Returns a tuple of the start and end indices of 'subseq' to 'seq'.
    Returns a list of two indices which can be formatted to a Boundary object.
    
    '''
    aligner.target_end_gap_score = aligner.query_end_gap_score = 0.0
    min_aln_score=len(subseq)*aln_percent

    if (orientation=='F'):
        aln=aligner.align(seq, subseq)
    elif (orientation=='R'):
        aln=aligner.align(seq.reverse_complement(),subseq)
    else:
        raise ValueError(f"Ambiguous alignment orientation.\n\tAccepted values: 'f,F', 'r,R'.\n\tActual value: {orientation}")
    
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
    alignment.barcode_boundaries=[Boundary('barcode','F' , 0, -1),Boundary('barcode','R' , 0, -1)]

    # testing method for setting attributes
    #print(alignment)


    barcode_element=ConstructElement('barcode', exact_aln, barcode)
    elementalignment=ConstructElementAlignment(seqrec1, barcode_element, aligner)
    elementalignment.align_ConstructElement('F')
    elementalignment.align_ConstructElement('R')
    #print(elementalignment)

if __name__=="__main__":
    main()