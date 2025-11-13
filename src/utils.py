__all__ = ["get_aln_boundaries","filter_aln_by_score","XOR_aln_boundaries","filter_seqs_by_single_subseq_validity","filter_seqs_by_multiple_subseq_validity","get_valid_seq_subseq_aln_boundaries"]
from Bio import Restriction
from Bio import Align

enzyme_lst=list(Restriction.__dict__)

def get_valid_seq_subseq_aln_boundaries(seq_lst, demux_df, subseq_name_lst, percent_match):
    '''
    Takes an seq_list, and two sets of unique barcodes (assumed to be fuzzy-matching), and a minimum match hit.
    Returns a list of all seqs where 1 of each subseq was found, along with their boundaries.
    '''
    boundaries_lst = []

    for col in demux_df[subseq_name_lst]:
        unique_subseq=demux_df[col].unique()
        seq_lst, seq_subseq_aln_boundaries_lst = filter_seqs_by_single_subseq_validity(seq_lst, col, unique_subseq, percent_match)
        boundaries_lst.append(seq_subseq_aln_boundaries_lst)
    print(len(boundaries_lst[0]),len(boundaries_lst[1]))
    boundaries_lst=list(zip(boundaries_lst[0],boundaries_lst[1]))

    # Now we need to do some kind of zip magic to make this look nice
    #print(len(seq_lst), len(boundaries_lst))
    seq_boundaries_lst=list(zip(seq_lst,boundaries_lst))
    print(seq_boundaries_lst[0])
    exit(0)
    return seq_boundaries_lst


def filter_seqs_by_single_subseq_validity(seq_lst, subseq_name, subseq_lst, percent_match):
    '''
    Takes a list of SeqRecord (`seq_lst`) objects, a list of Seqs (`subseq_lst`), and a float (`percent_match`)
        - `subseqs` in this context are assumed to be from demux file
        - Eg, `subseqs` are: 'full_idx', 'idx', 'full_barcode', 'barcode'
        - 'full_idx' and 'full_barcode' require fuzzy matches; the other two must be exact

    Function first aligns each combination of `seq` and `subseq` forward (`f`) and in reverse (`r`).
    Filters `f` and `r` alignments separately, marking them as invalid/not found if the alignment is low-qual.
    Four scenarios are possible:
        1. Valid `subseq` found in `f` and `r`
        2. Valid `subseq found in `f` only
        3. Valid `subseq found in `r` only
        4. No valid `subseq` found

    Scenarios 2,3 are valid; all others are removed.
    Finally, all valid SeqRecords and boundaries are returned as a tuple.
    '''
    print(f"Input seqs\t{len(seq_lst)}")
    valid_seq_lst=[]
    seq_subseq_aln_boundaries_lst=[]
    dup_lst=[]
    for i in seq_lst:
        for j in subseq_lst:
            if (i.name in dup_lst):
                continue
            else:
                seq_subseq_aln_boundaries=get_aln_boundaries(i, subseq_name, j, percent_match)
                if (XOR_aln_boundaries(seq_subseq_aln_boundaries[2], seq_subseq_aln_boundaries[3])):
                    valid_seq_lst.append(i)
                    seq_subseq_aln_boundaries_lst.append(seq_subseq_aln_boundaries)
                    dup_lst.append(i.name)
    print(f"nSeqs with valid substrings \t{len(valid_seq_lst)}")
    return(valid_seq_lst, seq_subseq_aln_boundaries_lst)

def XOR_aln_boundaries(f_idx, r_idx):
    '''
    Takes two barcode indices (`f` and `r`), which can be 'valid' or 'invalid' and returns a bool.
    Evaluated `f` and `r` with XOR logic:
        - `f` and `r` have alternate validity values, return False.
        - Otherwise return True.
    `f_idx` and `r_idx` should look like this:
    ['f',[int,int]],['r',[int,int]]
    '''
    if (((f_idx[1]==[-1,-1]) & (r_idx[1]==[-1,-1])) | ((f_idx[1]!=[-1,-1]) & (r_idx[1]!=[-1,-1]))):
        return False
    return True

def get_aln_boundaries(seq, subseq_name, subseq, percent_match=1):
    '''
    Takes a seq and a subseq, aligns them and returns a list of the alignments.
    Function searches both the forward and reverse complement.
    Returns an output list of alignment boundaries structured like: 
        [seq_name,idx,[f,(aln_start, aln_end)],[r,(aln_start, aln_end)]]
    '''
    # for every sequence/index combination, align and find the indices of all high-qual alignments
    max_aln_score=len(subseq)
    aln_f=align_target(seq, subseq, 'f')
    aln_r=align_target(seq, subseq, 'r')

    subseq_boundary_lst_f=filter_aln_by_score(aln_f, max_aln_score, percent_match)
    subseq_boundary_lst_r=filter_aln_by_score(aln_r, max_aln_score, percent_match)

    subseq_loc = ([seq.name, subseq_name, subseq.__str__(),['f',subseq_boundary_lst_f],['r',subseq_boundary_lst_r]])
    return(subseq_loc)

def filter_aln_by_score(aln, max_aln_score, match_percent=1):
    '''
    Takes an alignment, a max score and a minimum percent of that max score needed to pass.
    Returns a list - either the alignment indices or the invalid indices [-1,-1].
    '''
    min_aln_score=max_aln_score*match_percent
    if (aln.score < min_aln_score):
        return([-1,-1])
    else:
        # the alignment list is always 2 elements long.
        # the first element contains the indices corresponding to the sequence, which we want
        # taking the min and max gets us the boundary of where the index aligned plus whatever gaps were opened
        seq_aln_boundaries=aln[0].aligned[0].flatten()
        return(min(seq_aln_boundaries), max(seq_aln_boundaries)) # first element of the alignment

def align_target(seq, subseq, orientation):
    '''
    Helper function for `check_seq_for_full_index` and others.
    Takes a Bio.Record object and a Bio.Seq object, and aligns them.
    Returns a tuple of the start and end indices of 'subseq' to 'seq'.
    '''
    # We penalize opening gaps because our markers should theoretically be one group
    aligner=Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = aligner.query_end_gap_score = 0.0

    if (orientation=='f'):
        alignment=aligner.align(seq, subseq)
    elif (orientation=='r'):
        alignment=aligner.align(seq.reverse_complement(),subseq)
    else:
        raise ValueError(f"Ambiguous alignment orientation.\n\tAccepted values: 'f', 'r'.\n\tActual value: {orientation}")
    return(alignment)



def pair_lists_by_id(A, B):
    # Create dictionaries mapping ID to element for A and B
    dict_A = {item[0]: item for item in A}
    dict_B = {item[0]: item for item in B}

    # Find common IDs
    common_ids = set(dict_A.keys()) & set(dict_B.keys())

    # Build C: [ID, element_from_A, element_from_B]
    C = [[id_, dict_A[id_], dict_B[id_]] for id_ in common_ids]
    return C

# Example usage:
# A = [[1, 'foo'], [2, 'bar'], [3, 'baz']]
# B = [[2, 'qux'], [3, 'quux'], [4, 'corge']]
# C = pair_lists_by_id(A, B)
# print(C)






# DEPRECATED CODE
def filter_seqs_by_multiple_subseq_validity(seq_lst, boundary_lst):
    '''
    CURRENTLY DEPRECATED.

    Takes a list of SeqRecord (`seq_lst`) objects and alignment boundary indices for 2 subseqs.
    Uses XOR logic via `XOR_aln_boundaries` to remove seqs where both subseqs are in cis.
    All valid SeqRecords and boundaries are returned as a tuple.
    '''
    print(f"Input seqs\t{len(seq_lst)}")
    valid_seq_lst=[]
    seq_subseq_aln_boundaries_lst=[]

    for i,j in zip(seq_lst, boundary_lst):
        # We only test the seq1_f and seq2_r pair
        # Prior checks guarantee that one of them will always exist, and the only condition is that one, but not both, will exist
        # 
        subseq_1_boundaries_f = j[0][2] # element 2 is f, element 3 is r
        subseq_1_boundaries_r = j[0][3]

        subseq_2_boundaries_f = j[1][2] # element 2 is f, element 3 is r
        subseq_2_boundaries_r = j[1][3]
        #if ((XOR_aln_boundaries(subseq_1_boundaries_f, subseq_2_boundaries_r)) | (XOR_aln_boundaries(subseq_2_boundaries_f, subseq_1_boundaries_r))):
        if (XOR_aln_boundaries(subseq_2_boundaries_f, subseq_1_boundaries_r)):
            valid_seq_lst.append(i)
            seq_subseq_aln_boundaries_lst.append(j)
    print(f"nSeqs with valid substrings \t{len(valid_seq_lst)}")
    return(valid_seq_lst, seq_subseq_aln_boundaries_lst)

class DemuxConstruct:
    '''
    Container to represent a DNA sequence and the four regions of its demux construct.
    Contains one SeqRecord object and four DemuxSubseq objects.
    DemuxSubseq objects represent `index_full`	`index`	`barcode_full`	`barcode`
    '''
    
    def __init__(self, seq=None, index_full=None, index=None, barcode_full=None, barcode=None):
        self.seq=seq
        self.index_full=index_full
        self.index=index
        self.barcode_full=barcode_full
        self.barcode=barcode

class DemuxSubseq:
    '''
    Container for a demux subsequence.
    Represents one of `index_full` `index` `barcode_full` or `barcode`, where these are `types`.
    Initialized with one Seq object, a `seq_type`, and a `map_type`.
    Has slots for `seq_type`, `map_type`, `orientation`, `start_idx`, and `end_idx`.
    '''
    def __init__(self, seq=None, seq_type=None, map_type=None, orientation=None, start_idx=None, end_idx=None):
        self.seq=seq
        self.seq_type=seq_type
        self.map_type=map_type
        self.orientation=orientation
        self.start_idx=start_idx
        self.end_idx=end_idx
