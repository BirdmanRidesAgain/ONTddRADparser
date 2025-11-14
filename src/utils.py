__all__ = ["enzyme_lst"]
from Bio import Restriction

enzyme_lst=list(Restriction.__dict__)


#def get_valid_seq_subseq_aln_boundaries(seq_lst, construct_lst, subseq_name_lst, percent_match):
#    '''
#    Takes an seq_list, and two sets of unique barcodes (assumed to be fuzzy-matching), and a minimum match hit.
#    Returns a list of all seqs where 1 of each subseq was found, along with their boundaries.
#    '''
#    boundaries_lst = []
#
#    for col in construct_lst[subseq_name_lst]:
#        unique_subseq=construct_lst[col].unique()
#        seq_lst, seq_subseq_aln_boundaries_lst = filter_seqs_by_single_subseq_validity(seq_lst, col, unique_subseq, percent_match)
#        boundaries_lst.append(seq_subseq_aln_boundaries_lst)
#    print(len(boundaries_lst[0]),len(boundaries_lst[1]))
#    boundaries_lst=list(zip(boundaries_lst[0],boundaries_lst[1]))
#
#    # Now we need to do some kind of zip magic to make this look nice
#    #print(len(seq_lst), len(boundaries_lst))
#    seq_boundaries_lst=list(zip(seq_lst,boundaries_lst))
#    print(seq_boundaries_lst[0])
#    exit(0)
#    return seq_boundaries_lst
#
#
#def filter_seqs_by_single_subseq_validity(seq_lst, subseq_name, subseq_lst, percent_match):
#    '''
#    Takes a list of SeqRecord (`seq_lst`) objects, a list of Seqs (`subseq_lst`), and a float (`percent_match`)
#        - `subseqs` in this context are assumed to be from demux file
#        - Eg, `subseqs` are: 'full_idx', 'idx', 'full_barcode', 'barcode'
#        - 'full_idx' and 'full_barcode' require fuzzy matches; the other two must be exact
#
#    Function first aligns each combination of `seq` and `subseq` forward (`f`) and in reverse (`r`).
#    Filters `f` and `r` alignments separately, marking them as invalid/not found if the alignment is low-qual.
#    Four scenarios are possible:
#        1. Valid `subseq` found in `f` and `r`
#        2. Valid `subseq found in `f` only
#        3. Valid `subseq found in `r` only
#        4. No valid `subseq` found
#
#    Scenarios 2,3 are valid; all others are removed.
#    Finally, all valid SeqRecords and boundaries are returned as a tuple.
#    '''
#    print(f"Input seqs\t{len(seq_lst)}")
#    valid_seq_lst=[]
#    seq_subseq_aln_boundaries_lst=[]
#    dup_lst=[]
#    for i in seq_lst:
#        for j in subseq_lst:
#            if (i.name in dup_lst):
#                continue
#            else:
#                seq_subseq_aln_boundaries=get_aln_boundaries(i, subseq_name, j, percent_match)
#                if (XOR_aln_boundaries(seq_subseq_aln_boundaries[2], seq_subseq_aln_boundaries[3])):
#                    valid_seq_lst.append(i)
#                    seq_subseq_aln_boundaries_lst.append(seq_subseq_aln_boundaries)
#                    dup_lst.append(i.name)
#    print(f"nSeqs with valid substrings \t{len(valid_seq_lst)}")
#    return(valid_seq_lst, seq_subseq_aln_boundaries_lst)
#
#def XOR_aln_boundaries(f_idx, r_idx):
#    '''
#    Takes two barcode indices (`f` and `r`), which can be 'valid' or 'invalid' and returns a bool.
#    Evaluated `f` and `r` with XOR logic:
#        - `f` and `r` have alternate validity values, return False.
#        - Otherwise return True.
#    `f_idx` and `r_idx` should look like this:
#    ['f',[int,int]],['r',[int,int]]
#    '''
#    if (((f_idx[1]==[-1,-1]) & (r_idx[1]==[-1,-1])) | ((f_idx[1]!=[-1,-1]) & (r_idx[1]!=[-1,-1]))):
#        return False
#    return True
#
#def pair_lists_by_id(A, B):
#    # Create dictionaries mapping ID to element for A and B
#    dict_A = {item[0]: item for item in A}
#    dict_B = {item[0]: item for item in B}
#
#    # Find common IDs
#    common_ids = set(dict_A.keys()) & set(dict_B.keys())
#
#    # Build C: [ID, element_from_A, element_from_B]
#    C = [[id_, dict_A[id_], dict_B[id_]] for id_ in common_ids]
#    return C