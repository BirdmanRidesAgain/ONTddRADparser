from Bio.Seq import Seq

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

def main():
    # tests the class
    id='R10N00251'	
    index_full='CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAATT'
    index='CGTGAT'
    barcode_full='AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTACACCTTGCA'
    barcode='ACACCT'

    sample=Construct(id, index_full, index, barcode_full, barcode)
    print(sample.sample_id)
    print(sample.index_full)
    print(sample.index)
    print(sample.barcode_full)
    print(sample.barcode)

if __name__=="__main__":
    main()