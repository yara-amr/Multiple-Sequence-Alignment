#Merna Hesham Mahmoud
#Haidy Kamal Makram
#Yara Amr Ahmed
#Maram Nasser
#Maram Hatem
#Yasmine Mohamed El-Gazar


from PairwiseAlignment import PairwiseAlignment
from MyAlign import MyAlign
from MySeq import MySeq
from SubstMatrix import SubstMatrix
class MultipleAlignment():
 def _init_( self , sequences, alignmentsequence):
    self.sequences = sequences
    self.alignment_pars = alignmentsequence

def add_sequence_alignment(self, alignment, seq):
    res = []

    for i in range(len(alignment.listseqs) + 1):
        res.append("")
    cons = MySeq(alignment.consensus(), alignment.al_type)
    self.alignpars.needleman_Wunsch(cons, seq)
    alignment2 = self.alignment_pars.recover_align()
    orig = 0
    for i in range(len(alignment2)):

       if alignment2[0, i] == '-':
        for k in range(len(alignment.listseqs)):
            res[k] += "−"
        else:
            for k in range(len(alignment.listseqs)):
                res[k] += alignment[k, orig]
        orig += 1
        res[len(alignment.listseqs)] = alignment2.listseqs[1]
        return MyAlign(res, alignment.al_type)

def alignment_consensus(self):
     self.alignpars.needleman_Wunsch(self.seqs[0], self.seqs[1])

     res = self.alignpars.recover_align()
     for i in range(2, len(self.seqs)):
        res = self.add_sequence_alignment(res, self.seqs[i])
        return res

def test():
        first_sequence = MySeq("ATAGC")

        second_sequence = MySeq("AACC")
        third_sequence = MySeq("ATGAC")
        sm = SubstMatrix()
        sm.create_submat(1,-1 , "ACGT")
        aseq = PairwiseAlignment(sm,-1)
        ma = MultipleAlignment([first_sequence,second_sequence,third_sequence], aseq)
        al = ma.align_consensus()
        print(al)
        if __name__ == "_main_":
         test()
