LENGTH_CODON = 3
LENGTH_POLY_A_TAIL = 5
MRNA_DECODED_ERROR_RATE = 1e-5 # 1 mistake every 10,000 amino acids

class Ribosome:
    def __init__(self, codons2aminoacids_dict, aminoacids_dict):
        self.codons2aminoacids_dict = codons2aminoacids_dict
        self.aminoacids_dict = aminoacids_dict

    def translate(self, mrna_sequence):
        #mrna_sequence = self.activation(mrna_sequence)
        mrna_sequence = self.initialization(mrna_sequence)
        protein = self.elongation(mrna_sequence)
        
        return protein
    
    def activation(self): #TODO
        # required energy from adenosine triphosphate (ATP)
        pass

    def initialization(self, mrna_sequence):
        subunit = '' # init
        subunit = mrna_sequence[:LENGTH_CODON]
        i = LENGTH_CODON

        # find start codon
        while subunit not in self.codons2aminoacids_dict.keys():
            subunit = subunit[1:]
            subunit = subunit + mrna_sequence[i]
            i += 1
        
        return mrna_sequence[i:]

    def elongation(self, mrna_sequence):
        aminoacid = ''
        aminoacids_chain = ''
        i = 0
        end_mrna_length = LENGTH_CODON + LENGTH_POLY_A_TAIL

        while not self.termination(aminoacid) and i<(len(mrna_sequence)-end_mrna_length):
            codon = mrna_sequence[i:i+LENGTH_CODON]
            aminoacid = self.codons2aminoacids_dict[codon]
            aminoacids_chain = aminoacids_chain + self.aminoacids_dict[aminoacid]
            i += LENGTH_CODON
    
        return aminoacids_chain

    def termination(self, aminoacid):
        return aminoacid == 'Stop'