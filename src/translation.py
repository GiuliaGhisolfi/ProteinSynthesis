LENGTH_CODON = 3 # number of nucleotides that code for an amino acid
LENGTH_METHYL_CAP = 8 # length of 5'-methyl cap
LENGTH_POLY_A_TAIL = 5 # length of poly-A tail
MRNA_DECODED_ERROR_RATE = 1e-4 # 1 mistake every 10.000 amino acids

class Ribosome:
    def __init__(self, codons2aminoacids_dict, aminoacids_dict):
        # ribonucleoprotein complex in the cytoplasm
        self.codons2aminoacids_dict = codons2aminoacids_dict
        self.aminoacids_dict = aminoacids_dict

    def translate(self, mrna_sequence): # protein synthesis
        mrna_sequence = self.degradation(mrna_sequence)
        #mrna_sequence = self.activation(mrna_sequence)
        mrna_sequence = self.initialization(mrna_sequence)
        polypeptides_chain, polypeptides_chain_ext = self.elongation(mrna_sequence)
        #TODO: mechanism to correct transcription errors
        
        return polypeptides_chain, polypeptides_chain_ext
    
    def degradation(self, mrna_sequence):
        # degradation of the 5' cap and poly-A tail, enzime: exonuclease
        return mrna_sequence[:-LENGTH_POLY_A_TAIL]
    
    def activation(self): #TODO
        # required energy from adenosine triphosphate (ATP)
        pass

    def initialization(self, mrna_sequence):
        subunit = '' # init
        subunit = mrna_sequence[:LENGTH_CODON]
        i = LENGTH_CODON

        # find start codon: AUG (methionine)
        while subunit != 'AUG' and i<len(mrna_sequence):
            subunit = subunit[1:]
            subunit = subunit + mrna_sequence[i]
            i += 1

        return mrna_sequence[i-LENGTH_CODON:]

    def elongation(self, mrna_sequence):
        # enzima: aminoacyl-tRNA synthetases
        aminoacid = ''
        polypeptides_chain = 'H2N-' # amino group
        polypeptides_chain_ext = 'H2N-'

        i = LENGTH_CODON # not translation of the start codon
        end_mrna_length = LENGTH_CODON + LENGTH_POLY_A_TAIL

        while not self.termination(aminoacid) and i<(len(mrna_sequence)-end_mrna_length):
            codon = mrna_sequence[i:i+LENGTH_CODON]
            aminoacid = self.codons2aminoacids_dict[codon]
            polypeptides_chain = polypeptides_chain + self.aminoacids_dict[aminoacid]
            polypeptides_chain_ext = polypeptides_chain_ext + aminoacid + '-'
            i += LENGTH_CODON
        
        # add the carboxyl group to the polypeptide chain
        polypeptides_chain = polypeptides_chain + '-COOH'
        polypeptides_chain_ext = polypeptides_chain_ext.rstrip('-Stop-') + '-COOH'
    
        return polypeptides_chain, polypeptides_chain_ext

    def termination(self, aminoacid):
        return aminoacid == 'Stop'