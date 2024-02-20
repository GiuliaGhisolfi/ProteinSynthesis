import random
from Bio.Seq import Seq
from Bio import SeqUtils

LENGTH_CODON = 3 # number of nucleotides that code for an amino acid
LENGTH_METHYL_CAP = 8 # length of 5'-methyl cap
LENGTH_POLY_A_TAIL = 5 # length of poly-A tail
MRNA_DECODED_ERROR_RATE = 1e-4 # 1 mistake every 10.000 amino acids

class Ribosome:
    def __init__(self, environment, codons2aminoacids_dict, aminoacids_dict):
        # ribonucleoprotein complex in the cytoplasm
        self.env = environment

        self.codons2aminoacids_dict = codons2aminoacids_dict
        self.aminoacids_dict = aminoacids_dict

    def translate(self, mrna_sequences_list): # protein synthesis
        polypeptides_chain_list = []
        polypeptides_chain_ext_list = []

        for mrna_sequence in mrna_sequences_list: 
            #TODO: process mrna sequences in parallel

            mrna_sequence = self.degradation(mrna_sequence)
            #mrna_sequence = self.activation(mrna_sequence)
            mrna_sequence = self.initialization(mrna_sequence)
            polypeptides_chain, polypeptides_chain_ext = self.elongation(mrna_sequence)
            #TODO: mechanism to correct transcription errors

            polypeptides_chain_list.append(polypeptides_chain)
            polypeptides_chain_ext_list.append(polypeptides_chain_ext)
        
        return polypeptides_chain_list, polypeptides_chain_ext_list
    
    def degradation(self, mrna_sequence):
        # degradation of the 5' cap and poly-A tail, enzime: exonuclease
        return mrna_sequence[LENGTH_METHYL_CAP:-LENGTH_POLY_A_TAIL]
    
    def activation(self): #TODO
        # required energy from adenosine triphosphate (ATP) to activate tRNA
        # reaction catalyzed by the enzyme aminoacyl-tRNA synthetase
        # reaction: amino acid + ATP + tRNA -> aminoacyl-tRNA + AMP + PPi
        pass

    def initialization(self, mrna_sequence):
        start_codon = 'AUG' # start codon
        start_codon_position = str(mrna_sequence).find(start_codon)

        return mrna_sequence[start_codon_position:]

    def elongation(self, mrna_sequence):
        # enzima: aminoacyl-tRNA synthetases
        aminoacid = ''
        polypeptides_chain = 'NH2-' # amino group
        polypeptides_chain_ext = 'NH2-'

        i = LENGTH_CODON # not translation of the start codon
        end_mrna_length = LENGTH_CODON + LENGTH_POLY_A_TAIL

        while not self.termination(aminoacid) and i<(len(mrna_sequence)-end_mrna_length):
            codon = mrna_sequence[i:i+LENGTH_CODON]
            aminoacid = self.codons2aminoacids_dict[codon]
            polypeptides_chain = polypeptides_chain + self.aminoacids_dict[aminoacid]
            polypeptides_chain_ext = polypeptides_chain_ext + aminoacid + '-'
            i += LENGTH_CODON
            #FIXME: yield self.env.timeout(0.05) # 0.05 seconds to add each amino acid
        
        # add the carboxyl group to the polypeptide chain
        polypeptides_chain = polypeptides_chain + '-COOH'
        polypeptides_chain_ext = polypeptides_chain_ext.rstrip('-Stop-') + '-COOH'
    
    
        return polypeptides_chain, polypeptides_chain_ext

    def termination(self, aminoacid):
        return aminoacid == 'Stop'
        # TODO: add nucleotide degradation