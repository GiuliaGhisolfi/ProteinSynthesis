from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio import BiopythonWarning
import warnings
from src.resources.resource import EucaryotesCellResource

LENGTH_CODON = 3 # number of nucleotides that code for an amino acid
LENGTH_METHYL_CAP = 8 # length of 5'-methyl cap
LENGTH_POLY_A_TAIL = 5 # length of poly-A tail
MRNA_DECODED_ERROR_RATE = 1e-4 # 1 mistake every 10.000 amino acids
NUMBER_RIBOSOMES = 2
AMINO_GROUP = 'NH2-' # amino group
CARBOXYL_GROUP = '-COOH' # carboxyl group

class Ribosome:
    def __init__(self, environment, codons2aminoacids_dict, aminoacids_dict):
        # ribonucleoprotein complex in the cytoplasm
        self.env = environment

        self.codons2aminoacids_dict = codons2aminoacids_dict
        self.aminoacids_dict = aminoacids_dict

        #self.ribosomes = simpy.Resource(self.env, capacity=NUMBER_RIBOSOMES)
        self.ribosomes = EucaryotesCellResource(self.env, capacity=NUMBER_RIBOSOMES)

    def translate(self, mrna_sequence): # protein synthesis
        with self.ribosomes.request() as request:
            yield request # wait for a ribosome to be available
            polypeptides_chain, polypeptides_chain_ext = yield self.env.process(
                self.translation_process(mrna_sequence))
            
            self.ribosomes.release(request) # release the ribosome
            
        return polypeptides_chain, polypeptides_chain_ext
    
    def translation_process(self, mrna_sequence):
        mrna_sequence = self.degradation_cap_tail(mrna_sequence)
        #mrna_sequence = self.activation(mrna_sequence)
        mrna_sequence = self.initialization(mrna_sequence)
        polypeptides_chain, polypeptides_chain_ext = yield self.env.process(
            self.elongation(mrna_sequence))
        #TODO: mechanism to correct transcription errors + translation times

        return polypeptides_chain, polypeptides_chain_ext
    
    def degradation_cap_tail(self, mrna_sequence):
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
        mrna_sequence = Seq(mrna_sequence)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            polypeptides_chain = mrna_sequence.translate(stop_symbol='', to_stop=True)

        if len(polypeptides_chain) > 0:
            yield self.env.timeout(0.05* len(polypeptides_chain)) # 0.05 seconds to add each amino acid
            polypeptides_chain_ext = seq3(polypeptides_chain)

            polypeptides_chain = AMINO_GROUP + polypeptides_chain + CARBOXYL_GROUP
            polypeptides_chain_ext = AMINO_GROUP + polypeptides_chain_ext + CARBOXYL_GROUP
        
        return str(polypeptides_chain), str(polypeptides_chain_ext)
    
    def mrna_degredation(self, mrna_sequence):
        # enzima: ribonuclease
        pass