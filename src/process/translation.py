from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio import BiopythonWarning
import warnings
import random
from src.resources.resource import EucaryotesCellResource

START_CODON = 'AUG' # start codon
LENGTH_CODON = 3 # number of nucleotides that code for an amino acid
LENGTH_METHYL_CAP = 8 # length of 5'-methyl cap
LENGTH_POLY_A_TAIL = 5 # length of poly-A tail
AMINO_GROUP = 'NH2-' # amino group
CARBOXYL_GROUP = '-COOH' # carboxyl group
ELONGATION_TIME = 5e-2 # seconds to add each amino acid
MRNA_DECODED_ERROR_RATE = 1e-4 # 1 mistake every 10.000 amino acids

class Ribosome:
    def __init__(self, environment, number_ribosomes, nucleotides, amminoacids):
        # ribonucleoprotein complex in the cytoplasm
        self.env = environment
        self.ribosomes = EucaryotesCellResource(self.env, capacity=number_ribosomes)
        self.nucleotides = nucleotides
        self.amminoacids = amminoacids

    def translate(self, mrna_sequence, variables): # protein synthesis
        with self.ribosomes.request() as request:
            yield request # wait for a ribosome to be available
            
            polypeptides_chain, polypeptides_chain_ext = yield self.env.process(
                self.translation_process(mrna_sequence, variables.poly_adenine_tail_len))
            
        return polypeptides_chain, polypeptides_chain_ext
    
    def translation_process(self, mrna_sequence, poly_adenine_tail_len):
        self.ribosomes.available() # register the time when the resource is available

        # start translation
        mrna_sequence = self.degradation_cap_tail(mrna_sequence)
        initial_mrna_sequence = mrna_sequence

        #mrna_sequence = self.activation(mrna_sequence)
        mrna_sequence = self.initialization(mrna_sequence)
        polypeptides_chain, polypeptides_chain_ext = yield self.env.process(
            self.elongation(mrna_sequence))
        
        self.mrna_degredation(initial_mrna_sequence, poly_adenine_tail_len)
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
        start_codon = START_CODON # start codon
        start_codon_position = str(mrna_sequence).find(start_codon)

        return mrna_sequence[start_codon_position+LENGTH_CODON:]
    
    def elongation(self, mrna_sequence):
        mrna_sequence = Seq(mrna_sequence)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            polypeptides_chain_seq = mrna_sequence.translate(stop_symbol='', to_stop=True)
        
        # error in the translation process
        polypeptides_chain_list = list(polypeptides_chain_seq)
        polypeptides_chain = ''
        for amminoacid in polypeptides_chain_list:
            polypeptides_chain = polypeptides_chain + (random.choice([a for a in list(self.amminoacids) 
                 if a != amminoacid]) if random.random() <= MRNA_DECODED_ERROR_RATE else amminoacid)

        if len(polypeptides_chain) > 0:
            yield self.env.timeout(ELONGATION_TIME * len(polypeptides_chain)) # 0.05 seconds to add each amino acid
            polypeptides_chain_ext = seq3(polypeptides_chain)

            polypeptides_chain = AMINO_GROUP + polypeptides_chain + CARBOXYL_GROUP
            polypeptides_chain_ext = AMINO_GROUP + polypeptides_chain_ext + CARBOXYL_GROUP

            return str(polypeptides_chain), str(polypeptides_chain_ext)
        else:
            return None, None
        
    
    def mrna_degredation(self, mrna_sequence, poly_adenine_tail_len):
        # enzima: ribonuclease
        [self.release_nucleotide(nucleotide) for nucleotide in mrna_sequence]
        self.release_nucleotide('A', poly_adenine_tail_len)
        #TODO: yield self.env.timeout() degradation time
    
    def release_nucleotide(self, base, amount=1):
        self.nucleotides.release(base, amount)