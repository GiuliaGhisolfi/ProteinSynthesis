import simpy
import json
from src.transcription import Nucleus
from src.translation import Cytoplasm

DATA_PATH = 'data/'
CODONS_PATH = DATA_PATH + 'codons.json'
PEPTIDES_PATH = DATA_PATH + 'peptides.json'

class EucaryotesCell:
    def __init__(self, dna):
        self.dna = dna # template strand (3' to 5' direction)
        self.codons = json.load(open(CODONS_PATH))
        self.peptides = json.load(open(PEPTIDES_PATH))

        env = simpy.Environment
        #pipeline = BroadcastPipe(env)

        self.protein = self.synthesize_protein()

    def synthesize_protein(self):
        self.nucleus = Nucleus(
            intron_sequences_list=[], # TODO
            editing_sites_dict={}, #TODO
            promoters_sequence_list=['TATAAAA', 'TATAAAT', 'TATATAA', 'TATATAT'], # TATA box
            terminator_sequence=['UAA', 'UAG', 'UGA'] # TODO: valutare se cambiare stop codons
            )
        self.mrna = self.nucleus.transcript(self.dna)
        #TODO: Implement Translation
    
    def get_dna(self):
        return self.dna
    
    def get_mrna(self):
        return self.mrna
    
    def get_protein(self):
        return self.protein