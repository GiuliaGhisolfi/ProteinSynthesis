import simpy
import json
from src.transcription import Nucleus
from src.translation import Ribosome

DATA_PATH = 'data/'
CODONS_PATH = DATA_PATH + 'codons.json'
PEPTIDES_PATH = DATA_PATH + 'peptides.json'
PROMOTERS = ['TATAAAA', 'TATAAAT', 'TATATAA', 'TATATAT']
TERMINATORS = ['UAA', 'UAG', 'UGA'] # TODO: valutare se cambiare stop codons

class EucaryotesCell:
    def __init__(self, dna):
        self.dna = dna # template strand (3' to 5' direction)
        self.codons_dict = json.load(open(CODONS_PATH))
        self.peptides_dict = json.load(open(PEPTIDES_PATH))
        self.extron_list = self.codons_dict.keys()

        env = simpy.Environment
        #pipeline = BroadcastPipe(env)

        self.protein = self.synthesize_protein()

    def synthesize_protein(self):
        self.nucleus = Nucleus(
            extron_sequences_list=self.extron_list,
            editing_sites_dict={}, #TODO
            promoters_sequence_list=PROMOTERS, # TATA box
            terminator_sequence=TERMINATORS
            )
        self.mrna = self.nucleus.transcript(self.dna)
        #TODO: Implement Translation
    
    def get_dna(self):
        return self.dna
    
    def get_mrna(self):
        return self.mrna
    
    def get_protein(self):
        return self.protein