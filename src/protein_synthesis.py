import simpy
import json
from src.transcription import Nucleus
from src.translation import Ribosome

DATA_PATH = 'data/'
CODONS_PATH = DATA_PATH + 'codons.json'
PEPTIDES_PATH = DATA_PATH + 'peptides.json'
PROMOTERS = ['TATAAAA', 'TATAAAT', 'TATATAA', 'TATATAT']
TERMINATORS = ['UAA', 'UAG', 'UGA']

class EucaryotesCell:
    def __init__(self):
        self.codons2aminoacids_dict = json.load(open(CODONS_PATH))
        self.aminoacids_dict = json.load(open(PEPTIDES_PATH))
        self.extron_list = self.codons2aminoacids_dict.keys()


    def synthesize_protein(self, dna):
        self.dna = dna # template strand (3' to 5' direction)

        # transcription
        self.nucleus = Nucleus(
            extron_sequences_list=self.extron_list,
            editing_sites_dict={}, #TODO
            promoters_sequence_list=PROMOTERS, # TATA box
            terminator_sequence=TERMINATORS
            )
        self.mrna = self.nucleus.transcript(self.dna)

        if self.mrna is None: # no promoter found in the DNA
            self.protein, self.protein_extended_name = None, None
        else:
            # translation
            self.ribosome = Ribosome(
                codons2aminoacids_dict=self.codons2aminoacids_dict, 
                aminoacids_dict=self.aminoacids_dict,
            )
            self.protein, self.protein_extended_name = self.ribosome.translate(self.mrna)

            #TODO
            # protein folding: ordered three-dimensional structure
            # protein degradation: if the protein is not correctly folded, it is degraded by the proteasome
    
    def get_dna(self):
        return self.dna
    
    def get_mrna(self):
        return self.mrna
    
    def get_protein(self):
        return self.protein
    
    def get_extended_protein_name(self):
        return self.protein_extended_name