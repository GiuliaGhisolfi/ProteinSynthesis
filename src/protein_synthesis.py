import json
from src.transcription import Transcription
from src.translation import Translation

DATA_PATH = 'data/'
CODONS_PATH = DATA_PATH + 'codons.json'
PEPTIDES_PATH = DATA_PATH + 'peptides.json'

class ProteinSynthesis:
    def __init__(self, dna):
        self.dna = dna
        self.codons = json.load(open(CODONS_PATH))
        self.peptides = json.load(open(PEPTIDES_PATH))

        self.protein = self.synthesize_protein()

    def synthesize_protein(self):
        pass