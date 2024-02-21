import json
from src.transcription import Nucleus
from src.translation import Ribosome

DATA_PATH = 'data/'
CODONS_PATH = DATA_PATH + 'codons.json'
PEPTIDES_PATH = DATA_PATH + 'peptides.json'

class EucaryotesCell:
    def __init__(self, environment, verbose=False):
        self.env = environment
        self.verbose = verbose

        self.codons2aminoacids_dict = json.load(open(CODONS_PATH))
        self.aminoacids_dict = json.load(open(PEPTIDES_PATH))
        self.extron_list = self.codons2aminoacids_dict.keys()

        self.nucleus = Nucleus(
            environment=self.env,
            extron_sequences_list=self.extron_list,
            editing_sites_dict={}, #TODO
            )
        
        self.ribosome = Ribosome(
            environment=self.env,
            codons2aminoacids_dict=self.codons2aminoacids_dict, 
            aminoacids_dict=self.aminoacids_dict,
            )

    def synthesize_protein(self, dna):
        self.dna = dna # template strand (3' to 5' direction)

        # transcription
        if self.verbose:
            print(f'Time {self.env.now}: Transcription started')
        transcript_generator = yield self.env.process(self.nucleus.transcript(self.dna))
        #yield from transcript_generator # wait for the transcription to end
        if self.verbose:
            print(f'Time {self.env.now}: Transcription ended') 

        if transcript_generator is not None:
            self.mrna_list = []
            for result in transcript_generator: # list of process
                self.mrna_list.append(result.value)
            
            # translation
            if self.verbose: 
                print(f'mRNA synthesized: {len(self.mrna_list) if self.mrna_list is not None else 0}')
                print(f'Time {self.env.now}: Translation started')
            self.proteins, self.proteins_extended_name = self.ribosome.translate(self.mrna_list)
            #FIXME: i want this as a process
            if self.verbose: 
                print(f'Time {self.env.now}: Translation ended')
        else:
            self.mrna_list = None
            self.proteins, self.proteins_extended_name = None, None

            #TODO
            # protein folding: ordered three-dimensional structure
            # protein degradation: if the protein is not correctly folded, it is degraded by the proteasome
    
    def get_dna(self):
        return self.dna
    
    def get_mrna(self):
        return self.mrna_list
    
    def get_proteins(self):
        return self.proteins
    
    def get_extended_proteins_name(self):
        return self.proteins_extended_name