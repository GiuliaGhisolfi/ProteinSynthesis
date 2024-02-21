import json
import itertools
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

    def synthesize_protein(self, dna_sequence):
        self.dna = dna_sequence # template strand (3' to 5' direction)
        sequences_count = itertools.count()

        # start transcription
        if self.verbose:
            print(f'Time {self.env.now}: Transcription started')

        # detect promoter
        dna_sequences_to_transcript_list = self.nucleus.find_promoter(dna_sequence)
        
        promoters_count = len(dna_sequences_to_transcript_list) if dna_sequences_to_transcript_list is not None else 0
        if self.verbose:
            print(f'Promoters found: {promoters_count}')

        if dna_sequences_to_transcript_list is None:
            self.mrna_list = None
            self.proteins, self.proteins_extended_name = None, None
        else:
            self.mrna_list = []
            self.proteins = []
            self.proteins_extended_name = []

            while len(self.mrna_list) < promoters_count:
                # FIXME: DNA sequence non vengono processate consequenzialmente, ma aspettano che termini la precedente
                dna_sequence_to_transcript = dna_sequences_to_transcript_list[len(self.mrna_list)]
                seq_count = next(sequences_count)

                # transcription process
                if self.verbose:
                    print(f'Time {self.env.now}: Transcription started for mRNA sequence {seq_count}')

                mrna = yield self.env.process(self.nucleus.transcript(dna_sequence_to_transcript))
                self.mrna_list.append(mrna) # mature mRNA

                if self.verbose:
                    print(f'Time {self.env.now}: Transcription ended for mRNA sequence {seq_count}')
                    if len(self.mrna_list) == promoters_count:
                        print(f'Time {self.env.now}: Transcription ended') 

                # translation
                if self.verbose:
                    if len(self.mrna_list) == 0:
                        print(f'Time {self.env.now}: Translation started')
                    print(f'Time {self.env.now}: Translation started for mRNA sequence {seq_count}')
                
                protein, protein_extended_name = yield self.env.process(
                    self.ribosome.translate(mrna))
                self.proteins.append(protein) # polypeptides chain
                self.proteins_extended_name.append(protein_extended_name)

                if self.verbose:
                    print(f'Time {self.env.now}: Translation ended for mRNA sequence {seq_count}')

            if self.verbose:
                print(f'Time {self.env.now}: Translation ended')
    
    def get_dna(self):
        return self.dna
    
    def get_mrna(self):
        return self.mrna_list
    
    def get_proteins(self):
        return self.proteins
    
    def get_extended_proteins_name(self):
        return self.proteins_extended_name