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
        
    def synthesize_protein(self, variables):
        # start transcription
        if self.verbose:
            print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} start transcription process')

        # split the DNA sequence by promoter regions
        self.detect_promoter_process(variables)

        # continue with transcription and translation if promoters are found
        if variables.dna_sequences_to_transcript_list is None:
            variables.mrna_sequences_list = None
            variables.proteins_list, variables.proteins_extended_name_list = None, None
        else:
            sequences_count = itertools.count()
            variables.mrna_sequences_list = []
            variables.proteins_list, variables.proteins_extended_name_list = [], []
            
            # transcription and translation each promoter region
            for _ in variables.dna_sequences_to_transcript_list:
                seq_count = next(sequences_count)

                if self.verbose:
                    print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} '
                        f'(mRNA sequence {seq_count}) start transcription process')

                yield self.env.process(self.transcription_and_translation_process
                    (variables, seq_count=seq_count))                    

            if self.verbose:
                print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} '
                    f'end transcription and translation process')
    
    def detect_promoter_process(self, variables):
        # detect promoter
        variables.dna_sequences_to_transcript_list = self.nucleus.find_promoter(variables.dna_sequence)
        
        if variables.dna_sequences_to_transcript_list is not None:
            variables.promoters_count = len(variables.dna_sequences_to_transcript_list)  
        else: 
            variables.promoters_count = 0

        if self.verbose:
            print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} contains '
                f'{variables.promoters_count} promoters')
    
    def transcription_and_translation_process(self, variables, seq_count):
        # transcription process
        transcription_process = self.env.process(self.nucleus.transcript(variables.dna_sequence))
        variables.transcription_queue.append(transcription_process)
        
        yield transcription_process
        mrna = transcription_process.value
        variables.mrna_sequences_list.append(mrna)

        # wait for all the transcription process to be completed
        while variables.transcription_queue:
            variables.transcription_queue.pop(0)

        # translation
        if self.verbose:
            if seq_count == 0:
                print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} start translation process')
            print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} (mRNA sequence {seq_count}) '
                f'start translation process')
        
        yield self.env.process(self.translation_process(variables, mrna))

        if self.verbose:
            print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} (mRNA sequence {seq_count}) '
                f'end translation process')
    
    def translation_process(self, variables, mrna):
        # translation process
        translation_process = self.env.process(self.ribosome.translate(mrna))
        variables.translation_queue.append(translation_process)
        
        yield translation_process
        protein, protein_extended_name = translation_process.value
        variables.proteins_list.append(protein) # polypeptides chain
        variables.proteins_extended_name_list.append(protein_extended_name)

        while variables.translation_queue:
            variables.translation_queue.pop(0)