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

        # start transcription
        if self.verbose:
            print(f'Time {self.env.now:.4f}: Transcription started')

        # split the DNA sequence in the promoter regions
        dna_sequences_to_transcript_list, promoters_count = self.detect_promoter_process(dna_sequence)

        # continue with transcription and translation if promoters are found
        if dna_sequences_to_transcript_list is None:
            self.mrna_list = None
            self.proteins, self.proteins_extended_name = None, None
        else:
            self.init_queue()
            sequences_count = itertools.count()
            
            # transcription and translation each promoter region
            for dna_sequence_to_transcript in dna_sequences_to_transcript_list:
                seq_count = next(sequences_count)

                if self.verbose:
                    print(f'Time {self.env.now:.4f}: Transcription started for mRNA sequence {seq_count}')
                
                yield self.env.process(self.transcription_and_translation_process
                    (dna_sequence=dna_sequence_to_transcript, seq_count=seq_count))

                if self.verbose:
                    print(f'Time {self.env.now:.4f}: Transcription ended for mRNA sequence {seq_count}')
                    if seq_count == promoters_count:
                        print(f'Time {self.env.now:.4f}: Transcription ended')

            if self.verbose:
                print(f'Time {self.env.now:.4f}: Translation ended') 
    
    def detect_promoter_process(self, dna_sequence):
        # detect promoter
        dna_sequences_to_transcript_list = self.nucleus.find_promoter(dna_sequence)
        
        promoters_count = len(dna_sequences_to_transcript_list) if dna_sequences_to_transcript_list is not None else 0
        if self.verbose:
            print(f'Promoters found: {promoters_count}')

        return dna_sequences_to_transcript_list, promoters_count
    
    def init_queue(self):
        self.transcription_queue = []
        self.translation_queue = []

        self.mrna_list = []
        self.proteins = []
        self.proteins_extended_name = []
    
    def transcription_and_translation_process(self, dna_sequence, seq_count):
        # transcription process
        transcription_process = self.env.process(self.nucleus.transcript(dna_sequence))
        self.transcription_queue.append(transcription_process)
        
        yield transcription_process
        mrna = transcription_process.value
        self.mrna_list.append(mrna)

        # wait for all the transcription process to be completed
        while self.transcription_queue:
            yield self.transcription_queue.pop(0)

        # translation
        if self.verbose:
            if seq_count == 0:
                print(f'Time {self.env.now:.4f}: Translation started')
            print(f'Time {self.env.now:.4f}: Translation started for mRNA sequence {seq_count}')
        
        self.env.process(self.translation_process(self.mrna_list[seq_count]))

        if self.verbose:
            print(f'Time {self.env.now:.4f}: Translation ended for mRNA sequence {seq_count}')
    
    def translation_process(self, mrna):
        # translation process
        translation_process = self.env.process(self.ribosome.translate(mrna))
        self.translation_queue.append(translation_process)
        
        yield translation_process
        protein, protein_extended_name = translation_process.value
        self.proteins.append(protein) # polypeptides chain
        self.proteins_extended_name.append(protein_extended_name)

        while self.translation_queue:
            yield self.translation_queue.pop(0)
    
    def get_dna(self):
        return self.dna
    
    def get_mrna(self):
        return self.mrna_list
    
    def get_proteins(self):
        return self.proteins
    
    def get_extended_proteins_name(self):
        return self.proteins_extended_name