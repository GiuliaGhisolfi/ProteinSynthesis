import json
import itertools
import random
from src.process.transcription import Nucleus
from src.process.translation import Ribosome
from src.resources.nucleotides import Nucleotides

DATA_PATH = 'data/'
CODONS_PATH = DATA_PATH + 'codons.json'

class EucaryotesCell:
    def __init__(self, environment, number_rna_polymerases,number_ribosomes, number_rna_transfers_per_codon, 
            uracil_initial_amount, adenine_initial_amount, guanine_initial_amount,
            cytosine_initial_amount, random_seed, verbose=False):
        self.env = environment
        self.verbose = verbose
        random.seed(random_seed)

        self.extron_list = json.load(open(CODONS_PATH)).keys()
        self.amminoacids = json.load(open(CODONS_PATH)).values()

        self.nucleotides = Nucleotides(
            environment=self.env,
            uracil_initial_amount=uracil_initial_amount,
            adenine_initial_amount=adenine_initial_amount,
            guanine_initial_amount=guanine_initial_amount,
            cytosine_initial_amount=cytosine_initial_amount,
            random_seed=random_seed
            )

        self.nucleus = Nucleus(
            environment=self.env,
            extron_sequences_list=self.extron_list,
            editing_sites_dict={},
            number_rna_polymerases=number_rna_polymerases,
            nucleotides = self.nucleotides,
            random_seed=random_seed
            )
        
        self.ribosome = Ribosome(
            environment=self.env,
            number_ribosomes=number_ribosomes,
            number_rna_transfers_per_codon=number_rna_transfers_per_codon,
            codons_list=self.extron_list,
            nucleotides = self.nucleotides,
            amminoacids = self.amminoacids,
            random_seed=random_seed
            )
        
    def synthesize_protein(self, variables):
        # start transcription
        if self.verbose:
            print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} start transcription process')
        variables.found_promoter_time = self.env.now


        # split the DNA sequence by promoter regions
        self.detect_promoter_process(variables)

        # continue with transcription and translation if promoters are found
        if variables.dna_sequences_to_transcript_list is not None:
            sequences_count = itertools.count()
            variables.init_transcription_translation_var() # init variables
            
            # transcription and translation each promoter region
            for dna_sequence in variables.dna_sequences_to_transcript_list:
                seq_count = next(sequences_count)

                if self.verbose:
                    print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} '
                        f'(mRNA sequence {seq_count}) start transcription process')
                variables.start_transcription_time.append(self.env.now)

                yield self.env.process(self.transcription_and_translation_process
                    (dna_sequence, variables, seq_count=seq_count))                    

            if self.verbose:
                print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} '
                    f'end transcription and translation process')
    
    def detect_promoter_process(self, variables):
        # detect promoter
        variables.dna_sequences_to_transcript_list = self.nucleus.find_promoter(
            variables.dna_sequence, variables)
        
        if variables.dna_sequences_to_transcript_list is not None:
            variables.promoters_count = len(variables.dna_sequences_to_transcript_list)  
        else: 
            variables.promoters_count = 0

        if self.verbose:
            print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} contains '
                f'{variables.promoters_count} promoters')
    
    def transcription_and_translation_process(self, dna_sequence, variables, seq_count):
        # init list to store simpy processes for finding complement base
        variables.complement_base_queue_dict[seq_count] = []

        # transcription process
        transcription_process = self.env.process(
            self.nucleus.transcript(dna_sequence, variables, seq_count))
        variables.transcription_queue.append(transcription_process)
        
        yield transcription_process
        mrna = transcription_process.value
        variables.mrna_sequences_list[seq_count] = mrna

        # wait for all the transcription process to be completed
        while variables.transcription_queue:
            variables.transcription_queue.pop(0)

        # translation
        if self.verbose:
            if seq_count == 0:
                print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} start translation process')
            print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} (mRNA sequence {seq_count}) '
                f'start translation process')
        variables.start_translation_time.append(self.env.now)
    
        yield self.env.process(self.translation_process(variables, mrna, seq_count))

        if self.verbose:
            print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} (mRNA sequence {seq_count}) '
                f'end translation process')
        variables.end_translation_time.append(self.env.now)
    
    def translation_process(self, variables, mrna, seq_count):
        # translation process
        translation_process = self.env.process(self.ribosome.translate(mrna, variables, seq_count))
        variables.translation_queue.append(translation_process)
        
        yield translation_process
        protein, protein_extended_name, mrna_degradated = translation_process.value

        if variables.proteins_sintetized[seq_count] == 1:
            variables.proteins_list[seq_count] = protein # polypeptides chain
            variables.proteins_extended_name_list[seq_count] = protein_extended_name

        while variables.translation_queue:
            variables.translation_queue.pop(0)
        
        if not mrna_degradated:
            yield self.env.timeout(round(random.random()*10, ndigits=4)) # time to find the next ribosome
            yield self.env.process(self.translation_process(variables, mrna, seq_count))