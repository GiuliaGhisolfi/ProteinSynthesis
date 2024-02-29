import simpy
import random 
import itertools
import os
from src.process.protein_synthesis import EucaryotesCell
from src.variables.variables import EucaryotesCellVariables
from src.resources.resource import EucaryotesCellResource
from src.utils.utils import save_proteins_synthesized, post_processing_results

LENGTH_AMIO_GROUP = 4 # length of amino acid group
LENGTH_CARBOXYL_GROUP = 5 # length of carboxyl group
RESULTS_FOLDER = 'results/'

SIM_TIME = 1000
NUMBER_RESOURCES = 200
NUMBER_RNA_POLYMERASES = 3
NUMBER_RIBOSOMES = 2
NUMBER_RNA_TRANSFER = 2000
URACIL_INITIAL_AMOUNT = 5000
ADENINE_INITIAL_AMOUNT = 5000
GUANINE_INITIAL_AMOUNT = 5000
CYTOSINE_INITIAL_AMOUNT = 5000
RANDOM_SEED = None

class ProteinSinthesisProcess:
    def __init__(self, 
            dna_sequences_df, 
            number_resources=NUMBER_RESOURCES,
            number_rna_polymerases=NUMBER_RNA_POLYMERASES, 
            number_ribosomes=NUMBER_RIBOSOMES,
            number_rna_transfers_per_codon=NUMBER_RNA_TRANSFER, 
            uracil_initial_amount=URACIL_INITIAL_AMOUNT, 
            adenine_initial_amount=ADENINE_INITIAL_AMOUNT, 
            guanine_initial_amount=GUANINE_INITIAL_AMOUNT, 
            cytosine_initial_amount=CYTOSINE_INITIAL_AMOUNT,
            random_seed=RANDOM_SEED, verbose=False):
        self.dna_sequences_df = dna_sequences_df
        self.verbose = verbose

        # nucleotides
        self.uracil_initial_amount = uracil_initial_amount
        self.adenine_initial_amount = adenine_initial_amount
        self.guanine_initial_amount = guanine_initial_amount
        self.cytosine_initial_amount = cytosine_initial_amount

        # add columns to store the results
        columns = ['mrna_sequences', 'polypeptides_chains', 'polypeptides_chains_ext',
            'number_of_proteins_synthesized', 'protein_synthesized', 'request_start_process_time',
            'start_process_time', 'start_transcription_time', 'start_translation_time',
            'end_translation_time', 'end_process_time']
        for col in columns:
            self.dna_sequences_df[col] = None
        
        # initialize the simulation environment
        self.dna_sequences = self.dna_sequences_df['sequence'].values
        self.available =  {row['sequence']: True if row['protein_synthesized']==None else False
            for _, row in self.dna_sequences_df.iterrows()}
        
        random.seed(random_seed)
        self.env = simpy.Environment()
        self.resources = EucaryotesCellResource(
            self.env, capacity=number_resources, save_history=False) 
        #TODO: resources: enzimi, ATP
        self.env.process(self.setup_process())

        self.eucaryotes_cell = EucaryotesCell(
            environment=self.env, 
            number_rna_polymerases=number_rna_polymerases,
            number_ribosomes=number_ribosomes, 
            number_rna_transfers_per_codon=number_rna_transfers_per_codon,
            uracil_initial_amount=uracil_initial_amount, 
            adenine_initial_amount=adenine_initial_amount, 
            guanine_initial_amount=guanine_initial_amount,
            cytosine_initial_amount=cytosine_initial_amount, 
            random_seed=random_seed, 
            verbose=self.verbose
            )
        
        print('Simulation environment initialized, time unit: 0.0001 second.')
    
    def __str__(self):
        return (f'Protein Sinthesis Process:\n'
            f'{len(self.dna_sequences)} dna sequences to synthesize,\n'
            f'{self.resources.capacity} resources available,\n'
            f'{self.eucaryotes_cell.nucleus.rna_polymerase.capacity} RNA polymerases,\n'
            f'{self.eucaryotes_cell.ribosome.ribosomes.capacity} ribosomes,\n'
            f'{self.eucaryotes_cell.ribosome.rna_transfer.capacity} RNA transfer,\n'
            f'{self.uracil_initial_amount} uracil bases,\n'
            f'{self.adenine_initial_amount} adenine bases,\n'
            f'{self.guanine_initial_amount} guanine bases,\n'
            f'{self.cytosine_initial_amount} cytosine bases.')
    
    def run(self, simulation_time=SIM_TIME):
        print('Simulation started')
        self.env.run(until=simulation_time)

        # save simulation results
        proteins_number = self.dna_sequences_df[self.dna_sequences_df[
            'protein_synthesized'].notna()]['number_of_proteins_synthesized'].sum()
        dna_sequences_processed_number = self.dna_sequences_df[
            self.dna_sequences_df['protein_synthesized'].notna()].shape[0]
        
        print(f'End simulation: {proteins_number} proteins synthesized from '
            f'{dna_sequences_processed_number} DNA sequences.')
    
    def setup_process(self):
        process_queue = []
        sequences_count = itertools.count()

        while True:
            dna_sequence = random.choice(self.dna_sequences)
            if self.available[dna_sequence]:
                variables = EucaryotesCellVariables()
                variables.dna_sequence = dna_sequence
                variables.sequence_count = next(sequences_count)

                process_queue.append(self.env.process(self.process(variables)))
                
                self.available[dna_sequence] = False
                yield self.env.timeout(round(random.random()*100, ndigits=4)) # time between start of protein synthesis

                while process_queue: # wait for all the protein synthesis to be completed
                    process_queue.pop(0)
        
    def process(self, variables):
        # Synthesize dna sequences while the simulation is running       
        with self.resources.request() as request:
            if self.verbose:
                print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} requesting to start synthesis')
            variables.request_start_process_time = self.env.now

            yield request # wait for a cell be able to accepts dna sequence

            if self.verbose:
                print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} synthesize started')
            variables.start_process_time = self.env.now

            yield self.env.process(self.eucaryotes_cell.synthesize_protein(variables))
            variables.end_process_time = self.env.now

            # save the results
            self.save_proteins_synthesized_in_df(variables)

            if self.verbose:
                print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} synthetis ended')
        
    def save_proteins_synthesized_in_df(self, variables):
        self.dna_sequences_df = save_proteins_synthesized(
            dna_sequences_df=self.dna_sequences_df, 
            dna_sequence=variables.get_dna(),
            mrna_sequences=variables.get_mrna(),
            polypeptides_chain=variables.get_proteins(),
            polypeptides_chain_ext=variables.get_extended_proteins_name(),
            request_start_process_time=variables.request_start_process_time,
            start_process_time=variables.start_process_time,
            start_transcription_time=variables.start_transcription_time,
            start_translation_time=variables.start_translation_time,
            end_translation_time=variables.end_translation_time,
            end_process_time=variables.end_process_time,
            promoters_box=variables.promoters_box,
            )
        
    def save_process(self, folder_test_name=''):
        # create folder to save the results
        if folder_test_name != '':
            if not os.path.exists(RESULTS_FOLDER+folder_test_name):
                os.mkdir(RESULTS_FOLDER+folder_test_name)
            folder_test_name = folder_test_name + '/'

        # save dataframe        
        df_to_save = self.dna_sequences_df[self.dna_sequences_df['protein_synthesized'].notna()]
        df_to_save = df_to_save.apply(post_processing_results, axis=1)
        df_to_save.to_csv(RESULTS_FOLDER+folder_test_name+'results.csv')

        # save resources history 
        self.eucaryotes_cell.nucleus.rna_polymerase.save_history(
            RESULTS_FOLDER+folder_test_name+'rna_polymerase_history.json')
        
        self.eucaryotes_cell.ribosome.ribosomes.save_history(
            RESULTS_FOLDER+folder_test_name+'ribosome_history.json')
        
        if not os.path.exists(RESULTS_FOLDER+folder_test_name+'nucleotides'):
            os.mkdir(RESULTS_FOLDER+folder_test_name+'nucleotides')
        self.eucaryotes_cell.nucleotides.save_history(
            RESULTS_FOLDER+folder_test_name+'nucleotides/'+'nucleotides_history.json')

        if not os.path.exists(RESULTS_FOLDER+folder_test_name+'rna_transfer'):
            os.mkdir(RESULTS_FOLDER+folder_test_name+'rna_transfer')
        self.eucaryotes_cell.ribosome.rna_transfer.save_history(
            RESULTS_FOLDER+folder_test_name+'rna_transfer/'+'rna_transfer_history.json')
        
        print('Process saved.')
