import simpy
import random 
import itertools
from src.process.protein_synthesis import EucaryotesCell
from src.variables.variables import EucaryotesCellVariables
from src.resources.resource import EucaryotesCellResource
from src.utils import save_proteins_synthesized

LENGTH_AMIO_GROUP = 4 # length of amino acid group
LENGTH_CARBOXYL_GROUP = 5 # length of carboxyl group
RESULTS_FOLDER = 'results/'

SIM_TIME = 1000
NUMBER_RESOURCES = 5
NUMBER_RNA_POLYMERASES = 3
NUMBER_RIBOSOMES = 2
RANDOM_SEED = 42

class ProteinSinthesisProcess:
    def __init__(self, dna_sequences_df, number_resources=NUMBER_RESOURCES,
            number_rna_polymerases=NUMBER_RNA_POLYMERASES, number_ribosomes=NUMBER_RIBOSOMES,
            random_seed=RANDOM_SEED, verbose=False):
        self.dna_sequences_df = dna_sequences_df
        self.verbose = verbose

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
        self.resources = EucaryotesCellResource(self.env, capacity=number_resources) 
        #TODO: resources: enzimi, basi, ATP, tRNA, aminoacidi
        self.env.process(self.setup_process())

        self.eucaryotes_cell = EucaryotesCell(environment=self.env, number_rna_polymerases=number_rna_polymerases,
            number_ribosomes=number_ribosomes, random_seed=random_seed, verbose=self.verbose)
        
        print('Simulation environment initialized')
    
    def __str__(self):
        return (f'Protein Sinthesis Process:\n'
            f'{len(self.dna_sequences)} dna sequences to synthesize,\n'
            f'{self.resources.capacity} resources available,\n'
            f'{self.eucaryotes_cell.nucleus.rna_polymerase.capacity} RNA polymerases,\n'
            f'{self.eucaryotes_cell.ribosome.ribosomes.capacity} ribosomes.')
    
    def run(self, simulation_time=SIM_TIME):
        print('Simulation started')
        self.env.run(until=simulation_time)

        proteins_number = self.dna_sequences_df[self.dna_sequences_df[
            'protein_synthesized'].notna()]['number_of_proteins_synthesized'].sum()
        print(f'End simulation: {proteins_number} proteins synthesized.')
    
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
                yield self.env.timeout(random.random()) # time between start of protein synthesis

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
            self.save_protein_synthesized(variables)

            if self.verbose:
                print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} synthetis ended')
            yield self.resources.release(request)
        
    def save_protein_synthesized(self, variables):
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
            end_process_time=variables.end_process_time
            )
        
    def save_process(self):
        # save dataframe
        df_to_save = self.dna_sequences_df[self.dna_sequences_df['protein_synthesized'].notna()]
        df_to_save.to_csv(RESULTS_FOLDER+'results.csv')

        # save resources history 
        self.resources.save_history(RESULTS_FOLDER+'resources_history.csv')
        self.eucaryotes_cell.nucleus.rna_polymerase.save_history(
            RESULTS_FOLDER+'rna_polymerase_history.csv')
        self.eucaryotes_cell.ribosome.ribosomes.save_history(
            RESULTS_FOLDER+'ribosome_history.csv')
        
        print('Process saved.')
