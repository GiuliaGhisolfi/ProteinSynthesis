import simpy
import random 
import itertools
from Bio.Seq import Seq
from src.protein_synthesis import EucaryotesCell
from src.variables import EucaryotesCellVariables
from src.resources.resource import EucaryotesCellResource
from src.utils import save_proteins_synthesized

LENGTH_AMIO_GROUP = 4 # length of amino acid group
LENGTH_CARBOXYL_GROUP = 5 # length of carboxyl group
RANDOM_SEED = 42
SIM_TIME = 1000
NUMBER_RESOURCES = 5

class ProteinSinthesisProcess:
    def __init__(self, dna_sequences_df, verbose=False):
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
        
        random.seed(RANDOM_SEED)
        self.env = simpy.Environment()
        self.resources = EucaryotesCellResource(self.env, capacity=NUMBER_RESOURCES) #TODO: resources: enzimi, basi, ATP, tRNA, aminoacidi
        self.env.process(self.setup_process())

        self.eucaryotes_cell = EucaryotesCell(environment=self.env, verbose=self.verbose)
        
        #if self.verbose: 
        print('Simulation environment initialized')
    
    def run(self, simulation_time=SIM_TIME):
        #if self.verbose: 
        print('Simulation started:')
        self.env.run(until=simulation_time)
        #if self.verbose: 
        print('End simulation')
    
    def setup_process(self):
        process_queue = []
        sequences_count = itertools.count()

        while True:
            dna_sequence = random.choice(self.dna_sequences) #TODO: Seq(random.choice(self.dna_sequences))
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
            #if self.verbose:
            print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} requesting to start synthesis')
            variables.request_start_process_time = self.env.now

            yield request # wait for a cell be able to accepts dna sequence

            #if self.verbose:
            print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} synthesize started')
            variables.start_process_time = self.env.now

            yield self.env.process(self.eucaryotes_cell.synthesize_protein(variables))
            variables.end_process_time = self.env.now

            # save the results
            self.save_proteins_synthesized(variables)

            #if self.verbose:
            print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} synthetis ended')
            self.resources.release(request)
        
    def save_proteins_synthesized(self, variables):
        self.dna_sequences_df = save_proteins_synthesized(
            variables.get_dna(),
            variables.get_mrna(),
            variables.get_proteins(),
            variables.get_extended_proteins_name(),
            variables.request_start_process_time,
            variables.start_process_time,
            variables.start_transcription_time,
            variables.start_translation_time,
            variables.end_translation_time,
            variables.end_process_time
            )
