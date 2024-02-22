import simpy
import random 
import itertools
from Bio.Seq import Seq
from src.protein_synthesis import EucaryotesCell
from src.variables import DNASequenceVariables

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
            'number_of_proteins_synthesized', 'protein_synthesized']
        for col in columns:
            self.dna_sequences_df[col] = None
        
        # initialize the simulation environment
        self.dna_sequences = self.dna_sequences_df['sequence'].values
        self.available =  {row['sequence']: True if row['protein_synthesized']==None else False
            for _, row in self.dna_sequences_df.iterrows()}
        
        random.seed(RANDOM_SEED)
        self.env = simpy.Environment()
        self.resources = simpy.Resource(self.env, capacity=NUMBER_RESOURCES) #TODO: resources: enzimi, basi, ATP, tRNA, aminoacidi
        self.env.process(self.setup_process())

        self.eucaryotes_cell = EucaryotesCell(environment=self.env, verbose=self.verbose)
        
        if self.verbose: print('Simulation environment initialized \t')
    
    def run(self, simulation_time=SIM_TIME):
        if self.verbose: print('Simulation started: \n')
        self.env.run(until=simulation_time)
    
    def setup_process(self):
        process_queue = []
        sequences_count = itertools.count()

        while True:
            dna_sequence = random.choice(self.dna_sequences) #TODO: Seq(random.choice(self.dna_sequences))
            if self.available[dna_sequence]:
                variables = DNASequenceVariables()
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
            print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} requesting resources')
            yield request # wait for a cell be able to accepts dna sequence
            print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} got resource')
            
            if self.verbose:
                print(f'Time {self.env.now:.4f}: Protein synthesis started for sequence {variables.sequence_count}')
            yield self.env.process(self.eucaryotes_cell.synthesize_protein(variables))

            self.save_proteins_synthesized(
                variables.get_dna(),
                variables.get_mrna(),
                variables.get_proteins(),
                variables.get_extended_proteins_name()
            )

            print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} synthetis ended')
            self.resources.release(request)
            
    def save_proteins_synthesized(self, dna_sequence, mrna_sequences, polypeptides_chain, polypeptides_chain_ext):
        # TODO: gestire errori (dna_sequence non trovata nel dataframe)
        row_index = self.dna_sequences_df[self.dna_sequences_df['sequence'] == dna_sequence].index[0]

        results = {
            'ID': self.dna_sequences_df.iloc[row_index]['ID'],
            'sequence': dna_sequence,
            'category': self.dna_sequences_df.iloc[row_index]['category'],
            'mrna_sequences': mrna_sequences,
            'polypeptides_chains': polypeptides_chain,
            'polypeptides_chains_ext': polypeptides_chain_ext,
            'number_of_proteins_synthesized': len(mrna_sequences) if mrna_sequences else 0,
            'protein_synthesized': True if mrna_sequences else False
        }
        self.dna_sequences_df.iloc[row_index] = results

        if self.verbose:
            print(f'Protein synthesized: {len(mrna_sequences) if mrna_sequences else 0} \n')