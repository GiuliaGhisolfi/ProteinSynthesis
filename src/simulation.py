import simpy
import random 
from src.protein_synthesis import EucaryotesCell

LENGTH_AMIO_GROUP = 4 # length of amino acid group
LENGTH_CARBOXYL_GROUP = 5 # length of carboxyl group
RANDOM_SEED = 42
SIM_TIME = 1000


class ProteinSinthesisProcess():
    def __init__(self, dna_sequences_df, verbose=False):
        self.dna_sequences_df = dna_sequences_df
        self.verbose = verbose

        # add columns to store the results
        columns = ['mrna_sequences', 'polypeptides_chains', 'polypeptides_chains_ext',
            'number_of_proteins_synthesized', 'protein_synthesized']
        for col in columns:
            self.dna_sequences_df[col] = None
        
        self.initialize_environment()
        
    def initialize_environment(self):
        # initialize the simulation environment
        self.dna_sequences = self.dna_sequences_df['sequence'].values
        self.available =  {row['sequence']: True if row['protein_synthesized']==None else False
        for id, row in self.dna_sequences_df.iterrows()}
        
        #random.seed(RANDOM_SEED) TODO: remove or set the seed
        self.env = simpy.Environment()
        self.eucaryotes_cell = EucaryotesCell(environment=self.env, verbose=self.verbose)

        self.process = self.env.process(self.start())
        if self.verbose: print('Simulation environment initialized \t')
    
    def start(self):
        while True:
            yield self.env.timeout(random.random()*10) # time between one protein synthesis and another

            dna_sequence = random.choice(self.dna_sequences) # TODO: Seq object from biopython
            # var: enzimi, basi, ATP, tRNA, aminoacidi
            #atp = random.randint(1,6)

            if self.available[dna_sequence]:
                print(f'Time {self.env.now}: Protein synthesis started')
                #self.env.process(self.eucaryotes_cell.synthesize_protein(dna_sequence))
                self.eucaryotes_cell.synthesize_protein(dna_sequence)
                self.save_proteins_synthesized(
                    dna_sequence,
                    self.eucaryotes_cell.get_mrna(),
                    self.eucaryotes_cell.get_proteins(),
                    self.eucaryotes_cell.get_extended_proteins_name()
                )
                self.available[dna_sequence] = False
                
    def run(self, simulation_time=SIM_TIME):
        if self.verbose: print('Simulation started: \n')
        self.env.run(until=simulation_time)
    
    def save_proteins_synthesized(self, dna_sequence, mrna_sequences, polypeptides_chain, polypeptides_chain_ext):
        # TODO: gestire errori (dna_sequence non trovata nel dataframe)
        row_index = self.dna_sequences_df[self.dna_sequences_df['sequence'] == dna_sequence].index[0]

        results = {
            'dna_sequence': dna_sequence,
            'mrna_sequences': mrna_sequences,
            'polypeptides_chains': polypeptides_chain,
            'polypeptides_chains_ext': polypeptides_chain_ext,
            'number_of_proteins_synthesized': len(mrna_sequences) if mrna_sequences else 0,
            'protein_synthesized': True if mrna_sequences else False
        }
        self.dna_sequences_df.iloc[row_index] = results

        if self.verbose:
            print(f'Protein synthesized: {len(mrna_sequences) if mrna_sequences else 0} \n')

        """if polypeptides_chain:
            if self.verbose: print('Protein synthesized')
            self.dna_sequences_df.loc[row_index, 'protein_synthesized'] = True
            peptides = polypeptides_chain[LENGTH_AMIO_GROUP:-LENGTH_CARBOXYL_GROUP]
            self.dna_sequences_df.loc[row_index, 'peptides_cardinality'] = len(peptides)
        else:
            if self.verbose: print('Protein not synthesized')
            self.dna_sequences_df.loc[row_index, 'protein_synthesized'] = False
            self.dna_sequences_df.loc[row_index, 'peptides_cardinality'] = None
        """