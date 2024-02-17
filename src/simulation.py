import simpy
import random 
from src.protein_synthesis import EucaryotesCell
#from src.utils import load_data_from_pickle, load_data_from_csv

SAMPLES_DATA_PATH = 'data/samples/data.pkl'
HUMAN_GENOMA_DATA_PATH = 'data/human_genoma.csv'
LENGTH_AMIO_GROUP = 4 # length of amino acid group
LENGTH_CARBOXYL_GROUP = 5 # length of carboxyl group
RANDOM_SEED = 42
SIM_TIME = 1000


class ProteinSinthesisProcess():
    def __init__(self, dna_sequences_df, verbose=False):
        self.dna_sequences_df = dna_sequences_df
        self.dna_sequences = dna_sequences_df['sequence'].values
        self.verbose = verbose

        # initialize the simulation environment
        self.eucaryotes_cell = EucaryotesCell()
        if self.verbose: print('Eucaryotes cell initialized')
        self.available = {dna: True for dna in self.dna_sequences} #FIXME
        self.env = simpy.Environment()
        self.process = self.env.process(self.start())
        if self.verbose: print('Simulation environment initialized \t')
    
    def save_synthesize_protein(self, dna_sequence, mrna_sequence, polypeptides_chain, polypeptides_chain_ext):
        row_index = self.dna_sequences_df[self.dna_sequences_df['sequence'] == dna_sequence].index[0]
        
        self.dna_sequences_df.loc[row_index, 'mrna_sequence'] = mrna_sequence
        self.dna_sequences_df.loc[row_index, 'polypeptides_chain_synthetized'] = polypeptides_chain
        self.dna_sequences_df.loc[row_index, 'polypeptides_chain_extended'] = polypeptides_chain_ext

        if polypeptides_chain:
            if self.verbose: print('Protein synthesized')
            self.dna_sequences_df.loc[row_index, 'protein_synthesized'] = True
            peptides = polypeptides_chain[LENGTH_AMIO_GROUP:-LENGTH_CARBOXYL_GROUP]
            self.dna_sequences_df.loc[row_index, 'peptides_cardinality'] = len(peptides)
        else:
            if self.verbose: print('Protein not synthesized')
            self.dna_sequences_df.loc[row_index, 'protein_synthesized'] = False
            self.dna_sequences_df.loc[row_index, 'peptides_cardinality'] = None

    def start(self):
        while True:
            yield self.env.timeout(random.random()*10)

            dna_sequence = random.choice(self.dna_sequences)
            # var: enzimi, basi, ATP, tRNA, aminoacidi
            #atp = random.randint(1,6)

            if self.available[dna_sequence]:
                #self.env.process(self.eucaryotes_cell.synthesize_protein(dna_sequence))
                self.eucaryotes_cell.synthesize_protein(dna_sequence)
                self.save_synthesize_protein(
                    dna_sequence, 
                    self.eucaryotes_cell.get_mrna(),
                    self.eucaryotes_cell.get_protein(),
                    self.eucaryotes_cell.get_extended_protein_name()
                )
                
    def run(self, simulation_time=SIM_TIME):
        if self.verbose: print('Simulation started: \t')
        self.env.run(until=simulation_time)

"""if __name__ == '__main__':
    RESULTS_PATH = 'data/samples/results.pkl'
    SIM_TIME = 1000
    VERBOSE = True

    data = 'sample'
    verbose = VERBOSE

    protein_synthesis_process = ProteinSinthesisProcess(data, verbose)
    protein_synthesis_process.run() # run the simulation
    save_data(protein_synthesis_process.dna_sequences_df, RESULTS_PATH, verbose)
    if verbose: print('Simulation ended')"""