import simpy
import random 
import pickle
import pandas as pd
from src.protein_synthesis import EucaryotesCell

SAMPLES_DATA_PATH = 'data/samples/data.pkl'
HUMAN_GENOMA_DATA_PATH = 'data/human_genoma_rna.csv'
LENGTH_AMIO_GROUP = 4 # length of amino acid group
LENGTH_CARBOXYL_GROUP = 5 # length of carboxyl group
RANDOM_SEED = 42
SIM_TIME = 1000

def load_data_from_pickle(path):
    # Load the pickle file
    with open(path, 'rb') as f:
        data = pickle.load(f)

    # Convert the data to a pandas DataFrame
    df = pd.DataFrame(data.values(), columns=['sequence'])
    df['sequence'] = df['sequence'].apply(lambda x: x.upper())

    # add columns
    df['mrna_sequence'] = None
    df['polypeptides_chain_synthetized'] = None
    df['polypeptides_chain_extended'] = None
    df['protein_synthesized'] = None
    df['peptides_cardinality'] = None

    return df

def load_data_from_csv(path):
    # Load the csv file
    df = pd.read_csv(path, index_col=0)

    # convert columns to lowercase
    df.columns = df.columns.str.lower()
    # convert 'rna_sequence' column to 'sequence
    df = df.rename(columns={'rna_sequence': 'sequence'})

    # add columns
    df['mrna_sequence'] = None
    df['polypeptides_chain_synthetized'] = None
    df['polypeptides_chain_extended'] = None
    df['protein_synthesized'] = None
    df['peptides_cardinality'] = None

    return df

def save_data(results_df, path, verbose=False):
    # Save the DataFrame to a pickle file
    with open(path, 'wb') as f:
        pickle.dump(results_df, f)
    
    if verbose: print('Results saved')

class ProteinSinthesisProcess():
    def __init__(self, data='human_genoma', verbose=False):
        self.verbose = verbose

        # load dna sequences
        if data == 'sample':
            self.dna_sequences_df = load_data_from_pickle(SAMPLES_DATA_PATH)
        else:
            if data != 'human_genoma':
                print('Invalid data, using human genoma data as default')
            self.dna_sequences_df = load_data_from_csv(HUMAN_GENOMA_DATA_PATH)
        self.dna_sequences = self.dna_sequences_df['sequence'].tolist()
        if self.verbose: print('DNA sequences loaded')

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