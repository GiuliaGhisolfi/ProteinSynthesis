import simpy
import random 
import pickle
import pandas as pd
from src.protein_synthesis import EucaryotesCell

DATA_PATH = 'data/data.pkl'
LENGTH_AMIO_GROUP = 4 # length of amino acid group
LENGTH_CARBOXYL_GROUP = 5 # length of carboxyl group
RANDOM_SEED = 42
SIM_TIME = 1000

def load_data(path):
    # Load the pickle file
    with open(DATA_PATH, 'rb') as f:
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

class ProteinSinthesisProcess():
    def __init__(self):
        # load dna sequences
        self.dna_sequences_df = load_data(DATA_PATH)
        self.dna_sequences = self.dna_sequences_df['sequence'].tolist()
        print('DNA sequences loaded')

        # initialize the simulation environment
        self.eucaryotes_cell = EucaryotesCell()
        print('Eucaryotes cell initialized')
        self.available = {dna: True for dna in self.dna_sequences} #FIXME
        self.env = simpy.Environment()
        self.process = self.env.process(self.start())
        print('Simulation environment initialized \t')
    
    def save_synthesize_protein(self, dna_sequence, mrna_sequence, polypeptides_chain, polypeptides_chain_ext):
        row_index = self.dna_sequences_df[self.dna_sequences_df['sequence'] == dna_sequence].index[0]
        
        self.dna_sequences_df.loc[row_index, 'mrna_sequence'] = mrna_sequence
        self.dna_sequences_df.loc[row_index, 'polypeptides_chain_synthetized'] = polypeptides_chain
        self.dna_sequences_df.loc[row_index, 'polypeptides_chain_extended'] = polypeptides_chain_ext

        if polypeptides_chain:
            print('Protein synthesized')
            self.dna_sequences_df.loc[row_index, 'protein_synthesized'] = True
            peptides = polypeptides_chain[LENGTH_AMIO_GROUP:-LENGTH_CARBOXYL_GROUP]
            self.dna_sequences_df.loc[row_index, 'peptides_cardinality'] = len(peptides)
        else:
            print('Protein not synthesized')
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
                
    def run(self):
        print('Simulation started: \t')
        self.env.run(until=SIM_TIME)

if __name__ == '__main__':
    protein_synthesis_process = ProteinSinthesisProcess()
    protein_synthesis_process.run() # run the simulation
    print('Simulation ended')