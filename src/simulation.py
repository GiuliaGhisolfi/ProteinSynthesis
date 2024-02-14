import simpy 
import random 
import pickle
import pandas as pd
from src.protein_synthesis import EucaryotesCell
from collections import namedtuple

DATA_PATH = 'data/data.pkl'
LENGTH_AMIO_GROUP = 4 # length of amino acid group
LENGTH_CARBOXYL_GROUP = 5 # length of carboxyl group
RANDOM_SEED = 42
SIM_TIME = 100

def load_data(path):
    # Load the pickle file
    with open(DATA_PATH, 'rb') as f:
        data = pickle.load(f)

    # Convert the data to a pandas DataFrame
    df = pd.DataFrame(data.values(), columns=['sequence'])
    df['sequence'] = df['sequence'].apply(lambda x: x.upper())

    return df

def synthesize_protein(cell, df_row):
    cell.synthesize_protein(df_row['sequence'])

    polypeptides_chain = cell.get_protein()
    polypeptides_chain_ext = cell.get_extended_protein_name()

    df_row['polypeptides_chain_synthetized'] = polypeptides_chain
    df_row['polypeptides_chain_extended'] = polypeptides_chain_ext

    if polypeptides_chain:
        df_row['protein_synthesized'] = True
        peptides = polypeptides_chain[LENGTH_AMIO_GROUP:-LENGTH_CARBOXYL_GROUP]
        df_row['peptides_cardinality'] = len(peptides)
    else:
        df_row['protein_synthesized'] = False
        df_row['peptides_cardinality'] = None
    
    return df_row

def from_gene_to_protein(env, cell, df):
    while True:
        yield env.timeout(random.random()*10)

        dna_sequence = random.choice(cell.foods)
        # var: enzimi, basi, ATP, tRNA, aminoacidi
        #atp = random.randint(1,6)

        if cell.available[dna_sequence]:
            env.process(cell(env, dna_sequence, atp, cell))


# Set up and start the simulation
random.seed(RANDOM_SEED)
env = simpy.Environment()
eucaryotes_cell = EucaryotesCell()

dna_sequences_df = load_data(DATA_PATH)
dna_sequences = dna_sequences_df['sequence'].values

"""# Create environment and start processes
Nucleus = simpy.Resource(EucaryotesCell) #TODO: add capacity
Ribosome = simpy.Resource(EucaryotesCell)

gene_found = {dna_sequence: EucaryotesCell.event() for dna_sequence in dna_sequences}
gene_not_found = {dna_sequence: None for dna_sequence in dna_sequences}

Cell = namedtuple('Nucleus, Ribosome', 'gene_found, gene_not_found') #FIXME
    # not_enough_atp, not_enough_enzimes, not_enough_nucleotides')
cell = Cell(Nucleus, Ribosome, gene_found, gene_not_found)"""

# Start process and run
EucaryotesCell.process(from_gene_to_protein(env))
EucaryotesCell.run(until=SIM_TIME)