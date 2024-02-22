import pickle
import pandas as pd

SAMPLES_DATA_PATH = 'data/samples/data.pkl'
HUMAN_GENOMA_DATA_PATH = 'data/human_genoma.csv'

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

"""
# load dna sequences
    if data == 'sample':
        self.dna_sequences_df = load_data_from_pickle(SAMPLES_DATA_PATH)
    else:
        if data != 'human_genoma':
            print('Invalid data, using human genoma data as default')
        self.dna_sequences_df = load_data_from_csv(HUMAN_GENOMA_DATA_PATH)
    self.dna_sequences = self.dna_sequences_df['sequence'].tolist()
    if self.verbose: print('DNA sequences loaded')
"""

"""
if polypeptides_chain:
    if self.verbose: print('Protein synthesized')
    self.dna_sequences_df.loc[row_index, 'protein_synthesized'] = True
    peptides = polypeptides_chain[LENGTH_AMIO_GROUP:-LENGTH_CARBOXYL_GROUP]
    self.dna_sequences_df.loc[row_index, 'peptides_cardinality'] = len(peptides)
else:
    if self.verbose: print('Protein not synthesized')
    self.dna_sequences_df.loc[row_index, 'protein_synthesized'] = False
    self.dna_sequences_df.loc[row_index, 'peptides_cardinality'] = None
"""