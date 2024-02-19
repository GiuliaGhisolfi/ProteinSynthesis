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
        # add columns to store the results
        self.dna_sequences_df['mrna_sequences'] = None
        self.dna_sequences_df['polypeptides_chains'] = None
        self.dna_sequences_df['polypeptides_chains_ext'] = None
        self.dna_sequences_df['number_of_proteins_synthesized'] = None
        self.dna_sequences_df['protein_synthesized'] = None

        self.dna_sequences = dna_sequences_df['sequence'].values
        self.verbose = verbose

        # initialize the simulation environment
        self.eucaryotes_cell = EucaryotesCell()
        if self.verbose: print('Eucaryotes cell initialized')

        self.available = {dna: True for dna in self.dna_sequences} #FIXME
        self.env = simpy.Environment()

        self.process = self.env.process(self.start())
        if self.verbose: print('Simulation environment initialized \t')
    
    def save_synthesize_protein(self, dna_sequence, mrna_sequences, polypeptides_chain, polypeptides_chain_ext):
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
            print(f'Protein synthesized: {len(mrna_sequences) if mrna_sequences else 0}')

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
                    self.eucaryotes_cell.get_proteins(),
                    self.eucaryotes_cell.get_extended_proteins_name()
                )
                
    def run(self, simulation_time=SIM_TIME):
        if self.verbose: print('Simulation started: \t')
        self.env.run(until=simulation_time)