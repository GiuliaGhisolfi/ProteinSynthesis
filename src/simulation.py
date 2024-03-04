import simpy
import random 
import itertools
import json
import os
from src.process.protein_synthesis import EukaryoticCell
from src.variables.variables import EukaryoticCellVariables
from src.resources.resource import EukaryoticCellResource
from src.utils.utils import save_proteins_synthesized, post_processing_results

LENGTH_AMIO_GROUP = 4 # length of amino acid group
LENGTH_CARBOXYL_GROUP = 5 # length of carboxyl group
DATA_PATH = 'data/'
CODONS_PATH = DATA_PATH + 'codons.json'
RESULTS_FOLDER = 'results/'

SIM_TIME = 1000
NUMBER_RESOURCES = 200
NUMBER_RNA_POLYMERASES = 3
NUMBER_RIBOSOMES = 2
NUMBER_RNA_TRANSFER = 2000
URACIL_INITIAL_AMOUNT = 5000
ADENINE_INITIAL_AMOUNT = 5000
GUANINE_INITIAL_AMOUNT = 5000
CYTOSINE_INITIAL_AMOUNT = 5000
RANDOM_SEED = None

class ProteinSinthesisProcess:
    """
    Class to simulate the protein synthesis process in eukaryotic cells.

    Parameters
    ----------
    dna_sequences_df: pandas.DataFrame
        DataFrame containing a column 'sequence' with the DNA sequences to be synthesized.
    number_resources: int, optional
        Number of resources available in the simulation environment. The default is 200.
    number_rna_polymerases: int, optional
        Number of RNA polymerases in the simulation environment. The default is 3.
    number_ribosomes: int, optional
        Number of ribosomes in the simulation environment. The default is 2.
    number_rna_transfers_per_codon: int, optional
        Number of RNA transfer per codon in the simulation environment is a random integer in 
        [0.9*number_rna_transfers_per_codon, 1.1*number_rna_transfers_per_codon]. The default is 2000.
    uracil_initial_amount: int, optional
        Initial amount of uracil in the simulation environment. The default is 5000.
    adenine_initial_amount: int, optional
        Initial amount of adenine in the simulation environment. The default is 5000.
    guanine_initial_amount: int, optional
        Initial amount of guanine in the simulation environment. The default is 5000.
    cytosine_initial_amount: int, optional
        Initial amount of cytosine in the simulation environment. The default is 5000.
    random_seed: int, optional
        Random seed to reproduce the results. The default is None.
    verbose: bool, optional
        If True, print the simulation process. The default is False.

    Methods
    -------
    run(simulation_time)
        Run the simulation process.
    save_process(folder_test_name)
        Save the results of the simulation process.
    """
    def __init__(self, 
            dna_sequences_df, 
            number_resources=NUMBER_RESOURCES,
            number_rna_polymerases=NUMBER_RNA_POLYMERASES, 
            number_ribosomes=NUMBER_RIBOSOMES,
            number_rna_transfers_per_codon=NUMBER_RNA_TRANSFER, 
            uracil_initial_amount=URACIL_INITIAL_AMOUNT, 
            adenine_initial_amount=ADENINE_INITIAL_AMOUNT, 
            guanine_initial_amount=GUANINE_INITIAL_AMOUNT, 
            cytosine_initial_amount=CYTOSINE_INITIAL_AMOUNT,
            random_seed=RANDOM_SEED, verbose=False):
        self.dna_sequences_df = dna_sequences_df
        self.verbose = verbose

        # nucleotides
        self.uracil_initial_amount = uracil_initial_amount
        self.adenine_initial_amount = adenine_initial_amount
        self.guanine_initial_amount = guanine_initial_amount
        self.cytosine_initial_amount = cytosine_initial_amount

        # codons to translate into amino acids
        self.codons = list(json.load(open(CODONS_PATH)).keys())

        # add columns to store the results
        columns = ['mrna_sequences', 'polypeptides_chains', 'polypeptides_chains_ext',
            'length_mrna_sequences', 'number_of_proteins_synthesized_per_mrna', 
            'number_of_proteins_synthesized', 'protein_synthesized', 'request_start_process_time',
            'start_process_time', 'start_transcription_time', 'start_translation_time',
            'end_translation_time', 'end_process_time']
        for col in columns:
            self.dna_sequences_df[col] = None
        
        # initialize the simulation environment
        self.dna_sequences = self.dna_sequences_df['sequence'].values
        self.available =  {row['sequence']: True if row['protein_synthesized']==None else False
            for _, row in self.dna_sequences_df.iterrows()}
        
        random.seed(random_seed) # set the random seed
        self.env = simpy.Environment() # create the Simpy simulation environment
        self.resources = EukaryoticCellResource(
            self.env, capacity=number_resources, save_history=False) 
        self.env.process(self._setup_process())

        self.eukaryotic_cell = EukaryoticCell(
            environment=self.env, 
            number_rna_polymerases=number_rna_polymerases,
            number_ribosomes=number_ribosomes, 
            number_rna_transfers_per_codon=number_rna_transfers_per_codon,
            uracil_initial_amount=uracil_initial_amount, 
            adenine_initial_amount=adenine_initial_amount, 
            guanine_initial_amount=guanine_initial_amount,
            cytosine_initial_amount=cytosine_initial_amount, 
            random_seed=random_seed, 
            verbose=self.verbose
            )
        
        print('Simulation environment initialized, time unit: 0.0001 second.')
    
    def __str__(self):
        return (f'Protein Sinthesis Process:\n'
            f'{len(self.dna_sequences)} dna sequences to synthesize,\n'
            f'{self.resources.capacity} resources available,\n'
            f'{self.eukaryotic_cell.nucleus.rna_polymerase.capacity} RNA polymerases,\n'
            f'{self.eukaryotic_cell.ribosome.ribosomes.capacity} ribosomes,\n'
            f'{self.uracil_initial_amount} uracil bases,\n'
            f'{self.adenine_initial_amount} adenine bases,\n'
            f'{self.guanine_initial_amount} guanine bases,\n'
            f'{self.cytosine_initial_amount} cytosine bases.')
    
    def __repr__(self):
        trna_info = ''.join([f",\n{self.eukaryotic_cell.ribosome.rna_transfer.trna_resources_dict[codon].capacity} transfer RNA for {codon} codon"for codon in self.codons])
        return (f'Protein Sinthesis Process:\n'
            f'{len(self.dna_sequences)} dna sequences to synthesize,\n'
            f'{self.resources.capacity} resources available,\n'
            f'{self.eukaryotic_cell.nucleus.rna_polymerase.capacity} RNA polymerases,\n'
            f'{self.eukaryotic_cell.ribosome.ribosomes.capacity} ribosomes,\n'
            f'{self.uracil_initial_amount} uracil bases,\n'
            f'{self.adenine_initial_amount} adenine bases,\n'
            f'{self.guanine_initial_amount} guanine bases,\n'
            f'{self.cytosine_initial_amount} cytosine bases'
            f'{trna_info}.')
    
    def run(self, simulation_time=SIM_TIME):
        """
        Run the simulation process.

        Parameters
        ----------
        simulation_time: int
            Time to run the simulation process in seconds.
        """
        print('Simulation started')
        self.env.run(until=simulation_time)

        # save simulation results
        proteins_number = self.dna_sequences_df[self.dna_sequences_df[
            'protein_synthesized'].notna()]['number_of_proteins_synthesized'].sum()
        dna_sequences_processed_number = self.dna_sequences_df[
            self.dna_sequences_df['protein_synthesized'].notna()].shape[0]
        
        print(f'End simulation: {proteins_number} proteins synthesized from '
            f'{dna_sequences_processed_number} DNA sequences.')
    
    def _setup_process(self):
        """
        Setup the simulation process.
        Chose the dna sequences, check the availability of the dna sequences to 
        be synthesized and start the protein synthesis process.
        """
        process_queue = []
        sequences_count = itertools.count()

        while True:
            dna_sequence = random.choice(self.dna_sequences)
            if self.available[dna_sequence]:
                # initialize the variables related to the dna sequence
                variables = EukaryoticCellVariables()
                variables.dna_sequence = dna_sequence
                variables.sequence_count = next(sequences_count)

                process_queue.append(self.env.process(self._process(variables)))
                
                self.available[dna_sequence] = False
                # time between start of protein synthesis
                yield self.env.timeout(round(random.random()*10, ndigits=4))

                while process_queue: # wait for all the protein synthesis to be completed
                    process_queue.pop(0)
        
    def _process(self, variables):
        """
        Start the protein synthesis process.

        Parameters
        ----------
        variables: EukaryoticCellVariables
            Variables related to the dna sequence to be synthesized.
        """
        # Synthesize dna sequences while the simulation is running       
        with self.resources.request() as request:
            if self.verbose:
                print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} requesting to start synthesis')
            variables.request_start_process_time = self.env.now

            yield request # wait for a cell be able to accepts dna sequence

            if self.verbose:
                print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} synthesize started')
            variables.start_process_time = self.env.now

            yield self.env.process(self.eukaryotic_cell.synthesize_protein(variables))
            variables.end_process_time = self.env.now

            # save the results
            self._save_proteins_synthesized_in_df(variables)

            if self.verbose:
                print(f'Time {self.env.now:.4f}: DNA Sequence {variables.sequence_count} synthetis ended')
        
    def _save_proteins_synthesized_in_df(self, variables):
        """
        Save the results of the protein synthesis process in the dataframe.

        Parameters
        ----------
        variables: EukaryoticCellVariables
            Variables related to the dna sequence to be synthesized.
        """
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
            end_process_time=variables.end_process_time,
            promoters_box=variables.promoters_box,
            proteins_sintetized=variables.proteins_sintetized
            )
        
    def save_process(self, folder_test_name=''):
        """
        Save the results of the simulation process.

        Parameters
        ----------
        folder_test_name: str, optional
            By default, the results are saved in the root of the results folder.
        """
        # create folder to save the results
        if folder_test_name != '':
            if not os.path.exists(RESULTS_FOLDER+folder_test_name):
                os.mkdir(RESULTS_FOLDER+folder_test_name)
            folder_test_name = folder_test_name + '/'

        # save dataframe with the results    
        df_to_save = self.dna_sequences_df[self.dna_sequences_df['protein_synthesized'].notna()]
        df_to_save = df_to_save.apply(post_processing_results, axis=1)
        df_to_save.to_csv(RESULTS_FOLDER+folder_test_name+'results.csv')

        # save resources history 
        self.eukaryotic_cell.nucleus.rna_polymerase.save_history(
            RESULTS_FOLDER+folder_test_name+'rna_polymerase_history.json')
        
        self.eukaryotic_cell.ribosome.ribosomes.save_history(
            RESULTS_FOLDER+folder_test_name+'ribosome_history.json')
        
        if not os.path.exists(RESULTS_FOLDER+folder_test_name+'nucleotides'):
            os.mkdir(RESULTS_FOLDER+folder_test_name+'nucleotides')
        self.eukaryotic_cell.nucleotides.save_history(
            RESULTS_FOLDER+folder_test_name+'nucleotides/'+'nucleotides_history.json')

        if not os.path.exists(RESULTS_FOLDER+folder_test_name+'rna_transfer'):
            os.mkdir(RESULTS_FOLDER+folder_test_name+'rna_transfer')
        self.eukaryotic_cell.ribosome.rna_transfer.save_history(
            RESULTS_FOLDER+folder_test_name+'rna_transfer/'+'rna_transfer_history.json')
        
        print('Process saved.')
