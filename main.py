from src.simulation import ProteinSinthesisProcess
from HumanGenomeDataset.load_dataset import load_dataset

SIM_TIME = 1000
NUMBER_RESOURCES = 5
NUMBER_RNA_POLYMERASES = 3
NUMBER_RIBOSOMES = 2
NUMBER_RNA_TRANSFER = 800
URACIL_INITIAL_AMOUNT = 50000
ADENINE_INITIAL_AMOUNT = 50000
GUANINE_INITIAL_AMOUNT = 50000
CYTOSINE_INITIAL_AMOUNT = 50000
RANDOM_SEED = 42

if __name__ == '__main__':
    data_df = load_dataset('dna_protein_coding_sequences')

    protein_synthesis_process = ProteinSinthesisProcess(
        dna_sequences_df=data_df,
        number_resources=NUMBER_RESOURCES,
        number_rna_polymerases=NUMBER_RNA_POLYMERASES, 
        number_ribosomes=NUMBER_RIBOSOMES,
        number_rna_transfers_per_codon=NUMBER_RNA_TRANSFER,
        uracil_initial_amount=URACIL_INITIAL_AMOUNT,
        adenine_initial_amount=ADENINE_INITIAL_AMOUNT,
        guanine_initial_amount=GUANINE_INITIAL_AMOUNT,
        cytosine_initial_amount=CYTOSINE_INITIAL_AMOUNT,
        random_seed=RANDOM_SEED,
        verbose=True
        )
    protein_synthesis_process.run(simulation_time=SIM_TIME) # run the simulation
    protein_synthesis_process.save_process(folder_test_name='test')