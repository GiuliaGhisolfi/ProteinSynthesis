from src.simulation import ProteinSinthesisProcess
from HumanGenomeDataset.load_dataset import load_dataset

SIM_TIME = 36
NUMBER_RESOURCES = 5
NUMBER_RNA_POLYMERASES = 3
NUMBER_RIBOSOMES = 2
URACIL_INITIAL_AMOUNT = 10#5e100
ADENINE_INITIAL_AMOUNT = 5e100
GUANINE_INITIAL_AMOUNT = 5e100
CYTOSINE_INITIAL_AMOUNT = 5e100
RANDOM_SEED = None

if __name__ == '__main__':
    data_df = load_dataset('dna_protein_coding_sequences')

    protein_synthesis_process = ProteinSinthesisProcess(
        dna_sequences_df=data_df,
        number_resources=NUMBER_RESOURCES,
        number_rna_polymerases=NUMBER_RNA_POLYMERASES, 
        number_ribosomes=NUMBER_RIBOSOMES,
        uracil_initial_amount=URACIL_INITIAL_AMOUNT,
        adenine_initial_amount=ADENINE_INITIAL_AMOUNT,
        guanine_initial_amount=GUANINE_INITIAL_AMOUNT,
        cytosine_initial_amount=CYTOSINE_INITIAL_AMOUNT,
        random_seed=RANDOM_SEED,
        verbose=True
        )
    protein_synthesis_process.run(simulation_time=SIM_TIME) # run the simulation
    protein_synthesis_process.save_process()