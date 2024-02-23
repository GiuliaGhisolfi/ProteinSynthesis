from src.simulation import ProteinSinthesisProcess
from HumanGenomeDataset.load_dataset import load_dataset

SIM_TIME = 300
NUMBER_RESOURCES = 5
NUMBER_RNA_POLYMERASES = 3
NUMBER_RIBOSOMES = 2
RANDOM_SEED = 42

if __name__ == '__main__':
    data_df = load_dataset('dna_protein_coding_sequences')

    protein_synthesis_process = ProteinSinthesisProcess(
            dna_sequences_df=data_df,
            number_resources=NUMBER_RESOURCES,
            number_rna_polymerases=NUMBER_RNA_POLYMERASES, 
            number_ribosomes=NUMBER_RIBOSOMES,
            random_seed=RANDOM_SEED, 
            verbose=False
            )
    protein_synthesis_process.run(simulation_time=SIM_TIME) # run the simulation