from src.simulation import ProteinSinthesisProcess
from HumanGenomeDataset.load_dataset import load_dataset
SIM_TIME = 30

if __name__ == '__main__':
    data_df = load_dataset('dna_protein_coding_sequences')

    protein_synthesis_process = ProteinSinthesisProcess(dna_sequences_df=data_df, verbose=False)
    protein_synthesis_process.run(simulation_time=SIM_TIME) # run the simulation