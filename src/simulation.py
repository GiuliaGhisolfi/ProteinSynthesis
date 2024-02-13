import simpy 
import random 
from collections import namedtuple

RANDOM_SEED = 42
SIM_TIME = 100

def from_gene_to_protein(env, cell):
    while True:
        yield env.timeout(random.random()*10)

        dna_sequence = random.choice(cell.foods)
        # var: enzimi, basi, ATP, tRNA, aminoacidi
        atp = random.randint(1,6)

        if cell.available[dna_sequence]:
            env.process(cell(env, dna_sequence, atp, cell))

def extracted_gene(dna_sequences):
    return dna_sequences

# Set up and start the simulation
random.seed(RANDOM_SEED)
EucaryotesCell = simpy.Environment()

# Create environment and start processes
Nucleus = simpy.Resource(EucaryotesCell) #TODO: add capacity
Ribosome = simpy.Resource(EucaryotesCell)

dna_sequences = [] # TODO: read from file
genes = extracted_gene(dna_sequences)

gene_found = {gene: EucaryotesCell.event() for gene in genes}
gene_not_found = {gene: None for gene in genes}

Cell = namedtuple('Nucleus, Ribosome', 'gene_found, gene_not_found') #FIXME
    # not_enough_atp, not_enough_enzimes, not_enough_nucleotides')
cell = Cell(Nucleus, Ribosome, gene_found, gene_not_found)

# Start process and run
EucaryotesCell.process(from_gene_to_protein(EucaryotesCell, cell))
EucaryotesCell.run(until=SIM_TIME)