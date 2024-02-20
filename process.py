import simpy
import random

SEQUENCE = ['ATCG', 'AAAA', 'TTTT', 'GGGG', 'CCCC']
# call using yield to wait for the process to complete

class Cell():
    def __init__(self, env):
        self.env = env
        self.ribosome = simpy.Resource(self.env, capacity=2)

    def protein_synthesis(self, dna_sequence):
        with self.ribosome.request() as request:
            yield request
            print(f'Time {self.env.now}: Ribosome requested for protein synthesis')

            protein = dna_sequence
            yield self.env.timeout(4)  # Wait for 1 unit of time for protein synthesis
            print(f'Time {self.env.now}: Protein {protein} synthesis completed')
            return protein

def setup(env):
    cell = Cell(env)
    results = []
    while True:
        dna_sequence = random.choice(SEQUENCE)
        #print(f'Time {env.now}: Protein synthesis started')
        protein = env.process(cell.protein_synthesis(dna_sequence))
        results.append(protein)
        yield env.timeout(1)

if __name__ == '__main__':
    env = simpy.Environment()
    env.process(setup(env))
    print('Simulation started')
    env.run(until=15)
    print('Simulation ended')
