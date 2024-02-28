from src.resources.container import EucaryotesCellContainer
import os

NUCLEOTIDES_NAMES = ['uracil', 'adenine', 'guanine', 'cytosine']

class Nucleotides:
    def __init__(self, environment, uracil_initial_amount, adenine_initial_amount, 
            guanine_initial_amount, cytosine_initial_amount, random_seed):
        self.env = environment

        self.nucleotides_containers_dict = {
            'U': self._init_nucleotide(uracil_initial_amount, random_seed), # uracil
            'A': self._init_nucleotide(adenine_initial_amount, random_seed), # adenine
            'G': self._init_nucleotide(guanine_initial_amount, random_seed), # guanine
            'C': self._init_nucleotide(cytosine_initial_amount, random_seed), # cytosine
        }

    def _init_nucleotide(self, amount, random_seed):
        return EucaryotesCellContainer(
            self.env, capacity=float('inf'), init=amount, random_seed=random_seed)

    def request(self, nucleotide, amount):
        return self.nucleotides_containers_dict[nucleotide].get(amount)
    
    def release(self, nucleotide, amount):
        self.env.process(self.nucleotides_containers_dict[nucleotide].put(amount))
        #self.nucleotides_containers_dict[nucleotide].put(amount)
    
    def save_history(self, path_to_save):
        base_path, ext = os.path.splitext(path_to_save)
        for nucleotide, container in zip(NUCLEOTIDES_NAMES, self.nucleotides_containers_dict.values()):
            container.save_history(base_path + f'_{nucleotide}.json')