from src.resources.container import EucaryotesCellContainer

class Nucleotides:
    def __init__(self, environment, uracil_initial_amount, adenine_initial_amount, 
            guanine_initial_amount, cytosine_initial_amount):
        self.env = environment

        self.nucleotides_containers_dict = {
            'U': self._init_nucleotide(uracil_initial_amount), # uracil
            'A': self._init_nucleotide(adenine_initial_amount), # adenine
            'G': self._init_nucleotide(guanine_initial_amount), # guanine
            'C': self._init_nucleotide(cytosine_initial_amount), # cytosine
        }

    def _init_nucleotide(self, amount):
        return EucaryotesCellContainer(self.env, capacity=float('inf'), init=amount)

    def request(self, nucleotide, amount):
        return self.nucleotides_containers_dict[nucleotide].get(amount)
    
    def release(self, nucleotide, amount):
        return self.nucleotides_containers_dict[nucleotide].put(amount)
    
    def save_history(self, path_to_save):
        for nucleotide in self.nucleotides_containers_dict:
            self.nucleotides_containers_dict[nucleotide].save_history(
                path_to_save.replace('.json', f'_{nucleotide}.json'))