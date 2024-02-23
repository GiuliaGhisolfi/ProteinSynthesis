from src.resources.container import EucaryotesCellContainer

class Nucleotides:
    def __init__(self, environment, uracil, adenine, guanine, cytosine):
        self.env = environment

        self.nucleotides = {
            'U': self._init_nucleotide(uracil), # uracil
            'A': self._init_nucleotide(adenine), # adenine
            'G': self._init_nucleotide(guanine), # guanine
            'C': self._init_nucleotide(cytosine), # cytosine
        }

    def _init_nucleotide(self, amount):
        return EucaryotesCellContainer(self.env, capacity=float('inf'), init=amount)

    def request(self, nucleotide, amount):
        return self.nucleotides[nucleotide].get(amount)
    
    def release(self, nucleotide, amount):
        return self.nucleotides[nucleotide].put(amount)