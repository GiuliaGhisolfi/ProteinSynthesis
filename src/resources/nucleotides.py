from src.resources.container import EukaryoticCellContainer
import os
NUCLEOTIDES_NAMES = ['uracil', 'adenine', 'guanine', 'cytosine']

class Nucleotides:
    """
    This class represents the nucleotides in the cell. It has a container for
    each nucleotide. The containers are used to request and release the nucleotides
    in the cell. The containers are saved in a dictionary with the nucleotide name
    as the key.

    Parameters:
    -----------
    environment : simpy.Environment
        The simulation environment
    uracil_initial_amount : int
        The initial amount of uracil in the cell
    adenine_initial_amount : int
        The initial amount of adenine in the cell
    guanine_initial_amount : int
        The initial amount of guanine in the cell
    cytosine_initial_amount : int
        The initial amount of cytosine in the cell
    random_seed : int
        The random seed for the degradation time

    Attributes:
    -----------
    env : simpy.Environment
        The simulation environment
    nucleotides_containers_dict : dict
        The dictionary with the nucleotide name as the key and the container as the value
    
    Methods:
    --------
    request(nucleotide, amount)
        Request the amount of the nucleotide
    release(nucleotide, amount)
        Release the amount of the nucleotide
    save_history(path_to_save)
        Save the level history of the containers in a json file
    """
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
        return EukaryoticCellContainer(
            self.env, capacity=float('inf'), init=amount, random_seed=random_seed)

    def request(self, nucleotide, amount):
        return self.nucleotides_containers_dict[nucleotide].get(amount)
    
    def release(self, nucleotide, amount):
        self.env.process(self.nucleotides_containers_dict[nucleotide].put(amount))
    
    def save_history(self, path_to_save):
        base_path, ext = os.path.splitext(path_to_save)
        for nucleotide, container in zip(NUCLEOTIDES_NAMES, self.nucleotides_containers_dict.values()):
            container.save_history(base_path + f'_{nucleotide}.json')