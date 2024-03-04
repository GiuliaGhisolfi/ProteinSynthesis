from src.resources.resource import EukaryoticCellResource
import random
import os
DEV_AMOUNT_TRNA_PER_CODON = 0.1

class TransferRNA:
    """
    This class represents the transfer RNA in the cell. It has a resource for
    each coding codon. The resources are used to request and release the transfer RNA
    in the cell. The resources are saved in a dictionary with the codon as the key.

    Parameters:
    -----------
    environment : simpy.Environment
        The simulation environment
    amount : int
        The amount of transfer RNA for each codon
    codons_list : list
        The list with the coding codons
    random_seed : int
        The random seed for the degradation time
    
    Attributes:
    -----------
    env : simpy.Environment
        The simulation environment
    trna_resources_dict : dict
        The dictionary with the codon as the key and the resource as the value

    Methods:
    --------
    save_history(path_to_save)
        Save the level history of the resources in a json file
    """
    def __init__(self, environment, amount, codons_list, random_seed):
        self.env = environment
        random.seed(random_seed)

        self.trna_resources_dict = dict()
        for codon in codons_list:
            self.trna_resources_dict[codon] = self._init_nucleotide(amount)

    def _init_nucleotide(self, amount):
        trna_amount = random.randint(
            int(amount * (1 - DEV_AMOUNT_TRNA_PER_CODON)),
            int(amount * (1 + DEV_AMOUNT_TRNA_PER_CODON)))
        return EukaryoticCellResource(self.env, capacity=trna_amount)
    
    def save_history(self, path_to_save):
        base_path, ext = os.path.splitext(path_to_save)
        for codons, resource in self.trna_resources_dict.items():
            resource.save_history(base_path+f'_{codons}.json')