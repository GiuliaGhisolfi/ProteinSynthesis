from src.resources.resource import EucaryotesCellResource
import random
import os
DEV_AMOUNT_TRNA_PER_CODON = 0.1

class TransferRNA:
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
        return EucaryotesCellResource(self.env, capacity=trna_amount)
    
    def save_history(self, path_to_save):
        base_path, ext = os.path.splitext(path_to_save)
        for codons, resource in self.trna_resources_dict.items():
            resource.save_history(base_path+f'_{codons}.json')