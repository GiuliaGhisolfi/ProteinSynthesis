import random
import simpy
from src.nucleotides import NucleotidesSymbolsAllocations
from src.resources.resource import EucaryotesCellResource

BASE_COMPLEMENT_DNA2RNA = {
    'A': 'U', 
    'T': 'A', 
    'C': 'G', 
    'G': 'C',
}
BASE_COMPLEMENT_RNA2DNA = {
    'U': 'A', 
    'A': 'T', 
    'G': 'C', 
    'C': 'G'
}
PROMOTERS = [
    'TATAAAA', 'TATAAAT', 'TATATAA', 'TATATAT', # TATA box
]
"""    
    'CAAT', # CAAT box
    'GC', # GC box
    'GGGCGG' # Initiator element
]
"""
LENGTH_PROMOTER = 7
TERMINATORS = ['UAA', 'UAG', 'UGA']
RNA_POLYMERASE_ERROR_RATE = 10e-4 # 1 error per 10^4 nucleotides
LENGTH_EXTRON_SEQUENCE = 3 # length of extron sequence
LENGTH_METHYL_CAP = 8 # length of 5'-methyl cap
NUMBER_RNA_POLYMERASES = 3

class Nucleus:

    def __init__(self, environment, extron_sequences_list, editing_sites_dict):
        self.env = environment
        self.extron_sequences_list = extron_sequences_list

        self.editing_sites_dict = editing_sites_dict
        self.editing_sites_dict = dict(sorted(self.editing_sites_dict.items(), 
            key=lambda x: len(x[0]), reverse=False)) # sort by length of key
        
        self.rna_polymerase = EucaryotesCellResource(self.env, capacity=NUMBER_RNA_POLYMERASES)
        self.nucleotides = {'U': 0, 'A': 0, 'G': 0, 'C': 0} 
        # TODO: change e implementare il conteggio prima e dopo degradetion, nel ribosoma (o gestire dentro cell)
    
    def find_promoter(self, dna_sequence):
        # find promoter sequences in the DNA sequence
        promoter_positions_list = sorted([i for promoter in PROMOTERS for i, _ in 
            enumerate(dna_sequence) if dna_sequence[i:].startswith(promoter)])
        
        if len(promoter_positions_list) == 0:
            return None
        else:
            # split the DNA sequence in the promoter regions
            dna_sequences_list = []

            for i in range(len(promoter_positions_list)-1):
                dna_sequences_list.append(dna_sequence[promoter_positions_list[i]+LENGTH_PROMOTER
                    :promoter_positions_list[i+1]])
                
            dna_sequences_list.append(dna_sequence[promoter_positions_list[-1]+LENGTH_PROMOTER:])
            
            return dna_sequences_list

    def transcript(self, dna_sequence): # enzime: RNA polymerase
        with self.rna_polymerase.request() as request:
            yield request  # wait for RNA polymerase to be available

            # start transcript processes for DNA sequence
            messenger_rna_sequence = yield self.env.process(self.transcript_process(dna_sequence))

        return messenger_rna_sequence

    def transcript_process(self, dna_sequence):
        # make sequence univoque to transcript
        dna_sequence = ''.join([random.choice(NucleotidesSymbolsAllocations[n]) for n in dna_sequence])

        # transcript from gene to pre-mRNA
        messenger_rna_sequence = ''.join([BASE_COMPLEMENT_DNA2RNA[base] 
            if random.random() > RNA_POLYMERASE_ERROR_RATE 
            else random.choice([b for b in list(BASE_COMPLEMENT_DNA2RNA.values()) if b != BASE_COMPLEMENT_DNA2RNA[base]])
            for base in dna_sequence])

        messenger_rna_sequence = self.capping(messenger_rna_sequence)

        # elongation phase
        messenger_rna_sequence = self.splicing(messenger_rna_sequence)
        messenger_rna_sequence = self.editing(messenger_rna_sequence)
        
        # post-transcriptional modifications
        messenger_rna_sequence = self.cleavage(messenger_rna_sequence)
        messenger_rna_sequence = self.polyadenylation(messenger_rna_sequence)

        yield self.env.timeout(1) #TODO: implementare il tempo di trascrizione

        return messenger_rna_sequence
        
    def splicing(self, rna_sequence):
        # remove introns: non-coding regions
        i = LENGTH_METHYL_CAP # index
        while i+3 < len(rna_sequence):
            if rna_sequence[i:i+LENGTH_EXTRON_SEQUENCE] in self.extron_sequences_list:
                i += LENGTH_EXTRON_SEQUENCE
            else: # TODO: implementare il conteggio dei nucleotidi o togliere tutto
                # remouve nucleotide and save it
                self.nucleotides[rna_sequence[i]] += 1
                rna_sequence = rna_sequence[:i] + rna_sequence[i+1:]

        return rna_sequence

    def editing(self, rna_sequence):
        for editing_site in self.editing_sites_dict.keys():
            rna_sequence = rna_sequence.replace(editing_site, self.editing_sites_dict[editing_site])
        return rna_sequence

    def capping(self, rna_sequence):
        return 'CH3GPPP-{}'.format(rna_sequence) # Add 5'-methyl cap
    
    def cleavage(self, rna_sequence): #TODO
        # first step in adding a polyadenine tail to the pre-mRNA (post-transcriptional
        # modifications) it is necessary for producing a mature mRNA molecule
        return rna_sequence

    def polyadenylation(self, rna_sequence):
        return '{}-AAAA'.format(rna_sequence) # Add PolyA tail (250 nucleotides circa)
    