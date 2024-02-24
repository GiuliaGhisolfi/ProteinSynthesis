import random
from src.variables.nucleotides_allocations import NucleotidesSymbolsAllocations
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
LENGTH_PROMOTER = 7
TERMINATORS = ['UAA', 'UAG', 'UGA']
RNA_POLYMERASE_ERROR_RATE = 10e-4 # 1 error per 10^4 nucleotides
REPLICATION_TIME = 2e-2 # seconds to replicate a nucleotide
LENGTH_EXTRON_SEQUENCE = 3 # length of extron sequence
LENGTH_METHYL_CAP = 8 # length of 5'-methyl cap

class Nucleus:
    def __init__(self, environment, extron_sequences_list, editing_sites_dict, 
            number_rna_polymerases, nucleotides, random_seed):
        self.env = environment

        self.extron_sequences_list = extron_sequences_list
        self.editing_sites_dict = editing_sites_dict
        self.editing_sites_dict = dict(sorted(self.editing_sites_dict.items(), 
            key=lambda x: len(x[0]), reverse=False)) # sort by length of key
        
        self.rna_polymerase = EucaryotesCellResource(self.env, capacity=number_rna_polymerases)
        self.nucleotides = nucleotides

        self.random_seed = random_seed
    
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

    def transcript(self, dna_sequence, variables, sequence_count): # enzime: RNA polymerase
        with self.rna_polymerase.request() as request:
            yield request  # wait for RNA polymerase to be available

            # start transcript processes for DNA sequence
            messenger_rna_sequence = yield self.env.process(
                self.transcript_process(dna_sequence, variables, sequence_count))

        return messenger_rna_sequence

    def transcript_process(self, dna_sequence, variables, sequence_count):
        # make sequence univoque to transcript
        random.seed(self.random_seed)
        dna_sequence = ''.join([random.choice(NucleotidesSymbolsAllocations[n]) for n in dna_sequence])

        # transcript from gene to pre-mRNA
        messenger_rna_sequence = yield self.env.process(
            self.trascript_gene(dna_sequence, variables, sequence_count))

        # pre-transcriptional modifications
        messenger_rna_sequence = self.capping(messenger_rna_sequence)

        # elongation phase
        messenger_rna_sequence = self.splicing(messenger_rna_sequence)
        messenger_rna_sequence = self.editing(messenger_rna_sequence)
        
        # post-transcriptional modifications
        messenger_rna_sequence = self.cleavage(messenger_rna_sequence)
        messenger_rna_sequence = self.polyadenylation(messenger_rna_sequence, variables)

        yield self.env.timeout(1) #TODO: implementare il tempo di trascrizione

        return messenger_rna_sequence
    
    def trascript_gene(self, dna_sequence, variables, sequence_count):
        # transcript from gene to pre-mRNA
        #messenger_rna_sequence = ''.join([BASE_COMPLEMENT_DNA2RNA[base] for base in dna_sequence])
        messenger_rna_sequence = ''

        for base in dna_sequence:
            complement_base_process = yield self.env.process(self.find_complement_base(base))
            variables.complement_base_queue_dict[sequence_count].append(complement_base_process)

            yield complement_base_process
            complement_base = complement_base_process.value
            messenger_rna_sequence += complement_base

            # wait for all the transcription gene process to be completed
            while variables.complement_base_queue_dict[sequence_count]:
                variables.complement_base_queue_dict[sequence_count].pop(0)

        return messenger_rna_sequence
    
    def find_complement_base(self, base):
        if random.random() > RNA_POLYMERASE_ERROR_RATE:
            complement_base = BASE_COMPLEMENT_DNA2RNA[base]  
        else: 
            complement_base = random.choice([b for b in list(BASE_COMPLEMENT_DNA2RNA.values()) 
                if b != BASE_COMPLEMENT_DNA2RNA[base]])

        with self.nucleotides.request(complement_base, 1) as request:
            yield request
            
            return self.env.process(self.add_complement_base(complement_base))
    
    def add_complement_base(self, complement_base):
        yield self.env.timeout(REPLICATION_TIME) # time to replicate a nucleotide
        return complement_base
    
    def splicing(self, rna_sequence):
        # remove introns: non-coding regions
        i = LENGTH_METHYL_CAP # index
        while i+3 < len(rna_sequence):
            if rna_sequence[i:i+LENGTH_EXTRON_SEQUENCE] in self.extron_sequences_list:
                i += LENGTH_EXTRON_SEQUENCE
            else: 
                # TODO: use resources Nucleotides
                rna_sequence = rna_sequence[:i] + rna_sequence[i+1:]

        return rna_sequence

    def editing(self, rna_sequence):
        for editing_site in self.editing_sites_dict.keys():
            rna_sequence = rna_sequence.replace(editing_site, self.editing_sites_dict[editing_site])
            #TODO: use resources Nucleotides
        return rna_sequence

    def capping(self, rna_sequence):
        #TODO: implement and use resources METHYL_CAP
        return 'CH3GPPP-{}'.format(rna_sequence) # Add 5'-methyl cap
    
    def cleavage(self, rna_sequence): #TODO
        # first step in adding a polyadenine tail to the pre-mRNA (post-transcriptional
        # modifications) it is necessary for producing a mature mRNA molecule
        return rna_sequence

    def polyadenylation(self, rna_sequence, variables):
        variables.poly_adenine_tail_len = random.randint(230, 270)
        """FIXME
        with self.nucleotides.request('A', poly_adenine_tail_len) as request:
            yield request
        """

        return '{}-AAAA'.format(rna_sequence) # Add PolyA tail (250 nucleotides circa)
    