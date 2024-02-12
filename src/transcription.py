import random

BASE_COMPLEMENT_DNA2RNA = {
    'A': 'U', 
    'T': 'A', 
    'C': 'G', 
    'G': 'C'
}
BASE_COMPLEMENT_RNA2DNA = {
    'U': 'A', 
    'A': 'T', 
    'G': 'C', 
    'C': 'G'
}
RNA_POLYMERASE_ERROR_RATE = 10e-4 # 1 error per 10^4 nucleotides
LENGTH_EXTRON_SEQUENCE = 3 # length of extron sequence
LENGTH_METHYL_CAP = 8 # length of 5'-methyl cap

class Nucleus():

    def __init__(self, extron_sequences_list, editing_sites_dict,
            promoters_sequence_list, terminator_sequence):
        self.extron_sequences_list = extron_sequences_list

        self.editing_sites_dict = editing_sites_dict
        self.editing_sites_dict = dict(sorted(self.editing_sites_dict.items(), 
            key=lambda x: len(x[0]), reverse=False)) # sort by length of key

        self.promoters_sequence_list = promoters_sequence_list
        self.terminator_sequence = terminator_sequence

        self.nucleotides = {'U': 0, 'A': 0, 'G': 0, 'C': 0} # TODO: change e implementare il conteggio quando si usano

    def transcript(self, dna_sequence): # enzime: RNA polymerase
        # detect promoter
        dna_sequence_to_transcript = self.find_promoter(dna_sequence)
        
        # transcript from gene to pre-mRNA
        messenger_rna_sequence = ''.join([BASE_COMPLEMENT_DNA2RNA[base] 
            if random.random() > RNA_POLYMERASE_ERROR_RATE 
            else random.choice([b for b in list(BASE_COMPLEMENT_DNA2RNA.values()) if b != BASE_COMPLEMENT_DNA2RNA[base]])
            for base in dna_sequence_to_transcript])

        messenger_rna_sequence = self.capping(messenger_rna_sequence)

        # elongation phase
        messenger_rna_sequence = self.splicing(messenger_rna_sequence)
        messenger_rna_sequence = self.editing(messenger_rna_sequence)
        
        # post-transcriptional modifications
        messenger_rna_sequence = self.cleavage(messenger_rna_sequence)
        messenger_rna_sequence = self.polyadenylation(messenger_rna_sequence)

        return messenger_rna_sequence # mature mRNA
    
    def find_promoter(self, dna_sequence):
        # find promoter sequence
        positions_list = [dna_sequence.find(promoter) for promoter in self.promoters_sequence_list]
        if positions_list:
            promoter_posotion = min([pos for pos in positions_list if pos > 0])
        else:
            raise ValueError('No promoter found')

        promoter_string = self.promoters_sequence_list[positions_list.index(promoter_posotion)]
        len_promoter = len(promoter_string)

        return dna_sequence[promoter_posotion+len_promoter:]
    
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
    