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

def splicing(rna_sequence, intron_sequences):
    for intron_sequence in intron_sequences:
        rna_sequence = rna_sequence.replace(intron_sequence, '')
    return rna_sequence

def editing(rna_sequence, editing_sites_dict):
    for editing_site in editing_sites_dict:
        rna_sequence = rna_sequence.replace(editing_site, editing_sites_dict[editing_site])
    return rna_sequence

def capping(rna_sequence):
    return 'CH3GPPP-{}'.format(rna_sequence) # Add 5'-methyl cap

def polyadenylation(rna_sequence):
    return '{}-AAAA'.format(rna_sequence) # Add PolyA tail

class Nucleus():

    def __init__(self, intron_sequences_list, editing_sites_dict,
            promoters_sequence_list, terminator_sequence):
        self.intron_sequences_list = intron_sequences_list

        self.editing_sites_dict = editing_sites_dict
        self.editing_sites_dict = dict(sorted(self.editing_sites_dict.items(), 
            key=lambda x: len(x[0]), reverse=False)) # sort by length of key

        self.promoters_sequence_list = promoters_sequence_list
        self.terminator_sequence = terminator_sequence
        # enzime: RNA polymerase

        # type of RNA
        # - messenger RNA (mRNA): code for proteins
        # - transfer RNA (tRNA): transport amino acids to ribosomes
        # - ribosomal RNA (rRNA): form the core of a ribosome's functional center
        # - small nuclear RNA (snRNA): splicing of pre-mRNA
        # - small cajan RNA (scaRNA): modification of snRNA and snoRNA
        # - small nucleolar RNA (snoRNA): processing and modification of rRNA, tRNA, and snRNA
        # - microRNA (miRNA): regulate gene expression
        # - small interfering RNA (siRNA): regulate gene expression
    
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
    
    def find_terminator(self, rna_sequence):
        # find terminator sequence
        positions_list = [pos for pos in [rna_sequence.find(terminator) for terminator in 
            self.terminator_sequence] if pos > 0]
        terminator_position = -1 if not positions_list else max(positions_list) 
        
        return rna_sequence[:terminator_position]

    def transcript(self, dna_sequence):
        # detect promoter
        dna_sequence_to_transcript = self.find_promoter(dna_sequence)
        
        # transcription initiation
        # add error
        messenger_rna_sequence = ''.join([BASE_COMPLEMENT_DNA2RNA[base] 
            if random.random() > RNA_POLYMERASE_ERROR_RATE 
            else random.choice([b for b in list(BASE_COMPLEMENT_DNA2RNA.values()) if b != BASE_COMPLEMENT_DNA2RNA[base]])
            for base in dna_sequence_to_transcript])
        messenger_rna_sequence = self.find_terminator(messenger_rna_sequence)

        # transcription elongation
        messenger_rna_sequence = splicing(messenger_rna_sequence, self.intron_sequences_list)
        messenger_rna_sequence = editing(messenger_rna_sequence, self.editing_sites_dict)
        messenger_rna_sequence = capping(messenger_rna_sequence)
        messenger_rna_sequence = polyadenylation(messenger_rna_sequence)

        return messenger_rna_sequence
    
    