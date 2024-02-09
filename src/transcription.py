def splicing(rna_sequence, intron_sequences):
    rna_sequence = [rna_sequence.replace(intron_sequence, '') for intron_sequence in intron_sequences]

    return ''.join(rna_sequence)

def editing(rna_sequence, editing_sites_dict):
    rna_sequence = [rna_sequence.replace(key, value) for key, value in editing_sites_dict.items()]

    return ''.join(rna_sequence)

def capping(rna_sequence):
    return 'CH3GPPP{}'.format(rna_sequence) # Add 5'-methyl cap

def polyadenylation(rna_sequence):
    return '{}AAAA'.format(rna_sequence) # Add PolyA tail

BASE_COMPLEMENT = {
    'A': 'U', 
    'T': 'A', 
    'C': 'G', 
    'G': 'C'
}

class Transcription():

    def __init__(self, intron_sequences_list, editing_sites_dict):
        self.intron_sequences_list = intron_sequences_list
        self.editing_sites_dict = editing_sites_dict
    
    def transcript(self, dna_sequence):
        rna_sequence = [BASE_COMPLEMENT[base] for base in dna_sequence]
        rna_sequence = ''.join(rna_sequence)

        rna_sequence = splicing(rna_sequence, self.intron_sequences_list)

        rna_sequence = editing(rna_sequence, self.editing_sites_dict)

        rna_sequence = capping(rna_sequence)
        rna_sequence = polyadenylation(rna_sequence)

        return rna_sequence
    
    