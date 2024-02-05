def transcription(dna_sequence):
    rna_base_complement_dict = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
    rna_sequence = [rna_base_complement_dict[base] for base in dna_sequence]
    rna_sequence = ''.join(rna_sequence)

    intron_sequences = ['GU', 'AG'] # TODO: check and add more intron sequences
    rna_sequence = splicing(rna_sequence, intron_sequences)

    editing_sites_dict = {'A': 'I', 'C': 'U', 'G': 'U'} # TODO: check and add more editing sites
    rna_sequence = editing(rna_sequence, editing_sites_dict)

    rna_sequence = capping(rna_sequence)
    rna_sequence = polyadenylation(rna_sequence)

    return rna_sequence

def splicing(rna_sequence, intron_sequences):
    rna_sequence = [rna_sequence.replace(intron_sequence, '') for intron_sequence in intron_sequences]

    return ''.join(rna_sequence)

def editing(rna_sequence, editing_sites_dict):
    rna_sequence = [rna_sequence.replace(key, value) for key, value in editing_sites_dict.items()]

    return ''.join(rna_sequence)

def capping(rna_sequence): # TODO: change to add 5' cap
    return '5\'-{}-3\''.format(rna_sequence) # Add 5' cap

def polyadenylation(rna_sequence):
    return '{}-3\'-PolyA'.format(rna_sequence) # Add PolyA tail

def translation(rna_sequence):
    pass

def initialization():
    pass

def elongation():
    pass

def termination():
    pass