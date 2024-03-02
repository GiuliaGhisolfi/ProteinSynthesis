LENGTH_AMINO_CARBOXYL_GROUP = 9

def save_proteins_synthesized(dna_sequences_df, dna_sequence, mrna_sequences, polypeptides_chain, polypeptides_chain_ext,
    request_start_process_time, start_process_time, start_transcription_time, start_translation_time, 
    end_translation_time, end_process_time, promoters_box, proteins_sintetized):
    """
    Save the proteins synthesized in the dataframe of DNA sequences.

    Parameters
    ----------
    dna_sequences_df : pandas.DataFrame
        Dataframe of DNA sequences.
    dna_sequence : str
        DNA sequence.
    mrna_sequences : list
        List of mRNA sequences.
    polypeptides_chain : list
        List of polypeptides chains.
    polypeptides_chain_ext : list
        List of polypeptides chains extended.
    request_start_process_time : float
        Request start process time.
    start_process_time : float
        Start process time.
    start_transcription_time : float
        Start transcription time.
    start_translation_time : float
        Start translation time.
    end_translation_time : float
        End translation time.
    end_process_time : float
        End process time.
    promoters_box : list
        List of promoters box.
    proteins_sintetized : list
        List of proteins sintetized.

    Returns
    ----------
    pandas.DataFrame
        Dataframe of DNA sequences with the proteins synthesized.
    """
    row_index = dna_sequences_df[dna_sequences_df['sequence'] == dna_sequence].index[0]

    results = {
        'ID': dna_sequences_df.iloc[row_index]['ID'],
        'sequence': dna_sequence,
        'category': dna_sequences_df.iloc[row_index]['category'],
        'mrna_sequences': mrna_sequences,
        'length_mrna_sequences': [len(mrna) for mrna in mrna_sequences] if mrna_sequences else 0,
        'polypeptides_chains': polypeptides_chain,
        'polypeptides_chains_ext': polypeptides_chain_ext,
        'number_of_proteins_synthesized_per_mrna': proteins_sintetized,
        'number_of_proteins_synthesized': sum(proteins_sintetized) if proteins_sintetized else 0,
        'protein_synthesized': True if mrna_sequences else False, # boolean
        'request_start_process_time': request_start_process_time,
        'start_process_time': start_process_time,
        'start_transcription_time': start_transcription_time,
        'start_translation_time': start_translation_time,
        'end_translation_time': end_translation_time,
        'end_process_time': end_process_time, 
        'promoters_box': promoters_box
    }
    dna_sequences_df.iloc[row_index] = results

    return dna_sequences_df

def post_processing_results(row):
    """
    Compute post processing results.
    """
    if row['polypeptides_chains'] is not None:
        proteins = row['polypeptides_chains']
        proteins = [p for p in proteins if p != None]
        row['length_proteins'] = compute_length_proteins(proteins)
        row['number_different_proteins'] = compute_number_different_proteins(proteins)
    return row

def compute_length_proteins(proteins):
    """
    Compute the length of the proteins.
    """
    return [len(p)-LENGTH_AMINO_CARBOXYL_GROUP for p in proteins]

def compute_number_different_proteins(proteins):
    """
    Compute the number of different proteins.
    """
    return len(set(proteins))