
def save_proteins_synthesized(dna_sequences_df, dna_sequence, mrna_sequences, polypeptides_chain, polypeptides_chain_ext,
    request_start_process_time, start_process_time, start_transcription_time, start_translation_time, 
    end_translation_time, end_process_time, promoters_box):
    row_index = dna_sequences_df[dna_sequences_df['sequence'] == dna_sequence].index[0]

    results = {
        'ID': dna_sequences_df.iloc[row_index]['ID'],
        'sequence': dna_sequence,
        'category': dna_sequences_df.iloc[row_index]['category'],
        'mrna_sequences': mrna_sequences,
        'polypeptides_chains': polypeptides_chain,
        'polypeptides_chains_ext': polypeptides_chain_ext,
        'number_of_proteins_synthesized': len(mrna_sequences) if mrna_sequences else 0,
        'protein_synthesized': True if mrna_sequences else False,
        'request_start_process_time': request_start_process_time,
        'start_process_time': start_process_time,
        'start_transcription_time': start_transcription_time,
        'start_translation_time': start_translation_time,
        'end_translation_time': end_translation_time,
        'end_process_time': end_process_time, 
        'promoters_box': promoters_box,
    }
    dna_sequences_df.iloc[row_index] = results

    return dna_sequences_df

def post_processing_results(row):
    if row['polypeptides_chains'] is not None:
        proteins = row['polypeptides_chains']
        proteins = [p for p in proteins if p != None]
        row['length_proteins'] = compute_length_proteins(proteins)
        row['number_different_proteins'] = compute_number_different_proteins(proteins)
    return row

def compute_length_proteins(proteins):
    return [len(p)-9 for p in proteins]

def compute_number_different_proteins(proteins):
    return len(set(proteins))