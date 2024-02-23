
def save_proteins_synthesized(dna_sequences_df, dna_sequence, mrna_sequences, polypeptides_chain, polypeptides_chain_ext,
    request_start_process_time, start_process_time, start_transcription_time, start_translation_time, 
    end_translation_time, end_process_time):
    # TODO: gestire errori (dna_sequence non trovata nel dataframe)
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
        'end_process_time': end_process_time
    }
    dna_sequences_df.iloc[row_index] = results

    return dna_sequences_df