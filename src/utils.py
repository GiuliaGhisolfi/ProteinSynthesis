import matplotlib.pyplot as plt

def save_proteins_synthesized(dna_sequences_df, dna_sequence, mrna_sequences, polypeptides_chain, polypeptides_chain_ext,
    request_start_process_time, start_process_time, start_transcription_time, start_translation_time, 
    end_translation_time, end_process_time):

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

def barplot_proteins_number(results_df):
    plt.figure(figsize=(20, 5))
    plt.bar(results_df['number_of_proteins_synthesized'].value_counts().index,
            results_df['number_of_proteins_synthesized'].value_counts().values)
    plt.title('Number of proteins synthesized')
    plt.xlabel('Number of proteins')
    plt.ylabel('Number of DNA sequences')
    plt.show()

def plot_cumulative_proteins_number_over_time(results_df):
    number_of_proteins_synthesized = results_df[results_df['mrna_sequences'].notna()]['number_of_proteins_synthesized']
    time = results_df[results_df['mrna_sequences'].notna()]['end_process_time']

    time, number_of_proteins_synthesized = zip(*sorted(zip(time, number_of_proteins_synthesized)))
    cumulative_number_of_proteins_synthesized = [number_of_proteins_synthesized[0]]
    for i in range(len(number_of_proteins_synthesized)-1):
        cumulative_number_of_proteins_synthesized.append(
            cumulative_number_of_proteins_synthesized[i]+
            number_of_proteins_synthesized[i+1])

    plt.figure(figsize=(20, 5))
    plt.plot(time, cumulative_number_of_proteins_synthesized)
    plt.title('Cumulative number of proteins synthesized over time')
    plt.xlabel('Time')
    plt.ylabel('Number of proteins')
    plt.show()

def plot_proteins_number_over_time(results_df):
    number_of_proteins_synthesized = results_df[results_df['mrna_sequences'].notna()]['number_of_proteins_synthesized']
    time = results_df[results_df['mrna_sequences'].notna()]['end_process_time']

    plt.figure(figsize=(20, 5))
    plt.plot(time, number_of_proteins_synthesized, '.')
    plt.title('Number of proteins synthesized over time')
    plt.xlabel('Time')
    plt.ylabel('Number of proteins')
    plt.show()