import matplotlib.pyplot as plt
import ast

def barplot_proteins_number(results_df):
    plt.figure(figsize=(20, 5))
    plt.bar(results_df['number_of_proteins_synthesized'].value_counts().index,
            results_df['number_of_proteins_synthesized'].value_counts().values)
    plt.title('Number of proteins synthesized from one DNA sequence')
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

def barplot_proteins_length(results_df):
    protein_length = []
    protein_length_frequncy = []

    for lenght_list in [ast.literal_eval(x) if isinstance(x, str) else x for x in 
        results_df['length_proteins'].value_counts().index]:
        for element in lenght_list:
            protein_length.append(element)

    # agregate same number of proteins
    protein_length_final = list(set(protein_length))
    for value in protein_length_final:
        count = protein_length.count(value)
        protein_length_frequncy.append(count)
    
    plt.figure(figsize=(20, 5))
    plt.bar(protein_length_final, protein_length_frequncy);
    plt.title('Proteins length')
    plt.xlabel('Proteins length')
    plt.ylabel('Number of proteins')
    plt.show()