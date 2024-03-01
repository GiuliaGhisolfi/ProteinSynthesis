import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import ast
import json

TIME_UNIT = 0.0001
CODONS_PATH = 'data\codons.json'
CODONS = [
    'UUU', 'UUC', 'UUA', 'UUG', 'UCU', 'UCC', 'UCA', 'UCG', 'UAU', 'UAC', 'UAA', 'UAG', 'UGU', 'UGC', 'UGA', 'UGG',
    'CUU', 'CUC', 'CUA', 'CUG', 'CCU', 'CCC', 'CCA', 'CCG','CAU', 'CAC', 'CAA', 'CAG', 'CGU', 'CGC', 'CGA', 'CGG',
    'AUU', 'AUC', 'AUA', 'AUG', 'AGU', 'AGC', 'AGA', 'AGG', 'AAU', 'AAC', 'AAA', 'AAG', 'ACU', 'ACC', 'ACA', 'ACG',
    'GUU', 'GUC', 'GUA', 'GUG', 'GCU', 'GCC', 'GCA', 'GCG', 'GAU', 'GAC', 'GAA', 'GAG', 'GGU', 'GGC', 'GGA', 'GGG'
    ]

def barplot_proteins_number(results_df):
    plt.figure(figsize=(20, 5))
    plt.bar(results_df['number_of_proteins_synthesized'].value_counts().index,
            results_df['number_of_proteins_synthesized'].value_counts().values)
    plt.title('Number of proteins synthesized from one DNA sequence')
    plt.xlabel('Number of proteins')
    plt.ylabel('Number of DNA sequences')
    plt.show()

def barplot_number_proteins_per_mrna(results_df):
    number_of_proteins_synthesized_per_mrna = results_df[
        results_df['mrna_sequences'].notna()]['number_of_proteins_synthesized_per_mrna']
    number_of_proteins_synthesized_per_mrna = [ast.literal_eval(x) if isinstance(x, str) 
        else x for x in number_of_proteins_synthesized_per_mrna]
    number_of_proteins_synthesized_per_mrna = [int(item) for sublist in 
        number_of_proteins_synthesized_per_mrna for item in sublist]
    
    plt.figure(figsize=(20, 5))
    plt.bar(*np.unique(number_of_proteins_synthesized_per_mrna, return_counts=True))
    plt.title('Number of proteins synthesized per mRNA')
    plt.xlabel('Number of proteins synthesized per mRNA')
    plt.ylabel('Frequency')
    plt.show()

def plot_number_proteins_per_length_mrna(results_df):
    length_mrna = results_df[results_df['mrna_sequences'].notna()]['length_mrna_sequences']
    length_mrna = [ast.literal_eval(x) if isinstance(x, str) else x for x in length_mrna]
    length_mrna = [int(item) for sublist in length_mrna for item in sublist]

    number_of_proteins_synthesized_per_mrna = results_df[
        results_df['mrna_sequences'].notna()]['number_of_proteins_synthesized_per_mrna']
    number_of_proteins_synthesized_per_mrna = [ast.literal_eval(x) if isinstance(x, str) 
        else x for x in number_of_proteins_synthesized_per_mrna]
    number_of_proteins_synthesized_per_mrna = [int(item) for sublist in 
        number_of_proteins_synthesized_per_mrna for item in sublist]
    
    # sort the lists
    length_mrna, number_of_proteins_synthesized_per_mrna = zip(*sorted(
        zip(length_mrna, number_of_proteins_synthesized_per_mrna)))
    
    plt.figure(figsize=(20, 5))
    plt.plot(length_mrna, number_of_proteins_synthesized_per_mrna, '.--')
    plt.title('Number of proteins synthesized per mRNA length')
    plt.xlabel('mRNA length')
    plt.ylabel('Number of proteins')
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
    plt.xlabel('Time (s)')
    plt.ylabel('Number of proteins')
    plt.show()

def plot_proteins_number_over_time(results_df):
    number_of_proteins_synthesized = results_df[results_df['mrna_sequences'].notna()]['number_of_proteins_synthesized']
    time = results_df[results_df['mrna_sequences'].notna()]['end_process_time']
    time, number_of_proteins_synthesized = zip(*sorted(zip(time, number_of_proteins_synthesized)))

    plt.figure(figsize=(20, 5))
    plt.plot(time, number_of_proteins_synthesized, '.--')
    plt.title('Number of proteins synthesized over time')
    plt.xlabel('Time (s)')
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
    plt.bar(protein_length_final, protein_length_frequncy)
    plt.title('Proteins length')
    plt.xlabel('Proteins length')
    plt.ylabel('Number of proteins')
    plt.show()

def hist_process_time(results_df):
    data_df = results_df[results_df['mrna_sequences'].notna()]
    process_time = data_df['end_process_time'] - data_df['start_process_time']
    plt.figure(figsize=(20, 5))
    plt.hist(process_time, bins=100)
    plt.title('Process time')
    plt.xlabel('Process time (s)')
    plt.ylabel('Number of DNA sequences')
    plt.show()

def plot_process_time(results_df):
    data_df = results_df[results_df['mrna_sequences'].notna()]
    process_time = data_df['end_process_time'] - data_df['start_process_time']
    time = data_df['start_process_time']
    time, process_time = zip(*sorted(zip(time, process_time)))
    plt.figure(figsize=(20, 5))
    plt.plot(time, process_time, '.--')
    plt.title('Process time')
    plt.xlabel('Start process time (s)')
    plt.ylabel('Process time (s)')
    plt.show()

def level_series_over_time(nucleotide_dict, time_unit=TIME_UNIT):
    levels_list = nucleotide_dict['level']
    time_list = nucleotide_dict['time']

    levels = [levels_list[0]]
    current_time = 0
    level = levels_list[0]
    for level, time in zip(levels_list[1:], time_list[1:]):
        delta_t = np.round(time - current_time, 4)
        time_steps = int(delta_t / time_unit)
        levels.extend([level] * time_steps)
        if time_steps > 0:
            current_time = time

    return levels

def plot_nucleotide_level_over_time(
        uracil_dict, adenine_dict, guanine_dict, cytosine_dict, time_unit=TIME_UNIT):
    uracil_levels = level_series_over_time(uracil_dict, time_unit)
    adenine_levels = level_series_over_time(adenine_dict, time_unit)
    guanine_levels = level_series_over_time(guanine_dict, time_unit)
    cytosine_levels = level_series_over_time(cytosine_dict, time_unit)

    max_time = int(max([max(uracil_dict['time']), max(adenine_dict['time']), 
        max(guanine_dict['time']), max(cytosine_dict['time'])]))
    time = np.arange(0, max_time, time_unit)
    uracil_levels.extend([uracil_levels[-1]] * (len(time) - len(uracil_levels)))
    adenine_levels.extend([adenine_levels[-1]] * (len(time) - len(adenine_levels)))
    guanine_levels.extend([guanine_levels[-1]] * (len(time) - len(guanine_levels)))
    cytosine_levels.extend([cytosine_levels[-1]] * (len(time) - len(cytosine_levels)))

    plt.figure(figsize=(20, 5))
    plt.plot(time, uracil_levels, label='Uracil')
    plt.plot(time, adenine_levels, label='Adenine')
    plt.plot(time, guanine_levels, label='Guanine')
    plt.plot(time, cytosine_levels, label='Cytosine')
    plt.title('Nucleotides levels over time')
    plt.xlabel('Time (s)')
    plt.ylabel('Nucleotides level')
    plt.legend()
    plt.show()

def dict_to_dataframe(resources_dict):
    min_len = min([len(v) for v in resources_dict.values()])
    resources_dict = {k: v[:min_len] for k, v in resources_dict.items()}
    df = pd.DataFrame(resources_dict)
    return df

def resources_request_wait_time(rna_polymerase_df, ribosome_df):
    plt.figure(figsize=(10, 5))
    plt.plot(rna_polymerase_df['request_time'], rna_polymerase_df['wait_time'], label='RNA polymerase')
    plt.plot(ribosome_df['request_time'], ribosome_df['wait_time'], label='Ribosome')
    plt.title('Resources request time vs wait time')
    plt.xlabel('Request time (s)')
    plt.ylabel('Wait time (s)')
    plt.legend()
    plt.show()

def plot_codons_request(file_path, time_unit=TIME_UNIT):
    codon_dict_list = []
    for codon in CODONS:
        with open(file_path+f'rna_transfer_history_{codon}.json') as f:
            codon_dict_list.append(json.load(f))
    
    max_time = max([max(codon_dict['request_time']) for codon_dict in codon_dict_list])
    time = np.arange(0, max_time, time_unit)

    plt.figure(figsize=(18, 10))
    for codon_dict, codon in zip(codon_dict_list, CODONS):
        requestes = requestes_serie_over_time(codon_dict, time_unit)
        requestes.extend([0] * (len(time) - len(requestes)))
        plt.plot(time, requestes, '.--', alpha=0.5, label=codon)
    plt.title('Number of requests of tRNA')
    plt.ylabel('Number of requests')
    plt.xlabel('Request time (s)')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
    plt.show()

def plot_codons_request_per_aminoacid(file_path, time_unit=TIME_UNIT):
    codon_dict_list = []
    for codon in CODONS:
        with open(file_path+f'rna_transfer_history_{codon}.json') as f:
            codon_dict_list.append(json.load(f))
    
    with open(CODONS_PATH) as f:
        codons_dict = json.load(f)
    
    max_time = max([max(codon_dict['request_time']) for codon_dict in codon_dict_list])
    time = np.arange(0, max_time, time_unit)

    plt.subplots(4, 5, figsize=(20, 10))
    # subplot for each aminoacid (codons_dict keys: codon, value: aminoacid)
    for i, (aminoacid, codons) in enumerate(codons_dict.items()):
        for codon in codons:
            print(codon) # FIXME
            requestes = requestes_serie_over_time(codon_dict_list[CODONS.index(codon)], time_unit)
            requestes.extend([0] * (len(time) - len(requestes)))
            plt.subplot(4, 5, i+1)
            plt.plot(time, requestes, '.--', alpha=0.5, label=codon)
            plt.title(aminoacid)
            plt.ylabel('Number of requests')
            plt.xlabel('Request time (s)')
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)
    plt.show()

def requestes_serie_over_time(codons_dict, time_unit=TIME_UNIT):
    time_list = codons_dict['request_time']

    requestes = [0]
    current_time = 0
    req = 0
    for time in time_list:
        delta_t = np.round(time - current_time, 4)
        time_steps = int(delta_t / time_unit)
        if time_steps > 0:
            requestes.extend([req] * time_steps)
            req = 0
        else:
            req += 1
        current_time = time

    return requestes

##################### Compare models #####################
def create_model_df(parameters_dict_list):
    df = pd.DataFrame(parameters_dict_list)
    return(df)

def compare_wait_time(df_list):
    plt.figure(figsize=(20, 5))
    for i, df in enumerate(df_list):
        plt.plot(df['request_time'], df['wait_time'], '.--', label=f'Model {i}')
    plt.title('Resources request time vs wait time')
    plt.xlabel('Request time (s)')
    plt.ylabel('Wait time (s)')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()

def compare_proteins_number_over_time(results_df_list):
    plt.figure(figsize=(20, 5))
    for i, results_df in enumerate(results_df_list):
        number_of_proteins_synthesized = results_df[results_df['mrna_sequences'].notna()]['number_of_proteins_synthesized']
        time = results_df[results_df['mrna_sequences'].notna()]['end_process_time']
        time, number_of_proteins_synthesized = zip(*sorted(zip(time, number_of_proteins_synthesized)))
        plt.plot(time, number_of_proteins_synthesized, '.--', label=f'Model {i}')
    plt.title('Number of proteins synthesized over time')
    plt.xlabel('Time (s)')
    plt.ylabel('Number of proteins')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()

def compare_process_time(results_df_list):
    plt.figure(figsize=(20, 5))
    for i, results_df in enumerate(results_df_list):
        data_df = results_df[results_df['mrna_sequences'].notna()]
        process_time = data_df['end_process_time'] - data_df['start_process_time']
        time = data_df['start_process_time']
        time, process_time = zip(*sorted(zip(time, process_time)))
        plt.plot(time, process_time, '.--', label=f'Model {i}')
    plt.title('Process time')
    plt.xlabel('Start process time (s)')
    plt.ylabel('Process time (s)')
    plt.show()