import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import ast
import json

TIME_UNIT = 0.0001
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

def level_series_over_time(nucleotide_dict, time_unit=TIME_UNIT):
    levels_list = nucleotide_dict['level']
    time_list = nucleotide_dict['time']

    levels = [levels_list[0]]
    current_time = 0
    for level, time in zip(levels_list[1:], time_list[1:]):
        delta_t = np.round(time - current_time, 4)
        time_steps = int(delta_t / time_unit)
        levels.extend([level] * time_steps)
        current_time = time

    return levels

def barplot_nucleotide_level_over_time(
        uracil_dict, adenine_dict, guanine_dict, cytosine_dict, time_unit=TIME_UNIT):
    uracil_levels = level_series_over_time(uracil_dict, time_unit)
    adenine_levels = level_series_over_time(adenine_dict, time_unit)
    guanine_levels = level_series_over_time(guanine_dict, time_unit)
    cytosine_levels = level_series_over_time(cytosine_dict, time_unit)

    len_min = min(len(uracil_levels), len(adenine_levels), len(guanine_levels), len(cytosine_levels))

    df = pd.DataFrame({
        'time': np.arange(0, len_min*time_unit, time_unit),
        'Uracil': uracil_levels[:len_min],
        'Adenine': adenine_levels[:len_min],
        'Guanine': guanine_levels[:len_min],
        'Cytosine': cytosine_levels[:len_min]
        })
    
    fig = px.bar(
        df, 
        y=['Uracil', 'Adenine', 'Guanine', 'Cytosine'], #FIXME
        animation_frame='time',
        title='Nucleotides levels over time',
        labels={'value': 'Nucleotides level', 'index': 'Time'},
        template='plotly_white')
    fig.show()

def plot_nucleotide_level_over_time(
        uracil_dict, adenine_dict, guanine_dict, cytosine_dict, time_unit=TIME_UNIT):
    uracil_levels = level_series_over_time(uracil_dict)
    adenine_levels = level_series_over_time(adenine_dict)
    guanine_levels = level_series_over_time(guanine_dict)
    cytosine_levels = level_series_over_time(cytosine_dict)

    len_min = min(len(uracil_levels), len(adenine_levels), len(guanine_levels), len(cytosine_levels))

    plt.figure(figsize=(20, 5))
    plt.plot(np.arange(0, len_min*time_unit, time_unit), uracil_levels[:len_min], label='Uracil')
    plt.plot(np.arange(0, len_min*time_unit, time_unit), adenine_levels[:len_min], label='Adenine')
    plt.plot(np.arange(0, len_min*time_unit, time_unit), guanine_levels[:len_min], label='Guanine')
    plt.plot(np.arange(0, len_min*time_unit, time_unit), cytosine_levels[:len_min], label='Cytosine')
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

def plot_codons_request(file_path, len_min, time_unit=TIME_UNIT):
    codon_dict_list = []
    for codon in CODONS:
        with open(f'results/test/rna_transfer/rna_transfer_history_{codon}.json') as f:
            codon_dict_list.append(json.load(f))
    
    len_min = min([len(codon_dict['request_time']) for codon_dict in codon_dict_list])

    plt.figure(figsize=(20, 5))
    for codon_dict, codon in zip(codon_dict_list, CODONS):
        requestes = requestes_serie_over_time(codon_dict, time_unit)
        plt.plot(np.arange(0, len_min*time_unit, time_unit), requestes, label=codon)
    plt.title('Number of requests of tRNA')
    plt.xlabel('Number of requests')
    plt.ylabel('Request time (s)')
    plt.legend()
    plt.show()

def requestes_serie_over_time(codons_dict, len_min, time_unit=TIME_UNIT):
    time_list = codons_dict['request_time']

    requestes = []
    current_time = 0
    for time in time_list[:len_min]:
        delta_t = np.round(time - current_time, 4)
        time_steps = int(delta_t / time_unit)
        requestes.extend([1] * time_steps)
        current_time = time

    return requestes