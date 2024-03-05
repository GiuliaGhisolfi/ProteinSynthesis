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
AMINOACIDS = {
    'Trp': 'Tryptophan', 'Cys': 'Cysteine', 'Tyr': 'Tyrosine', 'Phe': 'Phenylalanine', 'Leu': 'Leucine',
    'Ser': 'Serine', 'Pro': 'Proline', 'His': 'Histidine', 'Gln': 'Glutamine', 'Arg': 'Arginine', 
    'Ile': 'Isoleucine', 'Met': 'Methionine', 'Thr': 'Threonine', 'Asn': 'Asparagine', 'Lys': 'Lysine', 
    'Val': 'Valine', 'Ala': 'Alanine', 'Asp': 'Aspartic acid', 'Glu': 'Glutamic acid', 'Gly': 'Glycine'
}

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
    
    plt.figure(figsize=(20, 5))
    plt.scatter(length_mrna, number_of_proteins_synthesized_per_mrna, marker='.')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Number of proteins synthesized per mRNA vs mRNAlength')
    plt.xlabel('mRNA length')
    plt.ylabel('Number of proteins synthesized per mrna')
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
    plt.plot(time, number_of_proteins_synthesized, '.', alpha=0.5)
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
    plt.ylabel('Frequency')
    plt.show()

def hist_process_time(results_df):
    data_df = results_df[results_df['mrna_sequences'].notna()]
    process_time = data_df['end_process_time'] - data_df['start_process_time']
    plt.figure(figsize=(20, 5))
    plt.hist(process_time, bins=100, edgecolor='black')
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
    plt.plot(time, process_time, '.--', alpha=0.5)
    plt.title('Process time')
    plt.xlabel('Start process time (s)')
    plt.ylabel('Process time (s)')
    plt.show()

def compute_mrna_lifetime(results_df):
    end_translation_time = series_to_list(results_df[results_df['mrna_sequences'].notna()]['end_translation_time'])
    start_transcription_time = series_to_list(results_df[results_df['mrna_sequences'].notna()]['start_translation_time'])

    mrna_lifetime = np.array(end_translation_time) - np.array(start_transcription_time)

    return mrna_lifetime

def series_to_list(series):
    series = [ast.literal_eval(x) if isinstance(x, str) else x for x in series]
    series = [int(item) for sublist in series for item in sublist]
    return series

def plot_mrna_lifetime(results_df):
    mrna_lifetime = compute_mrna_lifetime(results_df)

    length_mrna = results_df[results_df['mrna_sequences'].notna()]['length_mrna_sequences']
    length_mrna = [ast.literal_eval(x) if isinstance(x, str) else x for x in length_mrna]
    length_mrna = [item for sublist in length_mrna for item in sublist]
    
    plt.figure(figsize=(20, 5))
    plt.scatter(length_mrna, mrna_lifetime, marker='.')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Mature mRNA lifetime')
    plt.xlabel('mRNA length')
    plt.ylabel('mRNA lifetime (s)')
    plt.show()

    time = results_df[results_df['mrna_sequences'].notna()]['start_process_time']
    time, mrna_lifetime = zip(*sorted(zip(time, mrna_lifetime)))
    plt.figure(figsize=(20, 5))
    plt.plot(time, mrna_lifetime, '.--', alpha=0.5)
    plt.title('Mature mRNA lifetime')
    plt.xlabel('Start process time (s)')
    plt.ylabel('mRNA lifetime (s)')
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
    plt.figure(figsize=(20, 5))
    plt.plot(rna_polymerase_df['request_time'], rna_polymerase_df['wait_time'], 
        'm.', label='RNA polymerase')
    plt.plot(rna_polymerase_df['request_time'], rna_polymerase_df['wait_time'], 
        'm--', alpha=0.5)
    plt.plot(ribosome_df['request_time'], ribosome_df['wait_time'], 
        'g.', label='Ribosome')
    plt.plot(ribosome_df['request_time'], ribosome_df['wait_time'], 
        'g--', alpha=0.5)
    plt.title('Resources request time vs wait time')
    plt.xlabel('Request time (s)')
    plt.ylabel('Wait time (s)')
    plt.legend()
    plt.show()

def plot_codons_request_per_aminoacid(file_path, time_unit=TIME_UNIT):
    codon_dict_list = []
    for codon in CODONS:
        with open(file_path+f'rna_transfer_history_{codon}.json') as f:
            codon_dict_list.append(json.load(f))
    
    with open(CODONS_PATH) as f:
        codons_dict = json.load(f)

    # map each aminoacid to a number, exclude stop codons
    aminoacids = list(set(codons_dict.values()))
    aminoacids.remove('STOP')

    plt.figure(figsize=(20, 18))
    for codon, aminoacid in codons_dict.items():
        if aminoacid != 'STOP':
            i = aminoacids.index(aminoacid)
            time_list = codon_dict_list[CODONS.index(codon)]['request_time']
            plt.subplot(5, 4, i+1)
            plt.hist(time_list, bins=50, alpha=0.3, edgecolor='black', label=codon)
            plt.title(AMINOACIDS[aminoacid])
            plt.xlabel('Request time (s)')
            if i%4 == 0:
                plt.ylabel('Number of requests')
            plt.subplots_adjust(hspace=0.3)
            plt.legend()
    plt.show()

def plot_codons_request(file_path, time_unit=TIME_UNIT):
    codon_dict_list = []
    for codon in CODONS:
        with open(file_path+f'rna_transfer_history_{codon}.json') as f:
            codon_dict_list.append(json.load(f))
    
    max_time = max([max(codon_dict['request_time']) for codon_dict in codon_dict_list])
    time = np.arange(0, max_time, time_unit)

    with open(CODONS_PATH) as f:
        codons_dict = json.load(f)

    # map each aminoacid to a number, exclude stop codons
    aminoacids = list(set(codons_dict.values()))
    #aminoacids.remove('STOP')

    colors_list = ['b', 'm', 'g', 'r', 'c', 'y', 'k', 'orange', 'purple', 
        'brown', 'pink', 'lightblue', 'lightseagreen',
        'lightgreen', 'gray', 'lightcoral', 'lightcyan', 'lightyellow', 
        'lightgrey', 'lightpink', 'lightsteelblue', ]

    plt.figure(figsize=(20, 7))
    for codon_dict, codon in zip(codon_dict_list, CODONS):
        aminoacid = codons_dict[codon]
        i = aminoacids.index(aminoacid)

        requestes = requestes_serie_over_time(codon_dict, time_unit)
        requestes.extend([0] * (len(time) - len(requestes)))

        plt.plot(time, requestes, '.', color=colors_list[i], alpha=0.5, label=codon)
    plt.title('Number of requests of tRNA')
    plt.ylabel('Number of requests')
    plt.xlabel('Request time (s)')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=3)
    plt.show()

def requestes_serie_over_time(codons_dict, time_unit=TIME_UNIT):
    # return a serie of requests over time
    time_list = codons_dict['request_time']

    requestes = [0]
    current_time = 0
    req = 0
    for time in time_list:
        time_steps = int(time / time_unit) - int(current_time / time_unit)
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

def display_model_parameters(parameters_dict_list):
    parameters_df = create_model_df(parameters_dict_list)
    print('Parameters shared by all models:\n'
        f'simulation_time: {parameters_df["simulation_time"].unique()[0]}\n'
        f'number_resources: {parameters_df["number_resources"].unique()[0]}\n'
        f'number_rna_transfers_per_codon: {parameters_df["number_rna_transfers_per_codon"].unique()[0]}\n'
        f'uracil_initial_amount: {parameters_df["uracil_initial_amount"].unique()[0]}\n'
        f'adenine_initial_amount: {parameters_df["adenine_initial_amount"].unique()[0]}\n'
        f'guanine_initial_amount: {parameters_df["guanine_initial_amount"].unique()[0]}\n'
        f'cytosine_initial_amount: {parameters_df["cytosine_initial_amount"].unique()[0]}\n')

    print('Paramethers that differ between models:')
    print(parameters_df[['number_rna_polymerases', 'number_ribosomes']])

def compare_wait_time(df_list, df_list_seed, parameters_dict_list, resource_name):
    plt.subplots(3, 2, figsize=(20, 10))

    for i in range(len(df_list)):
        plt.subplot(3, 2, i+1)
        length = len(df_list[i]['wait_time'])
        plt.plot(df_list[i]['request_time'][:length], df_list[i]['wait_time'], '.--', alpha=0.5, label='Original')
        length = len(df_list_seed[i]['wait_time'])
        plt.plot(df_list_seed[i]['request_time'][:length], df_list_seed[i]['wait_time'], '.--', alpha=0.5, label='Seed')
        plt.xlabel('Request time')
        plt.ylabel('Wait time')
        plt.title(f'Test {i}: ribosome {parameters_dict_list[i]["number_ribosomes"]}, rna_polymerase {parameters_dict_list[i]["number_rna_polymerases"]}')
        plt.legend()
    plt.suptitle('Comparison of wait time for '+ resource_name)
    plt.tight_layout()
    plt.show()

def compare_proteins_number_over_time(results_df_list, df_list_seed, parameters_dict_list):
    plt.subplots(3, 2, figsize=(20, 10))
    for i in range(len(results_df_list)):
        number_of_proteins_synthesized = results_df_list[i][results_df_list[i]['mrna_sequences'
            ].notna()]['number_of_proteins_synthesized']
        time = results_df_list[i][results_df_list[i]['mrna_sequences'].notna()]['end_process_time']
        time, number_of_proteins_synthesized = zip(*sorted(zip(time, number_of_proteins_synthesized)))

        number_of_proteins_synthesized_seed = df_list_seed[i][df_list_seed[i]['mrna_sequences'
            ].notna()]['number_of_proteins_synthesized']
        time_seed = df_list_seed[i][df_list_seed[i]['mrna_sequences'].notna()]['end_process_time']
        time_seed, number_of_proteins_synthesized_seed = zip(*sorted(zip(time_seed, number_of_proteins_synthesized_seed)))

        plt.subplot(3, 2, i+1)
        plt.plot(time, number_of_proteins_synthesized, '.--', alpha=0.5, label='Original')
        plt.plot(time_seed, number_of_proteins_synthesized_seed, '.--', alpha=0.5, label='Seed')
        plt.xlabel('End process time')
        plt.ylabel('Number of proteins')
        plt.title(f'Test {i}: ribosome {parameters_dict_list[i]["number_ribosomes"]}, rna_polymerase {parameters_dict_list[i]["number_rna_polymerases"]}')
        plt.legend()
    plt.suptitle('Comparison of number of proteins synthesized over time')
    plt.tight_layout()
    plt.show()

def compare_number_proteins_per_length_mrna(results_df_list, df_list_seed, parameters_dict_list):
    plt.subplots(3, 2, figsize=(20, 10))
    for i in range(len(results_df_list)):
        length_mrna = results_df_list[i][results_df_list[i]['mrna_sequences'].notna()]['length_mrna_sequences']
        length_mrna = [ast.literal_eval(x) if isinstance(x, str) else x for x in length_mrna]
        length_mrna = [int(item) for sublist in length_mrna for item in sublist]

        number_of_proteins_synthesized_per_mrna = results_df_list[i][results_df_list[i]['mrna_sequences'
            ].notna()]['number_of_proteins_synthesized_per_mrna']
        number_of_proteins_synthesized_per_mrna = [ast.literal_eval(x) if isinstance(x, str) 
            else x for x in number_of_proteins_synthesized_per_mrna]
        number_of_proteins_synthesized_per_mrna = [int(item) for sublist in 
            number_of_proteins_synthesized_per_mrna for item in sublist]

        length_mrna_seed = df_list_seed[i][df_list_seed[i]['mrna_sequences'].notna()]['length_mrna_sequences']
        length_mrna_seed = [ast.literal_eval(x) if isinstance(x, str) else x for x in length_mrna_seed]
        length_mrna_seed = [item for sublist in length_mrna_seed for item in sublist]

        number_of_proteins_synthesized_per_mrna_seed = df_list_seed[i][df_list_seed[i]['mrna_sequences'
            ].notna()]['number_of_proteins_synthesized_per_mrna']
        number_of_proteins_synthesized_per_mrna_seed = [ast.literal_eval(x) if isinstance(x, str) 
            else x for x in number_of_proteins_synthesized_per_mrna_seed]
        number_of_proteins_synthesized_per_mrna_seed = [int(item) for sublist in 
            number_of_proteins_synthesized_per_mrna_seed for item in sublist]

        plt.subplot(3, 2, i+1)
        plt.scatter(length_mrna, number_of_proteins_synthesized_per_mrna, marker='.', label='Original')
        plt.scatter(length_mrna_seed, number_of_proteins_synthesized_per_mrna_seed, marker='.', label='Seed')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('mRNA length')
        plt.ylabel('Number of proteins synthesized per mRNA')
        plt.title(f'Test {i}: ribosome {parameters_dict_list[i]["number_ribosomes"]}, rna_polymerase {parameters_dict_list[i]["number_rna_polymerases"]}')
        plt.legend()
    plt.suptitle('Comparison of number of proteins synthesized per mRNA vs mRNA length')
    plt.tight_layout()
    plt.show()

def compare_process_time(results_df_list, df_list_seed, parameters_dict_list):
    plt.subplots(3, 2, figsize=(20, 10))
    for i in range(len(results_df_list)):
        data_df = results_df_list[i][results_df_list[i]['mrna_sequences'].notna()]
        process_time = data_df['end_process_time'] - data_df['start_process_time']
        time = data_df['start_process_time']
        time, process_time = zip(*sorted(zip(time, process_time)))
        plt.subplot(3, 2, i+1)
        plt.plot(time, process_time, '.--', alpha=0.5, label='Original')
        data_df_seed = df_list_seed[i][df_list_seed[i]['mrna_sequences'].notna()]
        process_time_seed = data_df_seed['end_process_time'] - data_df_seed['start_process_time']
        time_seed = data_df_seed['start_process_time']
        time_seed, process_time_seed = zip(*sorted(zip(time_seed, process_time_seed)))
        plt.plot(time_seed, process_time_seed, '.--', alpha=0.5, label='Seed')
        plt.xlabel('Start process time')
        plt.ylabel('Process time')
        plt.title(f'Test {i}: ribosome {parameters_dict_list[i]["number_ribosomes"]}, rna_polymerase {parameters_dict_list[i]["number_rna_polymerases"]}')
        plt.legend()
    plt.suptitle('Comparison of process time')
    plt.tight_layout()
    plt.show()

def compare_mrna_lifetime(results_df_list, results_df_list_seed, parameters_dict_list, scale='not_log'):
    plt.subplots(3, 2, figsize=(20, 10))
    for i in range(len(results_df_list)):
        mrna_lifetime = compute_mrna_lifetime(results_df_list[i])
        length_mrna = results_df_list[i][results_df_list[i]['mrna_sequences'].notna()]['length_mrna_sequences']
        length_mrna = [ast.literal_eval(x) if isinstance(x, str) else x for x in length_mrna]
        length_mrna = [item for sublist in length_mrna for item in sublist]

        mrna_lifetime_seed = compute_mrna_lifetime(results_df_list_seed[i])
        length_mrna_seed = results_df_list_seed[i][results_df_list_seed[i]['mrna_sequences'
            ].notna()]['length_mrna_sequences']
        length_mrna_seed = [ast.literal_eval(x) if isinstance(x, str) else x for x in length_mrna_seed]
        length_mrna_seed = [item for sublist in length_mrna_seed for item in sublist]

        plt.subplot(3, 2, i+1)
        plt.scatter(length_mrna, mrna_lifetime, marker='.', label='Original')
        plt.scatter(length_mrna_seed, mrna_lifetime_seed, marker='.', label='Seed')
        if scale == 'log':
            plt.xscale('log')
            plt.yscale('log')
        plt.xlabel('mRNA length')
        plt.ylabel('mRNA lifetime (s)')
        plt.title(f'Test {i}: ribosome {parameters_dict_list[i]["number_ribosomes"]}, rna_polymerase {parameters_dict_list[i]["number_rna_polymerases"]}')
        plt.legend()
    plt.suptitle('Comparison of mature mRNA lifetime')
    plt.tight_layout()
    plt.show()