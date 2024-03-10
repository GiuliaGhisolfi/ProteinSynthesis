# Protein Synthesis in Eukaryotic Cells 🧬🧫
This repository contains the code for a project aimed at modeling protein synthesis in eukaryotic cells as a discrete-time process.

## Modeling
The process is simulated using `SimPy` framework, with customized classes tailored to the specific requirements. 

The process is modeled by the `EukaryoticCell` class, which coordinates transcription and translation through its components: the `Nucleus` class handles transcription, while the `Ribosome` class manages translation. 

The throughput of DNA sequences processed concurrently during protein synthesis is regulated by the `Resource` object. Transcription relies on the availability of `RNA_polymerase` as a critical resource, while translation demands ribosomes and `RNA_transfer` molecules with the appropriate anticodons. Additionally, `Nucleotides` are crucial resources present throughout the process. 

All resources have been integrated into the simulation framework through extensions of `SimPy`'s resource or container classes.

For a visual representation of the [Protein Synthesis Process](sources/ProteinSynthesisProcess.png) modeling, refer to the diagram provided in the [sources folder](sources/).

## Experiments
To test the model, experiments were conducted and are detailed in the following notebooks:
1. `12h_simulation.ipynb`: This notebook presents a 12-hour simulation of protein synthesis aimed at testing the model and examining the resulting outcomes. [Available here](experiments/12h_simulation.ipynb).
2. `comparative_analysis.ipynb`: This notebook contains a comparative analysis of models utilizing different resources. The objective of these experiments is to evaluate how the model's performance, in terms of the number of synthesized proteins and execution times, varies with the number of available resources. [Available here](experiments/comparative_analysis.ipynb).

## Dataset
To conduct the experiments, a dataset containing sequences from the human genome was utilized. The dataset was obtained from the RefSeqGene section of the Reference Sequence Database for the [Homo sapiens gene](https://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/) from the [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) Database.

For more detailed information and access to the dataset, please refer to the [repository](https://github.com/GiuliaGhisolfi/HumanGenomeDataset) where it is stored.

## Repository structure (main elements)
```
ProteinSynthesis/
│
├─── HumanGenomeDataset/                 # Repository contains a dataset loaded from the RefSeq Database
│
├─── data/                               # Data files used in simulations and experiments' parameters
│    ├─── codons.json
│    ├─── parameters_ribosome.json
│    ├─── parameters_rna_polymerases.json
│    └─── peptides.json
│
├─── experiments/
│    ├─── 12h_simulation.ipynb                # 12-hour simulation of protein synthesis
│    └─── comparative_analysis.ipynb          # Comparative analysis of models utilizing different resources
│
├─── src/                                # Source code files
│    ├─── process/
│    │    ├─── protein_synthesis.py      # EukarioticCell class, simulates the protein synthesis
│    │    ├─── transcription.py          # Nucleus class, simulates the transcription process
│    │    └─── translation.py            # Ribosome class, simulates the translation process
│    │
│    ├─── resources/
│    │    ├─── container.py              # EukaryoticCellContainer class
│    │    ├─── nucleotides.py            # Nucleotides class
│    │    ├─── resource.py               # EukaryoticCellResource class
│    │    └─── transfer_mrna.py          # TransferRNA class
│    │
│    ├─── utils/
│    │    ├─── plot_utils.py             # Function to visualize simulations' results
│    │    └─── utils.py                  # Utility function
│    │
│    ├─── variables/
│    │    ├─── nucleotide_allocations.py
│    │    └─── variables.py              # Class to store the variables of the simulation
│    │
│    └─── simulation.py                  # Class to simulate the protein synthesis process
│
└─── main.py                             # Main script to run experiments
```
