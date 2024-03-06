# Protein Synthesis in Eukaryotic Cells 🧬🧫
This repository contains the code for a project aimed at modeling protein synthesis in eukaryotic cells as a discrete-time process.

## Modeling
The process is simulated using `SimPy` framework, with customized classes tailored to the specific requirements. 

The process is modeled by the `EukaryoticCell` class, which coordinates transcription and translation through its components: the `Nucleus` class handles transcription, while the `Ribosome` class manages translation. 

The throughput of DNA sequences processed concurrently during protein synthesis is regulated by the `Resource` object. Transcription relies on the availability of `RNA polymerase` as a critical resource, while translation demands ribosomes and `RNA transfer` molecules with the appropriate anticodons. Additionally, nucleotides are crucial resources present throughout the process. 

All resources have been integrated into the simulation framework through extensions of `SimPy`'s resource or container classes.

For a visual representation of the [Protein Synthesis Process](sources/ProteinSynthesisProcess.png) modeling, refer to the diagram provided in the [sources folder](sources/).

## Experiments
Al fine di testare il modello sono stati svolti degli esperimenti riportati in 

Notebook ['12h_simulation'](12h_simulation.ipynb) contains

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
├─── 12h_simulation.ipynb                # 12-hour simulation of protein synthesis
├─── comparative_analysis.ipynb          # Comparative analysis of models utilizing different resources
└─── main.py                             # Main script to run experiments
```
