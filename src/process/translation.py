from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio import BiopythonWarning
import warnings
import random
from src.resources.resource import EukaryoticCellResource
from src.resources.transfer_mrna import TransferRNA

START_CODON = 'AUG' # start codon
AMINO_GROUP = 'NH2-' # amino group
CARBOXYL_GROUP = '-COOH' # carboxyl group
LENGTH_CODON = 3 # number of nucleotides that code for an amino acid
LENGTH_METHYL_CAP = 8 # length of 5'-methyl cap
LENGTH_POLY_A_TAIL = 5 # length of poly-A tail
ATTIVATION_TIME = 1
TRANSFER_RNA_ATTACH_TIME = 1e-3 # seconds to attach each transfer RNA
ELONGATION_TIME = 5e-2 # seconds to add each amino acid
TRANSLATION_TIMEOUT = 10
MRNA_DECODED_ERROR_RATE = 1e-4 # 1 mistake every 10.000 amino acids

class Ribosome:
    """
    Ribosome class, this class models the ribosome and simulates the translation process.

    Parameters
    ----------
    environment : simpy.Environment
        The simulation environment.
    number_ribosomes : int
        The number of ribosomes in the cell.
    number_rna_transfers_per_codon : int
        Number of RNA transfer per codon in the simulation environment is a random integer in
        [0.9*number_rna_transfers_per_codon, 1.1*number_rna_transfers_per_codon].
    codons_list : list
        The list of codons.
    nucleotides : Nucleotides
        The nucleotides in the cell.
    amminoacids : list
        The list of amminoacids.
    random_seed : int
        The random seed for the simulation environment.

    Attributes
    ----------
    env : simpy.Environment
        The simulation environment.
    ribosomes : EukaryoticCellResource 
        The ribonucleoprotein complex in the cytoplasm.
    rna_transfer : TransferRNA
        The transfer RNA in the cell.
    nucleotides : Nucleotides
        The nucleotides in the cell.

    Methods
    -------
    translate(mrna_sequence, variables, seq_count)
        Start the translation process.
    translation_process(mrna_sequence, variables, seq_count)
        Start the translation process.
    degradation_cap_tail(mrna_sequence)
        Degradation of the 5' cap and poly-A tail, enzime: exonuclease.
    activation()
        Activation of the translation process.
    initialization(mrna_sequence)
        Initialization of the translation process.
    elongation(mrna_sequence)
        Elongation of the translation process.
    request_trna(codon)
        Request transfer RNA with the correct anticodon.
    compute_degradation_probability(mrna_sequence, mrna_degradation_rate)
        Compute the probability of the mRNA degradation.
    mrna_degradation(mrna_sequence, poly_adenine_tail_len)
        Degradation of the mRNA sequence, enzima: ribonuclease.
    release_nucleotide(base, amount=1)
        Release a nucleotide in the cell.
    """ 
    def __init__(self, environment, number_ribosomes, number_rna_transfers_per_codon,
            codons_list, nucleotides, amminoacids, random_seed):
        # ribonucleoprotein complex in the cytoplasm
        self.env = environment
        self.ribosomes = EukaryoticCellResource(self.env, capacity=number_ribosomes)
        self.rna_transfer = TransferRNA(self.env, amount=number_rna_transfers_per_codon,
            codons_list=codons_list, random_seed=random_seed)
        self.nucleotides = nucleotides
        self.amminoacids = amminoacids
        
        random.seed(random_seed)

    def translate(self, mrna_sequence, variables, seq_count): # protein synthesis
        """
        Start the translation process.
        """
        with self.ribosomes.request() as request:
            yield request # wait for a ribosome to be available
            
            polypeptides_chain, polypeptides_chain_ext, mrna_degradated = yield self.env.process(
                self.translation_process(mrna_sequence, variables, seq_count))
            
        return polypeptides_chain, polypeptides_chain_ext, mrna_degradated
    
    def translation_process(self, mrna_sequence, variables, seq_count):
        """
        Translation process.
        """
        self.ribosomes.available() # register the time when the resource is available
        variables.proteins_sintetized[seq_count] += 1

        # translation process
        mrna_sequence = self.degradation_cap_tail(mrna_sequence)
        initial_mrna_sequence = mrna_sequence

        yield self.env.process(self.activation())
        mrna_sequence = self.initialization(mrna_sequence)
        polypeptides_chain, polypeptides_chain_ext = yield self.env.process(
            self.elongation(mrna_sequence))
        
        yield self.env.timeout(TRANSLATION_TIMEOUT)

        # mRNA degradation
        if self.compute_degradation_probability(initial_mrna_sequence, 
            variables.mrna_degradation_rate[seq_count]) >= random.random():
            self.mrna_degradation(initial_mrna_sequence, variables.poly_adenine_tail_len[seq_count])
            mrna_degradated = True
        else:
            variables.mrna_degradation_rate[seq_count] += 1e-4
            variables.mrna_degradation_rate[seq_count] = min(
                variables.mrna_degradation_rate[seq_count], 1)
            mrna_degradated = False

        return polypeptides_chain, polypeptides_chain_ext, mrna_degradated
    
    def degradation_cap_tail(self, mrna_sequence):
        """
        Degradation of the 5' cap and poly-A tail, enzime: exonuclease.
        """
        return mrna_sequence[LENGTH_METHYL_CAP:-LENGTH_POLY_A_TAIL]
    
    def activation(self):
        """
        Activation of the translation process.
        """
        yield self.env.timeout(ATTIVATION_TIME)

    def initialization(self, mrna_sequence):
        """
        Initialization of the translation process, find the start codon for the translation.
        """
        start_codon = START_CODON # start codon
        start_codon_position = str(mrna_sequence).find(start_codon)

        return mrna_sequence[start_codon_position+LENGTH_CODON:]
    
    def elongation(self, mrna_sequence):
        """
        Elongation of the translation process.
        Request transfer RNA with the correct anticodon.
        """
        for i in range(int(len(mrna_sequence) // LENGTH_CODON)):
            j = i * LENGTH_CODON
            codon = mrna_sequence[j:j+LENGTH_CODON]
            with self.rna_transfer.trna_resources_dict[codon].request() as request:
                yield request
                yield self.env.process(self.request_trna(codon))

        # translation of the mRNA sequence
        mrna_sequence = Seq(mrna_sequence)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            polypeptides_chain_seq = mrna_sequence.translate(stop_symbol='', to_stop=True)
        
        # error in the translation process
        polypeptides_chain_list = list(polypeptides_chain_seq)
        polypeptides_chain = ''
        for amminoacid in polypeptides_chain_list:
            polypeptides_chain = polypeptides_chain + (random.choice([a for a in list(self.amminoacids) 
                 if a != amminoacid]) if random.random() <= MRNA_DECODED_ERROR_RATE else amminoacid)

        if len(polypeptides_chain) > 0:
            yield self.env.timeout(ELONGATION_TIME * len(polypeptides_chain)) # 0.05 seconds to add each amino acid
            polypeptides_chain_ext = seq3(polypeptides_chain)

            polypeptides_chain = AMINO_GROUP + polypeptides_chain + CARBOXYL_GROUP
            polypeptides_chain_ext = AMINO_GROUP + polypeptides_chain_ext + CARBOXYL_GROUP

            return str(polypeptides_chain), str(polypeptides_chain_ext)
        else:
            return None, None
    
    def request_trna(self, codon):
        """
        Request transfer RNA with the correct anticodon.
        """
        self.rna_transfer.trna_resources_dict[codon].available()
        yield self.env.timeout(TRANSFER_RNA_ATTACH_TIME)

    def compute_degradation_probability(self, mrna_sequence, mrna_degradation_rate):
        """
        Compute the probability of the mRNA degradation according to the length of the mRNA sequence.
        """
        return 1 - (1 - mrna_degradation_rate) ** len(mrna_sequence)

    def mrna_degradation(self, mrna_sequence, poly_adenine_tail_len):
        """
        Degradation of the mRNA sequence, enzime: ribonuclease.
        """
        # enzima: ribonuclease
        [self.release_nucleotide(nucleotide) for nucleotide in mrna_sequence]
        self.release_nucleotide('A', poly_adenine_tail_len)
    
    def release_nucleotide(self, base, amount=1):
        """
        Release a nucleotide in the cell.
        """
        self.nucleotides.release(base, amount)