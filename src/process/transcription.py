import random
from src.variables.nucleotides_allocations import NucleotidesSymbolsAllocations
from src.resources.resource import EukaryoticCellResource

BASE_COMPLEMENT_DNA2RNA = {
    'A': 'U', 
    'T': 'A', 
    'C': 'G', 
    'G': 'C',
}
BASE_COMPLEMENT_RNA2DNA = {
    'U': 'A', 
    'A': 'T', 
    'G': 'C', 
    'C': 'G'
}
PROMOTERS = {
    'TATAbox': ['TATAAAA', 'TATAAAT', 'TATATAA', 'TATATAT'], # TATA box
    'CAATbox': ['CCAAT'], # CAAT box
    'GCbox': ['GGGCGG'], # GC box
}
LENGTH_PROMOTER = {
    'TATAbox': 7,
    'CAATbox': 5,
    'GCbox': 6,
}
START_CODEN = 'AUG'
TERMINATORS = ['UAA', 'UAG', 'UGA']
LENGTH_EXTRON_SEQUENCE = 3 # length of extron sequence
LENGTH_METHYL_CAP = 8 # length of 5'-methyl cap
MIN_LENGTH_PROMOTER = 200 # minimum length between promoters
REPLICATION_TIME = 2e-2 # seconds to replicate a nucleotide
CLEAVAGE_TIME = 1e-2 # seconds to cleave the mRNA
TRANSCRIPTION_TIMEOUT = 1
RNA_POLYMERASE_ERROR_RATE = 10e-4 # 1 error per 10^4 nucleotides

class Nucleus:
    """
    Nucleus class, this class models the nucleus of the eukaryotic cell 
    and simulates the transcription process.

    Parameters
    ----------
    environment : simpy.Environment
        The simulation environment.
    extron_sequences_list : list
        The list of extron sequences.
    editing_sites_dict : dict
        The dictionary of editing sites.
    number_rna_polymerases : int
        The number of RNA polymerases in the cell.
    nucleotides : Nucleotides
        The nucleotides in the cell.
    random_seed : int
        The random seed for the simulation environment.

    Attributes
    ----------
    env : simpy.Environment
        The simulation environment.
    extron_sequences_list : list
        The list of extron sequences.
    editing_sites_dict : dict
        The dictionary of editing sites.
    rna_polymerase : EukaryoticCellResource
        The RNA polymerase resource.
    nucleotides : Nucleotides
        The nucleotides in the cell.

    Methods
    -------
    find_promoter(dna_sequence, variables)
        Find the promoter of a DNA sequence.
    find_promoters_positions(dna_sequence, promoters_list)
        Find the positions of the promoters in a DNA sequence.
    transcript(dna_sequence, variables, seq_count)
        Start the transcription of a DNA sequence.
    transcript_process(dna_sequence, variables, seq_count)
        Start the transcription process of a DNA sequence.
    trascript_gene(dna_sequence, variables, sequence_count)
        Start the transcription of a gene.
    find_complement_base(base)
        Find the complement base of a base.
    splicing(rna_sequence)
        Remove the introns from a RNA sequence.
    editing(rna_sequence)
        Edit a RNA sequence, replacing the editing sites with the edited sites.
    capping(rna_sequence)
        Add a 5'-methyl cap to a RNA sequence.
    cleavage()
        Cleavage a RNA sequence, start post-transcriptional modifications.
    polyadenylation(rna_sequence, variables, seq_count)
        Add a PolyA tail to a RNA sequence.
    find_adenosine_for_polyadenylation(amount)
        Find the adenosine for the polyadenylation.
    request_nucleotide(base, amount=1)
        Request a nucleotide.
    replicate_base(base)
        Replicate a nucleotide.
    release_nucleotide(base, amount=1)
        Release a nucleotide.
    """
    def __init__(self, environment, extron_sequences_list, editing_sites_dict, 
            number_rna_polymerases, nucleotides, random_seed):
        self.env = environment

        self.extron_sequences_list = extron_sequences_list
        self.editing_sites_dict = editing_sites_dict
        self.editing_sites_dict = dict(sorted(self.editing_sites_dict.items(), 
            key=lambda x: len(x[0]), reverse=False)) # sort by length of key
        
        self.rna_polymerase = EukaryoticCellResource(self.env, capacity=number_rna_polymerases)
        self.nucleotides = nucleotides

        random.seed(random_seed)
    
    def find_promoter(self, dna_sequence, variables):
        """
        Find the promoter of a DNA sequence, if present.
        The promoter is the region of a DNA sequence where the RNA polymerase binds to start the transcription.
        The DNA sequence is split in the promoter regions.
        """
        for promoter in PROMOTERS:
            promoter_positions_list = self.find_promoters_positions(dna_sequence, PROMOTERS[promoter])
            if len(promoter_positions_list) != 0:
                variables.promoters_box = promoter
                promoter_length = LENGTH_PROMOTER[promoter]
                break
        
        if len(promoter_positions_list) == 0:
            return None
        else:
            # split the DNA sequence in the promoter regions
            dna_sequences_list = []

            for i in range(len(promoter_positions_list)-1):
                dna_sequences_list.append(dna_sequence[promoter_positions_list[i]+promoter_length
                    :promoter_positions_list[i+1]])
                
            dna_sequences_list.append(dna_sequence[promoter_positions_list[-1]+promoter_length:])
            
            return dna_sequences_list
    
    def find_promoters_positions(self, dna_sequence, promoters_list):
        """
        Find the positions of the promoters in a DNA sequence.
        """
        # find promoter sequences in the DNA sequence
        promoter_positions_list = sorted([i for promoter in promoters_list for i, _ in 
            enumerate(dna_sequence) if dna_sequence[i:].startswith(promoter)])
        
        # if dist between promoters is less than 100, consider it as a single promoter
        i = 0
        while i < len(promoter_positions_list)-1: 
            if promoter_positions_list[i+1] - promoter_positions_list[i] < MIN_LENGTH_PROMOTER:
                del promoter_positions_list[i+1]  
            else: i += 1
        
        return promoter_positions_list

    def transcript(self, dna_sequence, variables, seq_count): # enzime: RNA polymerase
        """
        Start the transcription process of a DNA sequence.
        """
        with self.rna_polymerase.request() as request:
            yield request  # wait for RNA polymerase to be available

            # start transcript processes for DNA sequence
            messenger_rna_sequence = yield self.env.process(
                self.transcript_process(dna_sequence, variables, seq_count))

        return messenger_rna_sequence
    
    def transcript_process(self, dna_sequence, variables, seq_count):
        """
        Transcription process of a DNA sequence.
        """
        self.rna_polymerase.available() # register the time when the resource is available
        
        # make sequence univoque to transcript
        dna_sequence = ''.join([random.choice(NucleotidesSymbolsAllocations[n]) for n in dna_sequence])

        # transcript from gene to pre-mRNA
        messenger_rna_sequence = yield self.env.process(
            self.trascript_gene(dna_sequence, variables, seq_count))

        # pre-transcriptional modifications
        messenger_rna_sequence = self.capping(messenger_rna_sequence)

        # elongation phase
        messenger_rna_sequence = self.splicing(messenger_rna_sequence)
        messenger_rna_sequence = self.editing(messenger_rna_sequence)
        
        # post-transcriptional modifications
        yield self.env.process(self.cleavage())
        messenger_rna_sequence = yield self.env.process(
            self.polyadenylation(messenger_rna_sequence, variables, seq_count))

        yield self.env.timeout(TRANSCRIPTION_TIMEOUT)

        return messenger_rna_sequence
    
    def trascript_gene(self, dna_sequence, variables, sequence_count):
        """
        Transcription from a gene to a pre-mRNA, the first step of the transcription process.
        The DNA sequence is transcribed to a messenger RNA sequence.
        """
        messenger_rna_sequence = ''

        for base in dna_sequence:
            complement_base_process = yield self.env.process(self.find_complement_base(base))
            variables.complement_base_queue_dict[sequence_count].append(complement_base_process)

            yield complement_base_process
            complement_base = complement_base_process.value
            messenger_rna_sequence += complement_base

            # wait for all the transcription gene process to be completed
            while variables.complement_base_queue_dict[sequence_count]:
                variables.complement_base_queue_dict[sequence_count].pop(0)

        return messenger_rna_sequence
    
    def find_complement_base(self, base):
        """
        Find the complement base of a base.
        """
        if random.random() > RNA_POLYMERASE_ERROR_RATE:
            complement_base = BASE_COMPLEMENT_DNA2RNA[base]
        else: 
            complement_base = random.choice([b for b in list(BASE_COMPLEMENT_DNA2RNA.values()) 
                if b != BASE_COMPLEMENT_DNA2RNA[base]])

        return self.request_nucleotide(complement_base)
    
    def splicing(self, rna_sequence):
        """
        Remove the introns (non coding regions) from a RNA sequence.
        """
        # remove introns: non-coding regions
        i = LENGTH_METHYL_CAP # index
        while i+3 < len(rna_sequence):
            if rna_sequence[i:i+LENGTH_EXTRON_SEQUENCE] in self.extron_sequences_list:
                i += LENGTH_EXTRON_SEQUENCE
            else: 
                intron = rna_sequence[i]
                self.release_nucleotide(intron) # degrade intron
                rna_sequence = rna_sequence[:i] + rna_sequence[i+1:]

        return rna_sequence

    def editing(self, rna_sequence):
        """
        Edit a RNA sequence, replacing the editing sites with the edited sites.
        """
        for editing_site in self.editing_sites_dict.keys():
            # find the editing site in the rna sequence
            editing_site_count = rna_sequence.count(editing_site)
            for base in editing_site: # request the nucleotides to edit the rna sequence
                self.request_nucleotide(base, editing_site_count)
            
            # edit the rna sequence
            rna_sequence = rna_sequence.replace(editing_site, self.editing_sites_dict[editing_site])
            
            # degrade the editing site
            for base in editing_site:
                self.release_nucleotide(base, editing_site_count)

        return rna_sequence

    def capping(self, rna_sequence):
        """
        Add a 5'-methyl cap to a RNA sequence.
        """
        return 'CH3GPPP-{}'.format(rna_sequence) # Add 5'-methyl cap
    
    def cleavage(self):
        """
        Cleavage a RNA sequence, start post-transcriptional modifications.
        """
        yield self.env.timeout(CLEAVAGE_TIME)

    def polyadenylation(self, rna_sequence, variables, seq_count):
        """
        Add a PolyA tail to a RNA sequence.
        """
        variables.poly_adenine_tail_len[seq_count] = random.randint(230, 270)

        polyadenylation_process = self.env.process(
            self.find_adenosine_for_polyadenylation(variables.poly_adenine_tail_len[seq_count]))
        variables.polyadenylation_queue.append(polyadenylation_process)
        yield polyadenylation_process

        # wait for all the transcription gene process to be completed
        while variables.polyadenylation_queue:
                variables.polyadenylation_queue.pop(0)
        
        return '{}-AAAA'.format(rna_sequence) # Add PolyA tail (250 nucleotides circa)
    
    def find_adenosine_for_polyadenylation(self, amount):
        """
        Request the adenosine for the polyadenylation.
        """
        return self.request_nucleotide('A', amount)
    
    def request_nucleotide(self, base, amount=1):
        """
        Request a nucleotide.
        """
        with self.nucleotides.request(base, amount) as request:
            yield request
            
            return self.env.process(self.replicate_base(base))
    
    def replicate_base(self, base):
        """
        Replicate a nucleotide.
        """
        yield self.env.timeout(REPLICATION_TIME) # time to replicate a nucleotide
        return base
    
    def release_nucleotide(self, base, amount=1):
        """
        Release a nucleotide, degrade the nucleotide if not used.
        """
        self.nucleotides.release(base, amount)