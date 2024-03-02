
class EucaryotesCellVariables:
    """
    Class to store the variables of the eucaryotes cell simulation.

    Attributes
    ----------
    dna_sequence : str
        The DNA
    sequence_count : int
        Counter of the DNA sequence currently being processed
    dna_sequences_to_transcript_list : list
        List of DNA sequences to transcript to mRNA, each element found in the DNA sequencebetween two promoters
    promoters_count : int
        Number of the promoters found in the DNA sequence
    poly_adenine_tail_len : list
        Length of the poly-adenine tail for each mRNA sequence
    mrna_degradation_rate : list
        Rate of mRNA degradation for each mRNA sequence, initially set to 1e-4, 
        it increases with the number of ribosomes attached to the mRNA
    request_start_process_time : float
        Time when the simulation was requested to start
    start_process_time : float
        Time when the simulation started
    self.start_transcription_time: list
        Time when the transcription of each mRNA sequence started
    start_translation_time: list
        Time when the translation of each mRNA sequence started
    end_translation_time: list
        Time when the translation of each mRNA sequence ended
    end_process_time: float
        Time when the simulation ended
    promoters_box: list
        List of promoters found in the DNA sequence
    mrna_sequences_list: list
        List of mature mRNA sequences
    proteins_list: list
        List of polypeptides chains
    proteins_extended_name_list: list
        List of polypeptides chains extended name
    proteins_sintetized: list
        List of number of polypeptides chains sintetized for each mRNA sequence
    transcription_queue: list
        List of simpy processes for the transcription, 
        this list is used to store and relase the simpy processes for each mRNA sequence
    complement_base_queue_dict: dict
        Dictionary of list of simpy processes for each complement base,
        this dictionary is used to store and relase the simpy processes
    polyadenylation_queue: list
        List of simpy processes for the polyadenylation, 
        this list is used to store and relase the simpy processes for each mRNA sequence
    translation_queue: list
        List of simpy processes for the translation, 
        this list is used to store and relase the simpy processes for each mRNA sequence

    """
    def __init__(self):
        # input variables
        self.dna_sequence = None # template strand (3' to 5' direction)

        # intermediate variables
        self.sequence_count = None
        self.dna_sequences_to_transcript_list = None
        self.promoters_count = None
        self.poly_adenine_tail_len = None
        self.mrna_degradation_rate = None # 1e-4

        # time variables
        self.request_start_process_time = None
        self.start_process_time = None

        self.start_transcription_time = []
        self.start_translation_time = []
        self.end_translation_time = []

        self.end_process_time = None

        # others variables
        self.promoters_box = None # promoters found in the DNA sequence

        # output variables
        self.mrna_sequences_list = None
        self.proteins_list = None # polypeptides chains
        self.proteins_extended_name_list = None
        self.proteins_sintetized = None
        
        # simulation variables: list of simpy processes
        self.transcription_queue = []
        self.complement_base_queue_dict = dict()
        self.polyadenylation_queue = []
        self.translation_queue = []

    def init_transcription_translation_var(self):
        # variables
        self.mrna_sequences_list = [None] * len(self.dna_sequences_to_transcript_list)
        self.poly_adenine_tail_len = [0] * len(self.dna_sequences_to_transcript_list)
        self.proteins_list = [None] * len(self.dna_sequences_to_transcript_list)
        self.proteins_extended_name_list = [None] * len(self.dna_sequences_to_transcript_list)

        # time variables #TODO (?)
        """self.start_transcription_time = dict()
        self.start_translation_time = []
        self.end_translation_time = []"""

        self.proteins_sintetized = [0] * len(self.dna_sequences_to_transcript_list)
        self.mrna_degradation_rate = [1e-4] * len(self.dna_sequences_to_transcript_list)

    def get_dna(self):
        return self.dna_sequence
    
    def get_mrna(self):
        return self.mrna_sequences_list
    
    def get_proteins(self):
        return self.proteins_list
    
    def get_extended_proteins_name(self):
        return self.proteins_extended_name_list