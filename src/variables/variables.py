
class EucaryotesCellVariables:
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