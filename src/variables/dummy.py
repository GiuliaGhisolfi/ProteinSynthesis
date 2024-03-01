
class EucaryotesCellVariables:
    def __init__(self):
        # input variables
        self.dna_sequence = None # template strand (3' to 5' direction)
        self.promoters_box = None # promoters found in the DNA sequence
        self.dna_sequences_to_transcript_list = None

    def get_dna(self):
        return self.dna_sequence
    
    def get_mrna(self):
        return self.mrna_sequences_list
    
    def get_proteins(self):
        return self.proteins_list
    
    def get_extended_proteins_name(self):
        return self.proteins_extended_name_list
    
class mRNAVariables:
    def __init__(self):
        # output variables
        self.mrna_sequences = None
        self.proteins = None # polypeptides chains
        self.proteins_extended_name = None
        self.proteins_sintetized = 0

        # variables
        self.sequence_count = None
        self.promoters_count = None
        self.poly_adenine_tail_len = None
        self.mrna_degradation_rate = 1e-4
        
        # time variables
        self.request_start_process_time = None
        self.start_process_time = None
        self.start_transcription_time = None
        self.start_translation_time = []
        self.end_translation_time = []
        self.end_process_time = None
        
        # simulation variables: list of simpy processes
        self.transcription_queue = []
        self.complement_base_queue_dict = dict()
        self.polyadenylation_queue = []
        self.translation_queue = []