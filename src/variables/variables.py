
class EucaryotesCellVariables:
    def __init__(self):
        # input variables
        self.dna_sequence = None # template strand (3' to 5' direction)

        # intermediate variables
        self.sequence_count = None
        self.dna_sequences_to_transcript_list = None
        self.promoters_count = None

        # time variables
        self.request_start_process_time = None
        self.start_process_time = None
        self.start_transcription_time = []
        self.start_translation_time = []
        self.end_translation_time = []
        self.end_process_time = None

        # output variables
        self.mrna_sequences_list = []
        self.proteins_list = [] # polypeptides chains
        self.proteins_extended_name_list = []
        
        # simulation variables: list of simpy processes
        self.transcription_queue = []
        self.translation_queue = []
    
    def get_dna(self):
        return self.dna_sequence
    
    def get_mrna(self):
        return self.mrna_sequences_list
    
    def get_proteins(self):
        return self.proteins_list
    
    def get_extended_proteins_name(self):
        return self.proteins_extended_name_list