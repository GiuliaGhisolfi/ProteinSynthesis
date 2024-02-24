import simpy.resources.container as SimpyContainer
import json

class EucaryotesCellContainer(SimpyContainer.Container):
    def __init__(self, env, capacity, init):
        super().__init__(env, capacity, init)
        self._reset_queue_history()
        self._reset_amount_history(initial_amount=init)
        
    def get(self, *args, **kwargs):
        get = super().get(*args, **kwargs) # get the amount from the container
        
        # save the amount history
        self._amount_history['amount'].append(self.level)
        self._amount_history['time'].append(self._env.now)

        return get
    
    def put(self, *args, **kwargs):
        put = super().put(*args, **kwargs) # put the amount into the container

        # save the amount history
        self._amount_history['amount'].append(self.level)
        self._amount_history['time'].append(self._env.now)

        return put
    
    def queue_history(self):
        return self._queue_history
    
    def amount_history(self):
        return self._amount_history
    
    def save_history(self, path_to_save):
        with open(path_to_save, 'w') as outfile:
            json.dump(self.queue_history(), outfile)
        with open(path_to_save.replace('.json', '_amount.json'), 'w') as outfile:
            json.dump(self.amount_history(), outfile)
    
    def _reset_queue_history(self):
        self._queue_history = {
            'request_time': [],
            'get_time': [],
            'wait_time': []
            } #TODO: implement method to track the wait time for each request
        
    def _reset_amount_history(self, initial_amount):
        self._amount_history = {
            'amount': [initial_amount],
            'time': [self._env.now]
            }