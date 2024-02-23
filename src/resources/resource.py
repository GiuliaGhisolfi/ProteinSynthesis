import simpy.resources.resource as SimpyResource

class EucaryotesCellResource(SimpyResource.Resource):
    def __init__(self, env, number_ribosomes):
        super().__init__(env, capacity=number_ribosomes)
        self._reset_queue_history()

    def queue_history(self):
        return self._queue_history
    
    def request(self, *args, **kwargs):
        request = super().request(*args, **kwargs)
        self._queue_history['queue'].append(len(self.queue))
        self._queue_history['request_time'].append(self._env.now)
        return request
    
    def release(self, *args, **kwargs):
        release = super().release(*args, **kwargs)
        self._queue_history['release_time'].append(self._env.now)
        self._queue_history['wait_time'].append(
            self._queue_history['release_time'][-1] - self._queue_history['request_time'][-1])
        return release
    
    def _reset_queue_history(self):
        self._queue_history = {
            'queue': [],
            'request_time': [],
            'release_time': [], 
            'wait_time': []
            }
        
