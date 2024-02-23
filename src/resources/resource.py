import simpy.resources.resource as SimpyResource
import json

class EucaryotesCellResource(SimpyResource.Resource):
    def __init__(self, env, capacity):
        super().__init__(env, capacity=capacity)
        self._reset_queue_history()
    
    def request(self, *args, **kwargs):
        request = super().request(*args, **kwargs)
        self._queue_history['queue'].append(len(self.queue))
        self._queue_history['request_time'].append(self._env.now)
        return request

    #TODO: implement method to track the wait time for each request
    """self._queue_history['available_time'].append(self._env.now)
    self._queue_history['wait_time'].append(
        self._queue_history['available_time'][-1] - self._queue_history['request_time'][-1]) """
    
    """
    self._queue_history['end_time'].append(self._env.now)
    self._queue_history['usage_time'].append(
        self._queue_history['end_time'][-1] - self._queue_history['available_time'][-1])
        """
        
    def queue_history(self):
        return self._queue_history
    
    def save_history(self, path_to_save):
        with open(path_to_save, 'w') as outfile:
            json.dump(self.queue_history(), outfile)
    
    def _reset_queue_history(self):
        self._queue_history = {
            'queue': [],
            'request_time': [], # time when the request is made
            'available_time': [], # time when the resource is available
            'wait_time': [], # time the request waited in the queue
            'end_time': [], # time when the resource is released
            'usage_time': [] # time the resource is used
            }
        
