import simpy.resources.resource as SimpyResource
import json

class EukaryoticCellResource(SimpyResource.Resource):
    """
    This class extends the simpy Resource class to add the history of the queue
    and the usage of the resource. The history is saved in a dictionary with the
    queue, request time, available time, wait time, end time and usage time.

    Parameters:
    -----------
    env : simpy.Environment
        The simulation environment
    capacity : int
        The capacity of the resource
    save_history : bool
        If True, the history of the queue and the usage of the resource is saved

    Attributes:
    -----------
    _queue_history : dict
        The history of the queue and the usage of the resource

    Methods:
    --------
    request(*args, **kwargs)
        Request the resource
    available()
        Save the time when the resource is available
    release(*args, **kwargs)
        Release the resource
    queue_history()
        Return the queue history
    save_history(path_to_save)
        Save the queue history in a json file
    """
    def __init__(self, env, capacity, save_history=True):
        super().__init__(env, capacity=capacity)
        self.save_history_flag = save_history

        if self.save_history_flag:
            self._reset_queue_history()
    
    def request(self, *args, **kwargs):
        request = super().request(*args, **kwargs)

        if self.save_history_flag:
            self._queue_history['queue'].append(len(self.queue))
            self._queue_history['request_time'].append(self._env.now)

        return request
    
    def available(self):
        if self.save_history_flag:
            self._queue_history['available_time'].append(self._env.now)
            self._queue_history['wait_time'].append(
                self._queue_history['available_time'][-1] - self._queue_history['request_time'][-1])

    def release(self, *args, **kwargs):
        release = super().release(*args, **kwargs)

        if self.save_history_flag:
            self._queue_history['end_time'].append(self._env.now)
            self._queue_history['usage_time'].append(
                self._queue_history['end_time'][-1] - self._queue_history['available_time'][-1])
   
        return release

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
        
