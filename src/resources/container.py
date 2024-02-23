import simpy.resources.container as SimpyContainer
import json

class EucaryotesCellContainer(SimpyContainer.Container):
    def __init__(self, env, capacity, init):
        super().__init__(env, capacity, init)
        self._reset_queue_history()
    
    def put(self, *args, **kwargs):
        put = super().put(*args, **kwargs)
        self._queue_history['queue'].append(self.level)
        self._queue_history['put_time'].append(self._env.now)
        return put
    
    def get(self, *args, **kwargs):
        get = super().get(*args, **kwargs)
        self._queue_history['get_time'].append(self._env.now)
        self._queue_history['wait_time'].append(
            self._queue_history['get_time'][-1] - self._queue_history['put_time'][-1])
        return get
    
    def queue_history(self):
        return self._queue_history
    
    def save_history(self, path_to_save):
        with open(path_to_save, 'w') as outfile:
            json.dump(self.queue_history(), outfile)
    
    def _reset_queue_history(self):
        self._queue_history = {
            'queue': [],
            'put_time': [],
            'get_time': [],
            'wait_time': []
            }