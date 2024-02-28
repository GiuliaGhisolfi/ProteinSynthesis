import simpy.resources.container as SimpyContainer
import json
import random

MIN_DEGRADATION_TIME = 60 # seconds
MAX_DEGRADATION_TIME = 180 # seconds

class EucaryotesCellContainer(SimpyContainer.Container):
    def __init__(self, env, capacity, init, random_seed):
        super().__init__(env, capacity, init)
        random.seed(random_seed)
        self._reset_history()
        
    def get(self, *args, **kwargs):
        get = super().get(*args, **kwargs) # get the amount from the container

        # save level and time in the history
        self._history['level'].append(self.level)
        self._history['time'].append(self._env.now)

        return get
    
    def put(self, amount, **kwargs):
        time = self._env.now
        degradation_time = round(random.uniform(
            MIN_DEGRADATION_TIME*amount, MAX_DEGRADATION_TIME*amount), 4)
        yield self._env.timeout(degradation_time) # FIXME

        super().put(amount, **kwargs) # put the amount into the container

        # save level and time in the history
        self._history['level'].append(self.level)
        self._history['time'].append(self._env.now)
    
    def level_history(self):
        return self._history
    
    def save_history(self, path_to_save):
        with open(path_to_save, 'w') as outfile:
            json.dump(self.level_history(), outfile)
        
    def _reset_history(self):
        self._history = {
            'level': [self.level],
            'time': [self._env.now]
            }