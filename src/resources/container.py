import simpy.resources.container as SimpyContainer
import json

class EucaryotesCellContainer(SimpyContainer.Container):
    def __init__(self, env, capacity, init):
        super().__init__(env, capacity, init)
        self._reset_history()
        
    def get(self, *args, **kwargs):
        get = super().get(*args, **kwargs) # get the amount from the container

        # save level and time in the history
        self._history['level'].append(self.level)
        self._history['time'].append(self._env.now)

        return get
    
    def put(self, *args, **kwargs):
        put = super().put(*args, **kwargs) # put the amount into the container

        # save level and time in the history
        self._history['level'].append(self.level)
        self._history['time'].append(self._env.now)

        return put
    
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