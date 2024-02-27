import simpy.resources.container as SimpyContainer
import json
import os
TIME_UNIT = 1e-4 # 0.0001 seconds

class EucaryotesCellContainer(SimpyContainer.Container):
    def __init__(self, env, capacity, init):
        super().__init__(env, capacity, init)
        self._reset_history()
        env.process(self.monitor_container(env))
        
    def get(self, *args, **kwargs):
        get = super().get(*args, **kwargs) # get the amount from the container
        return get
    
    def put(self, *args, **kwargs):
        put = super().put(*args, **kwargs) # put the amount into the container
        return put
    
    def level_history(self):
        return self._history
    
    def save_history(self, path_to_save):
        with open(path_to_save, 'w') as outfile:
            json.dump(self.level_history(), outfile)
        
    def _reset_history(self):
        self._history = {
            'level': [],
            'time': []
            }

    def monitor_container(self, env): #TODO: change, computazionalmente troppo lento, fare preproc dopo
        while True:
            self._history['level'].append(self.level)
            self._history['time'].append(env.now)
            yield env.timeout(TIME_UNIT)