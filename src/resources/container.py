import simpy.resources.container as SimpyContainer
import json
import random

MIN_DEGRADATION_TIME = 60 # seconds
MAX_DEGRADATION_TIME = 180 # seconds

class EukaryoticCellContainer(SimpyContainer.Container):
    """
    This class extends the simpy Container class to add the degradation of the
    amount in the container. The degradation time is a random value between
    MIN_DEGRADATION_TIME and MAX_DEGRADATION_TIME. The degradation time is
    simulated using a timeout event in the simpy environment.
    The level and time when the amount is put or get from the container are
    saved in a history dictionary.
    
    Parameters:
    -----------
    env : simpy.Environment
        The simulation environment
    capacity : int
        The capacity of the container
    init : int
        The initial amount in the container
    random_seed : int
        The random seed for the degradation time

    Attributes:
    -----------
    _history : dict
        The history of the level and time in the container

    Methods:
    --------
    get(*args, **kwargs)
        Get the amount from the container
    put(amount)
        Put the amount into the container
    level_history()
        Return the level history
    save_history(path_to_save)
        Save the level history in a json file
    """
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
    
    def put(self, amount):
        degradation_time = round(random.uniform(
            MIN_DEGRADATION_TIME, MAX_DEGRADATION_TIME), 4)
        yield self._env.timeout(degradation_time)

        super().put(amount) # put the amount into the container

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