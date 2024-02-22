import simpy
import random

def process(env, sentence, t_inter, sentence_number, machines):
    with machines.request() as request:
        print(f'Sentence {sentence_number} requesting at {env.now}')
        yield request
        print(f'Sentence {sentence_number} got resource at {env.now}')

        words = sentence.split()
        for word in words:
            word = ''.join([letter for letter in word if letter not in 'aeiou'])

        yield env.timeout(t_inter*3)
        print(f'Sentence {sentence_number} splited at {env.now}')
        machines.release(request)

def setup(env, sentences, num_machines, t_inter):
    machines = simpy.Resource(env, num_machines)
    queue = []

    for sentence_number, sentence in enumerate(sentences):
        queue.append(env.process(process(env, sentence, t_inter, sentence_number, machines)))
        yield env.timeout(1)

    while queue:
        yield queue.pop(0)

# Setup and start the simulation
env = simpy.Environment()
sentences = ['Hello word', 'World cosa caso casa', 'Ciao come stai', 'Io sto bene', 'Tu come stai', 'Io sto'
    'bene', 'Tu come stai', 'Io sto bene', 'Tu come stai', 'Io sto bene', 'Tu come stai', 'Io sto bene', 'Tu come stai']
random.seed(42)
env.process(setup(env, sentences, num_machines=2, t_inter=1))
env.run(until=5)
