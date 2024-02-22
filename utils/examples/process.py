import simpy
import random

def process(env, sentence, t_inter, sentence_number):
    words = sentence.split()
    for word in words:
        word = ''.join([letter for letter in word if letter not in 'aeiou'])

    yield env.timeout(t_inter*3)
    print(f'Sentence {sentence_number} splited at {env.now}')
    

def setup(env, sentences, num_machines, t_inter):
    machines = simpy.Resource(env, num_machines)
    sentence_number = 0

    while True:
        sentence = sentences[sentence_number]

        while sentence_number < len(sentences):
            with machines.request() as request:
                print(f'Sentece {sentence_number} requested at {env.now}')
                yield request
                print(f'Sentece {sentence_number} got resource at {env.now}')

                env.process(process(env, sentence, t_inter, sentence_number))
                sentence_number += 1
                yield env.timeout(t_inter)

# Setup and start the simulation
env = simpy.Environment()
sentences = ['Hello word', 'World cosa caso casa', 'Ciao come stai', 'Io sto bene', 'Tu come stai', 'Io sto'
    'bene', 'Tu come stai', 'Io sto bene', 'Tu come stai', 'Io sto bene', 'Tu come stai', 'Io sto bene', 'Tu come stai']
random.seed(42)
env.process(setup(env, sentences, num_machines=1, t_inter=1))
env.run(until=5)