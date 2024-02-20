import simpy
import itertools
import random

def process(env, sentences):
    available = {sentence: True for sentence in sentences}
    print(f'Process started at {env.now} \n')

    while True:
        yield env.timeout(random.random())
        sentence = random.choice(sentences)

        if available[sentence]:
            print(f'\nTime {env.now}: {sentence}')
            env.process(splite_sentence(env, sentence))
            available[sentence] = False

def splite_sentence(env, sentence):
    words = sentence.split()
    print(f'Time {env.now:2}: {words}')
    yield env.timeout(2)

    for word in words:
        env.process(delete_vocal(env, word))
    
def delete_vocal(env, word):
    word = ''.join([letter for letter in word if letter not in 'aeiou'])
    print(f'Time {env.now}: {word}')
    yield env.timeout(1)

def setup(env, sentences, num_machines, t_inter):
    for _ in range(num_machines):
        env.process(process(env, sentences))

    while True:
        yield env.timeout(t_inter)
        env.process(process(env, sentences))

# Setup and start the simulation
env = simpy.Environment()
sentences = ['Hello word', 'World cosa caso casa']
random.seed(42)
env.process(setup(env, sentences, num_machines=2, t_inter=1))
env.run(until=10)