import os
import random
from tqdm import tqdm

hor_len_list = [10]
mono_diff_rate_list = [0.1, 0.15]
muta_rate_list = [0.005, 0.01, 0.015]
bases = ['A', 'C', 'G', 'T']
repeats = {'chicken':'GCGTTTTCTCTTCGCAAATCCCCCATTTAAGCGAAAATCAC',
            'human':'AATCTGCAAGTGGACATTTGGAGCGCTTTGAGGCCTATGGTGGAAAAGGAAATATCTTCACATAAAAACTAGACAGAAGCATTCTCAGAAACTTCTTTGTGATGTGTGCATTCAACTCACAGAGTTGAACCTTTCTTTTGATAGAGCAGTTTTGAAACACTCTTTTTGTAG'}

random_num = 3
hor_num = 500
mutation_type_list = ['insertion'] * 20 + ['deletion'] * 20 + ['substitution'] * 60

def mutate_sequence(sequence, mutation_rate):
    mutated_sequence = ''
    for base in sequence:
        if random.random() < mutation_rate:
            mutation_type = random.choice(mutation_type_list)
            if mutation_type == 'insertion':
                inserted_base = random.choice(bases)
                mutated_sequence += inserted_base + base
            elif mutation_type == 'deletion':
                pass
            else:
                new_base = random.choice(bases)
                while new_base == base:
                    new_base = random.choice(bases)
                mutated_sequence += new_base
        else:
            mutated_sequence += base
    return mutated_sequence

def main(horLen = 5, horNum = 100, continuous = True):
    for i in tqdm(range(random_num)):  # Number of random simulations
        for key in repeats:
            rep = repeats[key]
            for diff in mono_diff_rate_list:
                for hor_len in hor_len_list:
                    hor = ''
                    for j in range(hor_len):
                        muta_mono = mutate_sequence(rep, diff)
                        hor += muta_mono
                        # pass

                    for muta in muta_rate_list:
                        final_hor_array = ''
                        hor_array  = ''
                        for t in range(hor_num):
                            # pass
                            hor_array += hor

                        final_hor_array = mutate_sequence(hor_array, muta)
                        _file = open(f'{key}_{hor_len}_{diff}_{muta}_{i}.fa', 'w')
                        _file.write(f'>hor\n')
                        _file.write(f'{final_hor_array}\n')
                        _file.close()
    
if __name__ == '__main__':
    main()