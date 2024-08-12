import random

def generate_random_dna_sequence(length):
    return ''.join(random.choices('ACGT', k=length))

def generate_tandem_repeat_unit(length):
    return ''.join(random.choices('ACGT', k=length))

def mutate_sequence(seq, mutation_rate=0.01):
    seq_list = list(seq)
    for i in range(len(seq_list)):
        if random.random() < mutation_rate:
            seq_list[i] = random.choice('ACGT'.replace(seq_list[i], ''))
    return ''.join(seq_list)



def generate_tandem_repeat_region(min_length, max_length, min_unit_length, max_unit_length, mutation_rate=0.01):
    region_length = random.randint(min_length, max_length)
    unit_length = random.randint(min_unit_length, max_unit_length)
    num_units = region_length // unit_length
    repeat_unit = generate_tandem_repeat_unit(unit_length)
    region = ''.join(mutate_sequence(repeat_unit, mutation_rate=mutation_rate) for _ in range(num_units))
    return region[:region_length], unit_length

def generate_dna_with_repeats(total_length, num_repeats, repeat_min_length, repeat_max_length, unit_min_length, unit_max_length, mutation_rate=0.01):
    remaining_length = total_length
    dna_sequence = ""
    repeat_info = []
    for _ in range(num_repeats):
        repeat_region, unit_length = generate_tandem_repeat_region(repeat_min_length, repeat_max_length, unit_min_length, unit_max_length, mutation_rate = mutation_rate)
        non_repeat_length = (remaining_length - len(repeat_region)) // (num_repeats - _)
        if non_repeat_length > 0:
            dna_sequence += generate_random_dna_sequence(non_repeat_length)
            remaining_length -= non_repeat_length
        else:
            dna_sequence += generate_random_dna_sequence(10000)  # final region
        start_position = len(dna_sequence)
        dna_sequence += repeat_region
        end_position = len(dna_sequence) - 1
        remaining_length -= len(repeat_region)
        repeat_info.append((start_position, end_position, unit_length))
    if remaining_length > 0:
        dna_sequence += generate_random_dna_sequence(remaining_length)
    return dna_sequence, repeat_info

def write_fasta(filename, sequence):
    with open(filename, 'w') as f:
        f.write(">simulated_sequence\n")
        for i in range(0, len(sequence), 80):  
            f.write(sequence[i:i+80] + "\n")



def write_repeat_info(filename, repeat_info):
    with open(filename, 'w') as f:
        f.write("start_position,end_position,unit_length\n")
        for start, end, unit_length in repeat_info:
            f.write(f"{start},{end},{unit_length}\n")


estimated_total_length = 100000000  # The estimated length of simulated genome/assembly

num_repeats = 50  # 100 repeat regions

repeat_min_length = 50000  # The minimum length of repeat region

repeat_max_length = 1000000  # The maximum length of repeat region

unit_min_length = 40  # The minimum length of repeat unit

unit_max_length = 200  # The maximum length of repeat unit

mutation_rate = 0.01  # The mutation rate for given unit template

dna_sequence, repeat_info = generate_dna_with_repeats(estimated_total_length, num_repeats, repeat_min_length, repeat_max_length, unit_min_length, unit_max_length, mutation_rate)

write_fasta("simulated_genome.fasta", dna_sequence)

write_repeat_info("repeat_regions.csv", repeat_info)

print("Generated DNA sequence and wrote to 'simulated_genome.fasta'.")

print("Generated repeat region information and wrote to 'simulated_genome.fasta'.")

