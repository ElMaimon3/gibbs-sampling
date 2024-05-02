# Necessary imports
import numpy as np
import sys, os
import getopt
import random
import math

def calculate_frequency(sequences):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    frequency = {aa: 0 for aa in amino_acids}
    total = 0
    for sequence in sequences:
        for aa in sequence:
            if aa in frequency:
                frequency[aa] += 1
                total += 1
    for aa in frequency:
        frequency[aa] /= total
    return frequency

def gibbs_sample(infile, temp,decay,number):
    fIn = open(infile)
    sequences = []
    for line in fIn:
        if line[0] != '>':
            sequences.append(line.strip())
    fIn.close()

    # Calculate the frequency of each amino acid
    frequency = calculate_frequency(sequences)

    # Start with random positions
    positions = [random.randint(0, len(seq) - 10) for seq in sequences]

    # Initialize temperature T
    T = temp

    for iteration in range(10000):
        for i in range(len(sequences)):
            # Build a position-specific count matrix
            matrix = {aa: [frequency[aa]] * 10 for aa in frequency}

            # Exclude the current sequence
            for j, seq in enumerate(sequences):
                if i != j:
                    for k in range(10):
                        matrix[seq[positions[j] + k]][k] += 1

            # Calculate the log-likelihood for each entry in the matrix
            pwm = {aa: [np.log2((count / sum(matrix[aa2][pos] for aa2 in frequency)) / frequency[aa]) for pos, count in enumerate(counts)] for aa, counts in matrix.items()}


            # Pick random position
            new_position = random.randint(0, len(sequences[i]) - 10)

            # Calculate the mean log-likelihood of residues in the new window
            old_score = sum(pwm[aa][pos] for pos, aa in enumerate(sequences[i][positions[i]:positions[i] + 10]))
            new_score = sum(pwm[aa][pos] for pos, aa in enumerate(sequences[i][new_position:new_position + 10]))

            # Accept change if likelihood score of new position (Sj) is greater than last position (Si) or with probability e(Sj - Si)/T
            if new_score > old_score or random.random() < math.exp((new_score - old_score) / T):
                positions[i] = new_position

                # Update the count matrix and PWM with the new position
                for k in range(10):
                    matrix[sequences[i][new_position + k]][k] += 1
                pwm = {aa: [np.log2((count / sum(matrix[aa2][pos] for aa2 in frequency)) / frequency[aa]) for pos, count in enumerate(counts)] for aa, counts in matrix.items()}

        # Lower temperature T
        T *= decay
        if iteration % 500 == 0:
            print("Iteration "+str(iteration))

    IC = calculate_IC(pwm)
    n = infile[9]
    outfile1 = f"positions{n}_{temp}_{decay}_{number}.tsv"
    outfile2 = f"pwm{n}_{temp}_{decay}_{number}.tsv"
    write_outputs(outfile1, outfile2, pwm, positions)
    return (old_score,IC)

def write_outputs(outfile1, outfile2, matrix, positions):
    # Write positions output file
    with open(outfile1, 'w') as f:
        for i in range(len(positions)):
            f.write(f'seq{i}\t{positions[i]}\n')
    # Write pwm file
    with open(outfile2, 'w') as f:
        for aa in matrix:
            f.write(f"{aa}\t")
        f.write('\n')
        for i in range(10):
            for aa in matrix:
                f.write(f"{matrix[aa][i]}\t")
            f.write('\n')

def calculate_IC(pwm):
    H = np.log2(20)  # Maximum entropy for amino acids
    IC = 0
    for pos in range(10):
        h = sum(p * np.log2(p) if p > 0 else 0 for aa in pwm for p in [np.power(2,pwm[aa][pos])])
        IC += H + h
    return IC


def main():
    infile = sys.argv[1]
    parameters = sys.argv[2:]
    results = []
    i = 1
    for parameter in parameters:
        temp,decay = map(float,parameter.split(','))
        print(f"Running setting {i}")
        results.append(gibbs_sample(infile, temp,decay,i))
        i += 1
    with open('scores.txt', 'w') as f:
        for r in results:
            f.write(f"score:{r[0]}, IC:{r[1]}\n")

if __name__ == '__main__':
    main()