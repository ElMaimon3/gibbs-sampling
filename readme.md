# Gibbs Sampling for Protein Sequences
This Python script performs Gibbs sampling on a set of protein sequences to find regions of high similarity, which can indicate functionally important motifs.

# Usage
The script can be run from the command line with the following syntax:

python gibbs_sampling.py <input_file> <settings>...

<input_file>: A FASTA file containing the protein sequences.
<settings>: Any number of settings in the format {T},{decay}. Each setting will run the Gibbs sampling with the specified initial temperature T and decay factor decay.
For example, to run the algorithm on the second set of sequences, with T=60 and a decay of 0.99925, doing 4 runs, the command would be:

python gibbs_sampling.py sequences2.fa 60,0.99925 60,0.99925 60,0.99925 60,0.99925

# Output
The script saves a positions file and a PWM (Position Weight Matrix) file for each run. The positions file contains the final positions of the motifs in each sequence. The PWM file contains the PWM representing the motif.

# Visualizing Position Weight Matrices with Sequence Logos
The makelogo.py script allows you to visualize Position Weight Matrices (PWMs) as sequence logos.
Please note that you need to have the logomaker library installed. You can install it using pip:
pip install logomaker
Make sure to change the filename to the appropriate PWM file from the gibbs sampling algorithm.
This script will generate a sequence logo from the PWM file, which provides a visual representation of the motif found by the Gibbs sampling algorithm. The height of each letter in the sequence logo corresponds to the frequency of the corresponding amino acid at that position in the motif. This can help you to understand the properties of the motif and to compare it with other motifs.

# Pseudocode
Here’s a high-level overview of the algorithm implemented in the script:

Procedure calculate_frequency(sequences)
  Initialize frequency dictionary with amino acids as keys and 0 as values
  For each sequence, increment the count of each amino acid in frequency
  Normalize frequency by total count of all amino acids
  Return frequency
End Procedure

Procedure gibbs_sample(infile, temp, decay, number)
  Read sequences from infile, ignoring lines starting with '>'
  Calculate frequency of amino acids in sequences
  Initialize positions with random integers for each sequence
  Initialize temperature T

  For a set number of iterations
    For each sequence
      Initialize matrix with frequency of each amino acid at each position
      Update matrix based on other sequences' current positions
      Calculate position-weight matrix (pwm) from matrix
      Choose new position for current sequence based on pwm
      If new position improves score or meets a random threshold, update position and pwm
    Reduce T by decay factor
    Print iteration number every 500 iterations

  Calculate information content (IC) from final pwm
  Write positions and pwm to output files
  Return final score and IC
End Procedure

Procedure write_outputs(outfile1, outfile2, matrix, positions)
  Write positions to outfile1 with each line in the format "seq[i]\tpositions[i]"
  Write pwm to outfile2 with each column containing an amino acid and its frequencies at each position
End Procedure

# Dependencies
This script requires Python 3 and the following Python libraries:

numpy
sys
os
getopt
random
math
# Limitations
This script assumes that the input file and parameters are provided in the correct format and that the input file exists and can be read. The script also assumes that the sequences in the input file are long enough for the motif length and that the sequences contain only the 20 standard amino acids. The script does not check for these conditions and may produce incorrect results or crash if the conditions are not met.

# Future Improvements
Future versions of this script could include error handling and input validation, performance optimization, parallelization, checkpointing and resuming, progress monitoring and visualization, result analysis and interpretation, comparison with other methods, validation of results, handling of biases and errors in the input data, handling of uncertainties and variabilities in the sequences, handling of the biological context or function of the sequences, handling of the evolutionary relationships or conservation of the sequences, handling of the structural or physicochemical properties of the sequences, handling of the experimental or computational methods used to obtain the sequences, handling of the limitations or assumptions of the Gibbs sampling or the PWM, handling of the limitations or assumptions of the scoring system or the IC, handling of the limitations or assumptions of the temperature and decay parameters, handling of the limitations or assumptions of the random number generation, handling of the limitations or assumptions of the numerical computations, handling of the limitations or assumptions of the file formats or the command-line arguments.
