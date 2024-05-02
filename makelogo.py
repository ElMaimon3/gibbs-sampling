import pandas as pd
import logomaker
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

with open('pwm4_60.0_0.99925_1.tsv', 'r') as f:
    lines = f.readlines()

# Get the amino acids from the first line
amino_acids = lines[0].strip().split('\t')

# Get the log likelihoods from the rest of the lines
log_likelihoods = [list(map(float, line.strip().split('\t'))) for line in lines[1:]]

# Create a dictionary with amino acids as keys and log likelihoods as values
data = {amino_acids[i]: [log_likelihoods[j][i] for j in range(len(log_likelihoods))] for i in range(len(amino_acids))}

df = pd.DataFrame(data)
df2 = df.applymap(np.exp2)
# Convert the PWM to a logomaker.Logo object
logo = logomaker.Logo(df2)

# Display the sequence logo
logo.draw()
plt.show()