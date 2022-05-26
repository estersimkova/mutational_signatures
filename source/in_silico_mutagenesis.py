# imports
import json
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path  # to save the mutational table 


# the different paths to the files
PATH_TO_HUMAN_COUNTS = "../data/codon_counts_GRCh37.json"
PATH_TO_SIGNATURES = "../data/external/COSMIC_v3.2_SBS_GRCh37.txt"
PATH_TO_SIGNATURES_NORMALIZED = "../data/COSMIC_v3.2_SBS_GRCh37_normalized.txt"
PATH_TO_REF_SEQ = "../data/sars-ref-seq.txt"

# 1. Normalise the COSMIC signatures to human trinucleotide content

# code from normalize_cosmic_signatures.py

def load_triplet_counts(path: str):
    """ read and collapse raw trinucleotide counts """
    with open(path) as fin:
        counts = json.load(fin)

    new_counts = defaultdict(int)
    for trinuc, num in counts.items():
        standart_trinuc = trinuc.upper()
        if len(set(standart_trinuc).difference("ATGC")) == 0:
            new_counts[standart_trinuc] += num
    return new_counts


def load_reference_sig(path: str):
    df = pd.read_csv(path, sep="\t")
    return df


def normalize(triplets, signatures):
    assert "Type" in signatures.columns
    ref_triplet = signatures.Type.apply(lambda s: "".join([s[0], s[2], s[-1]]))
    counts_raw = ref_triplet.map(triplets)
    counts = counts_raw / counts_raw.sum()

    signatures = signatures.set_index("Type")
    signatures_normed = signatures.div(counts.values, axis=0)
    signatures_normed = signatures_normed.div(
        signatures_normed.sum(axis=1).values, axis=0)
    return signatures_normed.reset_index()


# normalise the signatures to human trinucleotide content
triplets = load_triplet_counts(PATH_TO_HUMAN_COUNTS)
signatures = load_reference_sig(PATH_TO_SIGNATURES)
    
signatures_normed = normalize(triplets, signatures)
signatures_normed.to_csv(PATH_TO_SIGNATURES_NORMALIZED, index=None, sep="\t") # output path

# assuming symmetry, we create 192 probabilities from the 96 originally present by taking the complementary context
signatures_normed_copy = signatures_normed.copy()

def replace_nucleotides(s):
    # function that replaces the type of mutation with the complementary
    return s.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()

signatures_normed_copy["Type"] = signatures_normed_copy["Type"].apply(replace_nucleotides)

# we do the reverse complement now, 
def shift_nucleotides(s):
    # function that shifts the first and last nucleotide of the type column
    if len(s) <= 1:
        return s
 
    mid = s[1:len(s) - 1]
    return s[len(s) - 1] + mid + s[0]

signatures_normed_copy["Type"] = signatures_normed_copy["Type"].apply(shift_nucleotides)

signatures_normed = pd.concat([signatures_normed, signatures_normed_copy], ignore_index=True)

# get the nucleotides of the Type column
signatures_normed["Context"] = signatures_normed.Type.str.slice(start=0, stop=7, step=2) 

# drop the third value to get the trinucleotide context before each mutation 

for i, strObj in enumerate(signatures_normed["Context"]):
    index = 2
    if len(strObj) > index:
        strObj = strObj[0 : index : ] + strObj[index + 1 : :]
    signatures_normed["Context"][i] = strObj
    
    
# Get a chosen signature, here SBS20 = symmetric and not that common, SBS25 = more uniform, SBS34 = very caracteristic and different from the others 

#SBS = "SBS20"
#SBS = "SBS25"
SBS = "SBS34"

sbs = signatures_normed.loc[:,["Type", SBS, "Context"]]
sbs = sbs.rename(columns={SBS: "SBS"})

# 2. Opening the SARS-cov2 ref seq 
with open(PATH_TO_REF_SEQ) as f:
    ref_seq_raw = f.readlines()

ref_seq = " ".join(ref_seq_raw) # concatenate the string
ref_seq.replace(" ", "")

lst = list(ref_seq) # convert to vector of char to be able to sample fast

buffer = []
for c in lst:
    if c != "\n" and c != ' ':
        buffer.append(c)
    

# 3. Mutate the SARS-cov2 reference sequence with respect to the chosen signature
# Do the loop until we get n mutations

n = 1000
i = 0

# create the table in which to add the mutations
mutation_table = sbs.loc[:,["Type"]]
mutation_table["Number"] = 0


while i < n:
    # Sample a random nucleotide from the sequence and get the context
    seq = buffer
    nt = np.random.randint(1, len(buffer)-1, size=1)

    list_context = [seq[nt[0]-1],seq[nt[0]],seq[nt[0]+1]]
    context = ''.join(list_context)

    # Check the probabilities in the sbs signature
    subdf = sbs.loc[sbs['Context'] == context]

    # Draw random number between 0 and 1
    rand_nb = np.random.random_sample()

    # Sum all three mutational probabilities, if higher than rand_nb: mutate, otherwise: not
    sum_probas = subdf['SBS'].values.sum()
    mutate = False

    if sum_probas > rand_nb:
        mutate = True
    else:
        mutate = False

    # Get the normalised probas, that sum up to 1
    subdf_normed = subdf.copy()
    subdf_normed['SBS']= subdf['SBS'].values/sum_probas

    if mutate == True: # so if proba bigger than random number
        # Choose randomly which type of mutation it will be
        chosen_mut = subdf_normed.sample(n=1, weights='SBS')
        # Extract the nucleotide to mutate into
        mut_nt = chosen_mut.Type.str.slice(start=4, stop=5, step=1)
        i = i+1
        print(i, "out of", n)
        
        # Write the mutation type
        mutation = seq[nt[0]-1] + "[" + seq[nt[0]] + ">" + mut_nt.values[0] + "]" + seq[nt[0]+1]
        #print(mutation)
        
        # Add one to count if this mutation occured
        mutation_table["Number"].loc[mutation_table['Type'] == mutation] += 1

        # Mutate the reference sequence
        buffer[nt[0]] = mut_nt.values[0]

        
# Print the mutation table obtained and save it to tsv format to be opened in the jupyter notebook        
print("The mutation spectrum obtained is the following:", mutation_table)
filepath = Path('../results/results/{}_mut/mut_table_{}.tsv'.format(n,n))  
filepath.parent.mkdir(parents=True, exist_ok=True)  
mutation_table.to_csv(filepath, sep="\t") 

# 4. Save the resulting sequencing in another file 

f = open("../results/results/{}_mut/results_{}.txt".format(n,n), "w")
string_output = ' '.join(map(str, buffer))
f.write(string_output)
f.close()