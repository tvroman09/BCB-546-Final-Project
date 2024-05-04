Puja_markdown_file_fig4
## Unix Command
```
mkdir cleaned_data
# Loop through all bd*.txt files except bd1.txt
for file in bd{2..78}.txt; do
    # Remove the header using tail and save the result to a new file in the new_folder
    tail -n +2 "$file" > "new_folder/$file"
done
cp bd1.txt combined.txt
#comined all the bd.txt files and saved into new combined.txt
for file in bd*.txt; do
    if [ "$file" != "bd1.txt" ]; then
        tail -n +2 "$file" >> combined.txt
    fi
done
```
## R script to filter the data with only the targeted 16 species and removing all the other rows.
```
species_to_keep <- c("Sphyraena putnamae", "Scomberomorus commerson", "Scomberoides tol",
  "Saurida tumbil", "Rhabdosargus sarba", "Rastrelliger kanagurta",
  "Psettodes erumei", "Platycephalus indicus", "Nemipterus japonicus",
  "Megalaspis cordyla", "Lethrinus nebulosus", "Epinephelus chlorostigma",
  "Decapterus russelli", "Chirocentrus dorab", "Chanos chanos", "Argyrops spinifer")
filtered_df <- barcode %>% filter(species_name %in% species_to_keep)
write.table(filtered_df, "filtered_combined.txt", sep="\t", quote=FALSE, row.names=FALSE)
```

# Data Inspection. 
## Read the TSV file into a pandas DataFrame
```
import pandas as pd
file_path = "~/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"
try:
    barcode_2 = pd.read_csv(file_path, sep="\t", encoding='utf-8')
except UnicodeDecodeError:

    barcode_2 = pd.read_csv(file_path, sep="\t", encoding='latin-1')
print("Column Names:")
print(barcode_2.columns)
print("\nFirst few rows:")
print(barcode_2.head())
(barcode_2)
```

## Looking for the number of rows for species 'Sphyraena putnamae'
```
import pandas as pd
file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

barcode_2 = pd.read_csv(file_path, sep="\t", encoding='latin1')

# Count the number of rows for the target species
target_species = "Sphyraena putnamae"
num_rows = len(barcode_2[barcode_2['species_name'] == target_species])

print(f"Number of rows for species '{target_species}': {num_rows}")

```
## Similarly, looking for the sequences of the species 'Sphyraena putnamae'
```
from Bio.Seq import Seq

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species name
target_species = "Sphyraena putnamae"

# Initialize a list to store sequences of the target species
target_sequences = []

# Open the TSV file with Latin-1 encoding
with open(file_path, "r", encoding="latin1") as file:
    for line in file:
        cols = line.strip().split("\t")
        species_name = cols[21]  # Assuming species name is in column 22 (0-indexed)
        nucleotides = cols[71]  # Assuming nucleotides are in column 72 (0-indexed)
        if species_name == target_species:
            target_sequences.append(Seq(nucleotides))

# Print the nucleotide sequences of the target species
for seq in target_sequences:
    print(seq)

```
## Calculating the genetic distance for for species 'Sphyraena putnamae' globally.
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations

# Define the file path
file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species name
target_species = "Sphyraena putnamae"

# Initialize a list to store sequences of the target species
target_sequences = []

# Open the TSV file with 
with open(file_path, "r", encoding="latin1") as file:
    for line in file:
        cols = line.strip().split("\t")
        species_name = cols[21]  
        nucleotides = cols[71]  
        if species_name == target_species:
            target_sequences.append(Seq(nucleotides))

# Calculate genetic distances among sequences of the target species
if target_sequences:
    aligner = PairwiseAligner()
    aligner.mode = 'global'  # Use global alignment
    aligner.match_score = 1  # Match score
    aligner.mismatch_score = -1  # Mismatch score
    aligner.open_gap_score = -1  # Open gap penalty
    aligner.extend_gap_score = -1  # Extend gap penalty

    num_sequences = len(target_sequences)
    total_dist = 0
    dists = []

    for seq1, seq2 in combinations(target_sequences, 2):
        alignments = aligner.align(seq1, seq2)
        top_alignment = alignments[0]
        aligned_seq1 = top_alignment[0]  # Get the aligned sequence from the alignment
        aligned_seq2 = top_alignment[1]  # Get the aligned sequence from the alignment
        
        similarity = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2)) / len(aligned_seq1)
        genetic_dist = 1 - similarity
        
        total_dist += genetic_dist
        dists.append(genetic_dist)

    mean_dist = total_dist / (num_sequences * (num_sequences - 1) / 2)
    sd_dist = (sum((dist - mean_dist) ** 2 for dist in dists) / len(dists)) ** 0.5

    print(f"Mean genetic distance: {mean_dist}")
    print(f"Standard deviation of genetic distance: {sd_dist}")
else:
    print(f"No sequences found for species '{target_species}'.")

```
## Calculating the mean and SD of genetic distance for all the 16 target species globally.
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store genetic distances for each species
species_distances = {}

# Open the TSV file 
with open(file_path, "r", encoding="latin1") as file:
    # Read the lines from the file
    lines = file.readlines()

# Parse sequences for each target species
for target_species in target_species_list:
    # Initialize a list to store sequences of the target species
    target_sequences = []

    for line in lines:
        cols = line.strip().split("\t")
        species_name = cols[21]  
        nucleotides = cols[71]  
        if species_name == target_species:
            target_sequences.append(Seq(nucleotides))

    # Calculate genetic distances among sequences of the target species
    if target_sequences:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Use global alignment
        aligner.match_score = 1  # Match score
        aligner.mismatch_score = -1  # Mismatch score
        aligner.open_gap_score = -1  # Open gap penalty
        aligner.extend_gap_score = -1  # Extend gap penalty

        num_sequences = len(target_sequences)
        total_dist = 0
        dists = []

        for seq1, seq2 in combinations(target_sequences, 2):
            alignments = aligner.align(seq1, seq2)
            top_alignment = alignments[0]
            aligned_seq1 = top_alignment[0]  # Get the aligned sequence from the alignment
            aligned_seq2 = top_alignment[1]  # Get the aligned sequence from the alignment
            
            similarity = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2)) / len(aligned_seq1)
            genetic_dist = 1 - similarity
            
            total_dist += genetic_dist
            dists.append(genetic_dist)

        mean_dist = total_dist / (num_sequences * (num_sequences - 1) / 2)
        sd_dist = (sum((dist - mean_dist) ** 2 for dist in dists) / len(dists)) ** 0.5

        # Store mean and standard deviation in the dictionary
        species_distances[target_species] = (mean_dist, sd_dist)
    else:
        species_distances[target_species] = (None, None)  # No sequences found

# Print the mean and standard deviation for each species
for species, (mean, sd) in species_distances.items():
    if mean is not None and sd is not None:
        print(f"Species: {species}")
        print(f"Mean global genetic distance: {mean}")
        print(f"Standard deviation of global genetic distance: {sd}")
        print()
    else:
        print(f"No sequences found for species '{species}'.")

```
## Determining MOTUs for each target species.
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from collections import defaultdict

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the genetic distance threshold for clustering MOTUs
threshold = 0.1  # Adjust this threshold as needed

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store MOTUs for each species
species_motus = defaultdict(list)

# Define the PairwiseAligner
aligner = PairwiseAligner()
aligner.mode = 'global'  # Use global alignment
aligner.match_score = 1  # Match score
aligner.mismatch_score = -1  # Mismatch score
aligner.open_gap_score = -1  # Open gap penalty
aligner.extend_gap_score = -1  # Extend gap penalty

# Open the TSV file with Latin-1 encoding
with open(file_path, "r", encoding="latin1") as file:
    # Read the lines from the file
    lines = file.readlines()

# Parse sequences and cluster MOTUs for each target species
for line in lines:
    cols = line.strip().split("\t")
    species_name = cols[21]  
    nucleotides = cols[71]  

    if species_name in target_species_list:
        # Check if any MOTU contains similar sequences
        found_motu = False
        for motu_seqs in species_motus[species_name]:
            for seq in motu_seqs:
                alignments = aligner.align(Seq(seq), Seq(nucleotides))
                top_alignment = alignments[0]
                aligned_seq1 = top_alignment[0]  # Get the aligned sequence from the alignment
                aligned_seq2 = top_alignment[1]  # Get the aligned sequence from the alignment
                similarity = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2)) / len(aligned_seq1)
                if similarity >= threshold:
                    motu_seqs.append(nucleotides)
                    found_motu = True
                    break
            if found_motu:
                break

        if not found_motu:
            species_motus[species_name].append([nucleotides])

# Count the number of MOTUs for each species
species_motu_counts = {species: len(motus) for species, motus in species_motus.items()}

# Print the number of MOTUs for each species
for species, motu_count in species_motu_counts.items():
    print(f"Species: {species}")
    print(f"Number of MOTUs: {motu_count}")
    print()
```
## Determining the average genetic distance of all the 16 species in Iran.
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store genetic distances for each species within Iran
species_distances_iran = {}

# Open the TSV file with Latin-1 encoding
with open(file_path, "r", encoding="latin1") as file:
    # Read the lines from the file
    lines = file.readlines()

# Parse sequences for each target species within Iran
for target_species in target_species_list:
    # Initialize a list to store sequences of the target species within Iran
    target_sequences_iran = []

    for line in lines:
        cols = line.strip().split("\t")
        species_name = cols[21]  
        nucleotides = cols[71]  
        country = cols[54]  

        if species_name == target_species and country == "Iran":
            target_sequences_iran.append(Seq(nucleotides))

    # Calculate genetic distances among sequences of the target species within Iran
    if target_sequences_iran:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Use global alignment
        aligner.match_score = 1  # Match score
        aligner.mismatch_score = -1  # Mismatch score
        aligner.open_gap_score = -1  # Open gap penalty
        aligner.extend_gap_score = -1  # Extend gap penalty

        num_sequences_iran = len(target_sequences_iran)
        total_dist_iran = 0
        dists_iran = []

        for seq1, seq2 in combinations(target_sequences_iran, 2):
            alignments = aligner.align(seq1, seq2)
            top_alignment = alignments[0]
            aligned_seq1 = top_alignment[0]  # Get the aligned sequence from the alignment
            aligned_seq2 = top_alignment[1]  # Get the aligned sequence from the alignment
            
            similarity = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2)) / len(aligned_seq1)
            genetic_dist = 1 - similarity
            
            total_dist_iran += genetic_dist
            dists_iran.append(genetic_dist)

        mean_dist_iran = total_dist_iran / (num_sequences_iran * (num_sequences_iran - 1) / 2)
        sd_dist_iran = (sum((dist - mean_dist_iran) ** 2 for dist in dists_iran) / len(dists_iran)) ** 0.5

        # Store mean and standard deviation in the dictionary
        species_distances_iran[target_species] = (mean_dist_iran, sd_dist_iran)
    else:
        species_distances_iran[target_species] = (None, None)  # No sequences found within Iran

# Print the mean and standard deviation for each species within Iran
for species, (mean_iran, sd_iran) in species_distances_iran.items():
    if mean_iran is not None and sd_iran is not None:
        print(f"Species: {species}")
        print(f"Mean genetic distance within Iran: {mean_iran}")
        print(f"Standard deviation of genetic distance within Iran: {sd_iran}")
        print()
    else:
        print(f"No sequences found for species '{species}' within Iran.")
```
## Calculating average pair-wise genetic distance for all targeted species in India.
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store genetic distances for each species
species_distances = {}

# Open the TSV file with Latin-1 encoding
with open(file_path, "r", encoding="latin1") as file:
    # Read the lines from the file
    lines = file.readlines()

# Parse sequences for each target species within India (column 54)
for target_species in target_species_list:
    # Initialize a list to store sequences of the target species within India
    target_sequences_india = []

    for line in lines:
        cols = line.strip().split("\t")
        species_name = cols[21]  
        country = cols[54]  
        nucleotides = cols[71]  
        if species_name == target_species and country == "India":
            target_sequences_india.append(Seq(nucleotides))

    # Calculate genetic distances among sequences of the target species within India
    if target_sequences_india:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Use global alignment
        aligner.match_score = 1  # Match score
        aligner.mismatch_score = -1  # Mismatch score
        aligner.open_gap_score = -1  # Open gap penalty
        aligner.extend_gap_score = -1  # Extend gap penalty

        num_sequences = len(target_sequences_india)
        total_dist = 0
        dists = []

        for seq1, seq2 in combinations(target_sequences_india, 2):
            alignments = aligner.align(seq1, seq2)
            top_alignment = alignments[0]
            aligned_seq1 = top_alignment[0]  # Get the aligned sequence from the alignment
            aligned_seq2 = top_alignment[1]  # Get the aligned sequence from the alignment
            
            similarity = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2)) / len(aligned_seq1)
            genetic_dist = 1 - similarity
            
            total_dist += genetic_dist
            dists.append(genetic_dist)

        mean_dist = total_dist / (num_sequences * (num_sequences - 1) / 2)
        sd_dist = (sum((dist - mean_dist) ** 2 for dist in dists) / len(dists)) ** 0.5

        # Store mean and standard deviation in the dictionary
        species_distances[target_species] = (mean_dist, sd_dist)
    else:
        species_distances[target_species] = (None, None)  # No sequences found within India

# Print the mean and standard deviation for each species within India
for species, (mean, sd) in species_distances.items():
    if mean is not None and sd is not None:
        print(f"Species: {species}")
        print(f"Mean genetic distance within India: {mean}")
        print(f"Standard deviation of genetic distance within India: {sd}")
        print()
    else:
        print(f"No sequences found for species '{species}' within India.")
```
## Looking for the number of rows `Argyrops spinifer` in India
```
import pandas as pd

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Read the TSV file into a pandas DataFrame with specified encoding
barcode_2 = pd.read_csv(file_path, sep="\t", encoding='latin1')

# Filter the DataFrame for Argyrops spinifer in India
argyrops_spinifer_india = barcode_2[(barcode_2['species_name'] == 'Sphyraena putnamae') & (barcode_2['country'] == 'India')]

# Count the number of rows
num_argyrops_spinifer_india = len(argyrops_spinifer_india)

print(f"Number of rows for species 'Argyrops spinifer' in India: {num_argyrops_spinifer_india}")

```
## Calculating the average pair-wise genetic distance for all the targeted species in China
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store genetic distances for each species within China
species_distances_china = {}

# Open the TSV file with Latin-1 encoding
with open(file_path, "r", encoding="latin1") as file:
    # Read the lines from the file
    lines = file.readlines()

# Parse sequences for each target species within China
for target_species in target_species_list:
    # Initialize a list to store sequences of the target species within China
    target_sequences_china = []

    for line in lines:
        cols = line.strip().split("\t")
        species_name = cols[21]  
        nucleotides = cols[71]  
        country = cols[54]  

        if species_name == target_species and country == "China":
            target_sequences_china.append(Seq(nucleotides))

    # Calculate genetic distances among sequences of the target species within China
    if target_sequences_china:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Use global alignment
        aligner.match_score = 1  # Match score
        aligner.mismatch_score = -1  # Mismatch score
        aligner.open_gap_score = -1  # Open gap penalty
        aligner.extend_gap_score = -1  # Extend gap penalty

        num_sequences_china = len(target_sequences_china)
        total_dist_china = 0
        dists_china = []

        for seq1, seq2 in combinations(target_sequences_china, 2):
            alignments = aligner.align(seq1, seq2)
            top_alignment = alignments[0]
            aligned_seq1 = top_alignment[0]  # Get the aligned sequence from the alignment
            aligned_seq2 = top_alignment[1]  # Get the aligned sequence from the alignment
            
            similarity = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2)) / len(aligned_seq1)
            genetic_dist = 1 - similarity
            
            total_dist_china += genetic_dist
            dists_china.append(genetic_dist)

        mean_dist_china = total_dist_china / (num_sequences_china * (num_sequences_china - 1) / 2)
        sd_dist_china = (sum((dist - mean_dist_china) ** 2 for dist in dists_china) / len(dists_china)) ** 0.5

        # Store mean and standard deviation in the dictionary
        species_distances_china[target_species] = (mean_dist_china, sd_dist_china)
    else:
        species_distances_china[target_species] = (None, None)  # No sequences found within China

# Print the mean and standard deviation for each species within China
for species, (mean_china, sd_china) in species_distances_china.items():
    if mean_china is not None and sd_china is not None:
        print(f"Species: {species}")
        print(f"Mean genetic distance within China: {mean_china}")
        print(f"Standard deviation of genetic distance within China: {sd_china}")
        print()
    else:
        print(f"No sequences found for species '{species}' within China.")

```
## Determining the average pair-wise genetic distance for all the targeted species in Australia
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store genetic distances for each species within Australia
species_distances_australia = {}

# Open the TSV file with Latin-1 encoding
with open(file_path, "r", encoding="latin1") as file:
    # Read the lines from the file
    lines = file.readlines()

# Parse sequences for each target species within Australia
for target_species in target_species_list:
    # Initialize a list to store sequences of the target species within Australia
    target_sequences_australia = []

    for line in lines:
        cols = line.strip().split("\t")
        species_name = cols[21]  
        nucleotides = cols[71]  
        country = cols[54]  

        if species_name == target_species and country == "Australia":
            target_sequences_australia.append(Seq(nucleotides))

    # Calculate genetic distances among sequences of the target species within Australia
    if target_sequences_australia:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Use global alignment
        aligner.match_score = 1  # Match score
        aligner.mismatch_score = -1  # Mismatch score
        aligner.open_gap_score = -1  # Open gap penalty
        aligner.extend_gap_score = -1  # Extend gap penalty

        num_sequences_australia = len(target_sequences_australia)
        total_dist_australia = 0
        dists_australia = []

        for seq1, seq2 in combinations(target_sequences_australia, 2):
            alignments = aligner.align(seq1, seq2)
            top_alignment = alignments[0]
            aligned_seq1 = top_alignment[0]  # Get the aligned sequence from the alignment
            aligned_seq2 = top_alignment[1]  # Get the aligned sequence from the alignment

            similarity = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2)) / len(aligned_seq1)
            genetic_dist = 1 - similarity

            total_dist_australia += genetic_dist
            dists_australia.append(genetic_dist)

        if num_sequences_australia > 1:
            mean_dist_australia = total_dist_australia / (num_sequences_australia * (num_sequences_australia - 1) / 2)
            sd_dist_australia = (sum((dist - mean_dist_australia) ** 2 for dist in dists_australia) / len(dists_australia)) ** 0.5

            # Store mean and standard deviation in the dictionary
            species_distances_australia[target_species] = (mean_dist_australia, sd_dist_australia)
        else:
            species_distances_australia[target_species] = (None, None)  # Only one sequence found within Australia
    else:
        species_distances_australia[target_species] = (None, None)  # No sequences found within Australia

# Print the mean and standard deviation for each species within Australia
for species, (mean_australia, sd_australia) in species_distances_australia.items():
    if mean_australia is not None and sd_australia is not None:
        print(f"Species: {species}")
        print(f"Mean genetic distance within Australia: {mean_australia}")
        print(f"Standard deviation of genetic distance within Australia: {sd_australia}")
        print()
    else:
        print(f"No sequences found for species '{species}' within Australia.")

```
## Calculating the average pair-wise genetic distance for all 16 targeted species in South Africa.
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store genetic distances for each species within South Africa
species_distances_south_africa = {}

# Open the TSV file with Latin-1 encoding
with open(file_path, "r", encoding="latin1") as file:
    # Read the lines from the file
    lines = file.readlines()

# Parse sequences for each target species within South Africa
for target_species in target_species_list:
    # Initialize a list to store sequences of the target species within South Africa
    target_sequences_south_africa = []

    for line in lines:
        cols = line.strip().split("\t")
        species_name = cols[21]  
        nucleotides = cols[71]  
        country = cols[54]  

        if species_name == target_species and country == "South Africa":
            target_sequences_south_africa.append(Seq(nucleotides))

    # Calculate genetic distances among sequences of the target species within South Africa
    if target_sequences_south_africa:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Use global alignment
        aligner.match_score = 1  # Match score
        aligner.mismatch_score = -1  # Mismatch score
        aligner.open_gap_score = -1  # Open gap penalty
        aligner.extend_gap_score = -1  # Extend gap penalty

        num_sequences_south_africa = len(target_sequences_south_africa)
        total_dist_south_africa = 0
        dists_south_africa = []

        for seq1, seq2 in combinations(target_sequences_south_africa, 2):
            alignments = aligner.align(seq1, seq2)
            top_alignment = alignments[0]
            aligned_seq1 = top_alignment[0]  # Get the aligned sequence from the alignment
            aligned_seq2 = top_alignment[1]  # Get the aligned sequence from the alignment

            similarity = sum(a == b for a, b in zip(aligned_seq1, aligned_seq2)) / len(aligned_seq1)
            genetic_dist = 1 - similarity

            total_dist_south_africa += genetic_dist
            dists_south_africa.append(genetic_dist)

        if num_sequences_south_africa > 1:
            mean_dist_south_africa = total_dist_south_africa / (num_sequences_south_africa * (num_sequences_south_africa - 1) / 2)
            sd_dist_south_africa = (sum((dist - mean_dist_south_africa) ** 2 for dist in dists_south_africa) / len(dists_south_africa)) ** 0.5

            # Store mean and standard deviation in the dictionary
            species_distances_south_africa[target_species] = (mean_dist_south_africa, sd_dist_south_africa)
        else:
            species_distances_south_africa[target_species] = (None, None)  # Only one sequence found within South Africa
    else:
        species_distances_south_africa[target_species] = (None, None)  # No sequences found within South Africa

# Print the mean and standard deviation for each species within South Africa
for species, (mean_south_africa, sd_south_africa) in species_distances_south_africa.items():
    if mean_south_africa is not None and sd_south_africa is not None:
        print(f"Species: {species}")
        print(f"Mean genetic distance within South Africa: {mean_south_africa}")
        print(f"Standard deviation of genetic distance within South Africa: {sd_south_africa}")
        print()
    else:
        print(f"No sequences found for species '{species}' within South Africa.")
```
## Determining the average pair-wise genetic distance of all targeted species between Iran and India.
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations
import pandas as pd

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store genetic distances for each species between Iran and India
species_distances_iran_india = {}

# Read the TSV file into a pandas DataFrame with specified encoding
barcode_df = pd.read_csv(file_path, sep="\t", encoding='latin1')

# Filter rows based on country 
iran_df = barcode_df[barcode_df['country'] == 'Iran']
india_df = barcode_df[barcode_df['country'] == 'India']

# Iterate through each target species
for target_species in target_species_list:
    # Initialize lists to store sequences from Iran and India for the target species
    sequences_iran = iran_df.loc[iran_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()
    sequences_india = india_df.loc[india_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()

    # Calculate genetic distances between Iran and India sequences for the target species
    if sequences_iran and sequences_india:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Use global alignment
        aligner.match_score = 1  # Match score
        aligner.mismatch_score = -1  # Mismatch score
        aligner.open_gap_score = -1  # Open gap penalty
        aligner.extend_gap_score = -1  # Extend gap penalty

        total_dist = 0
        dists = []
        num_sequences = len(sequences_iran)

        for seq_iran, seq_india in zip(sequences_iran, sequences_india):
            alignments = aligner.align(seq_iran, seq_india)
            top_alignment = alignments[0]
            aligned_seq_iran = top_alignment[0]  # Get the aligned sequence from the alignment
            aligned_seq_india = top_alignment[1]  # Get the aligned sequence from the alignment

            similarity = sum(a == b for a, b in zip(aligned_seq_iran, aligned_seq_india)) / len(aligned_seq_iran)
            genetic_dist = 1 - similarity

            total_dist += genetic_dist
            dists.append(genetic_dist)

        mean_dist = total_dist / num_sequences
        sd_dist = (sum((dist - mean_dist) ** 2 for dist in dists) / len(dists)) ** 0.5

        # Store mean and standard deviation in the dictionary
        species_distances_iran_india[target_species] = (mean_dist, sd_dist)
    else:
        species_distances_iran_india[target_species] = (None, None)  # No sequences found for the target species

# Print the mean and standard deviation for each species between Iran and India
for species, (mean_dist, sd_dist) in species_distances_iran_india.items():
    if mean_dist is not None and sd_dist is not None:
        print(f"Species: {species}")
        print(f"Mean genetic distance between Iran and India: {mean_dist}")
        print(f"Standard deviation of genetic distance between Iran and India: {sd_dist}")
        print()
    else:
        print(f"No sequences found for species '{species}' between Iran and India.")
```
## Calculating the interregional average pair wise genetic distance of all targeted species between Iran and China.
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations
import pandas as pd

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store genetic distances for each species between Iran and China
species_distances_iran_china = {}

# Read the TSV file into a pandas DataFrame with specified encoding
barcode_df = pd.read_csv(file_path, sep="\t", encoding='latin1')

# Filter rows based on country 
iran_df = barcode_df[barcode_df['country'] == 'Iran']
china_df = barcode_df[barcode_df['country'] == 'China']

# Iterate through each target species
for target_species in target_species_list:
    # Initialize lists to store sequences from Iran and China for the target species
    sequences_iran = iran_df.loc[iran_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()
    sequences_china = china_df.loc[china_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()

    # Calculate genetic distances between Iran and China sequences for the target species
    if sequences_iran and sequences_china:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Use global alignment
        aligner.match_score = 1  # Match score
        aligner.mismatch_score = -1  # Mismatch score
        aligner.open_gap_score = -1  # Open gap penalty
        aligner.extend_gap_score = -1  # Extend gap penalty

        total_dist = 0
        dists = []
        num_sequences = len(sequences_iran)

        for seq_iran, seq_china in zip(sequences_iran, sequences_china):
            alignments = aligner.align(seq_iran, seq_china)
            top_alignment = alignments[0]
            aligned_seq_iran = top_alignment[0]  # Get the aligned sequence from the alignment
            aligned_seq_china = top_alignment[1]  # Get the aligned sequence from the alignment

            similarity = sum(a == b for a, b in zip(aligned_seq_iran, aligned_seq_china)) / len(aligned_seq_iran)
            genetic_dist = 1 - similarity

            total_dist += genetic_dist
            dists.append(genetic_dist)

        mean_dist = total_dist / num_sequences
        sd_dist = (sum((dist - mean_dist) ** 2 for dist in dists) / len(dists)) ** 0.5

        # Store mean and standard deviation in the dictionary
        species_distances_iran_china[target_species] = (mean_dist, sd_dist)
    else:
        species_distances_iran_china[target_species] = (None, None)  # No sequences found for the target species

# Print the mean and standard deviation for each species between Iran and China
for species, (mean_dist, sd_dist) in species_distances_iran_china.items():
    if mean_dist is not None and sd_dist is not None:
        print(f"Species: {species}")
        print(f"Mean genetic distance between Iran and China: {mean_dist}")
        print(f"Standard deviation of genetic distance between Iran and China: {sd_dist}")
        print()
    else:
        print(f"No sequences found for species '{species}' between Iran and China.")
```
## Calculating the interregional divergence with average pairwise genetic distance between Iran and Australia.
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations
import pandas as pd

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store genetic distances for each species between Iran and Australia
species_distances_iran_australia = {}

# Read the TSV file into a pandas DataFrame with specified encoding
barcode_df = pd.read_csv(file_path, sep="\t", encoding='latin1')

# Filter rows based on country 
iran_df = barcode_df[barcode_df['country'] == 'Iran']
australia_df = barcode_df[barcode_df['country'] == 'Australia']

# Iterate through each target species
for target_species in target_species_list:
    # Initialize lists to store sequences from Iran and Australia for the target species
    sequences_iran = iran_df.loc[iran_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()
    sequences_australia = australia_df.loc[australia_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()

    # Calculate genetic distances between Iran and Australia sequences for the target species
    if sequences_iran and sequences_australia:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Use global alignment
        aligner.match_score = 1  # Match score
        aligner.mismatch_score = -1  # Mismatch score
        aligner.open_gap_score = -1  # Open gap penalty
        aligner.extend_gap_score = -1  # Extend gap penalty

        total_dist = 0
        dists = []
        num_sequences = len(sequences_iran)

        for seq_iran, seq_australia in zip(sequences_iran, sequences_australia):
            alignments = aligner.align(seq_iran, seq_australia)
            top_alignment = alignments[0]
            aligned_seq_iran = top_alignment[0]  # Get the aligned sequence from the alignment
            aligned_seq_australia = top_alignment[1]  # Get the aligned sequence from the alignment

            similarity = sum(a == b for a, b in zip(aligned_seq_iran, aligned_seq_australia)) / len(aligned_seq_iran)
            genetic_dist = 1 - similarity

            total_dist += genetic_dist
            dists.append(genetic_dist)

        mean_dist = total_dist / num_sequences
        sd_dist = (sum((dist - mean_dist) ** 2 for dist in dists) / len(dists)) ** 0.5

        # Store mean and standard deviation in the dictionary
        species_distances_iran_australia[target_species] = (mean_dist, sd_dist)
    else:
        species_distances_iran_australia[target_species] = (None, None)  # No sequences found for the target species

# Print the mean and standard deviation for each species between Iran and Australia
for species, (mean_dist, sd_dist) in species_distances_iran_australia.items():
    if mean_dist is not None and sd_dist is not None:
        print(f"Species: {species}")
        print(f"Mean genetic distance between Iran and Australia: {mean_dist}")
        print(f"Standard deviation of genetic distance between Iran and Australia: {sd_dist}")
        print()
    else:
        print(f"No sequences found for species '{species}' between Iran and Australia.")
```
## Determining the interregional divergence with the pairwise genetic distance of all targeted species beyween Iran and South Africa.
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations
import pandas as pd

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store genetic distances for each species between Iran and South Africa
species_distances_iran_south_africa = {}

# Read the TSV file into a pandas DataFrame with specified encoding
barcode_df = pd.read_csv(file_path, sep="\t", encoding='latin1')

# Filter rows based on country 
iran_df = barcode_df[barcode_df['country'] == 'Iran']
south_africa_df = barcode_df[barcode_df['country'] == 'South Africa']

# Iterate through each target species
for target_species in target_species_list:
    # Initialize lists to store sequences from Iran and South Africa for the target species
    sequences_iran = iran_df.loc[iran_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()
    sequences_south_africa = south_africa_df.loc[south_africa_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()

    # Calculate genetic distances between Iran and South Africa sequences for the target species
    if sequences_iran and sequences_south_africa:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Use global alignment
        aligner.match_score = 1  # Match score
        aligner.mismatch_score = -1  # Mismatch score
        aligner.open_gap_score = -1  # Open gap penalty
        aligner.extend_gap_score = -1  # Extend gap penalty

        total_dist = 0
        dists = []
        num_sequences = len(sequences_iran)

        for seq_iran, seq_south_africa in zip(sequences_iran, sequences_south_africa):
            alignments = aligner.align(seq_iran, seq_south_africa)
            top_alignment = alignments[0]
            aligned_seq_iran = top_alignment[0]  # Get the aligned sequence from the alignment
            aligned_seq_south_africa = top_alignment[1]  # Get the aligned sequence from the alignment

            similarity = sum(a == b for a, b in zip(aligned_seq_iran, aligned_seq_south_africa)) / len(aligned_seq_iran)
            genetic_dist = 1 - similarity

            total_dist += genetic_dist
            dists.append(genetic_dist)

        mean_dist = total_dist / num_sequences
        sd_dist = (sum((dist - mean_dist) ** 2 for dist in dists) / len(dists)) ** 0.5

        # Store mean and standard deviation in the dictionary
        species_distances_iran_south_africa[target_species] = (mean_dist, sd_dist)
    else:
        species_distances_iran_south_africa[target_species] = (None, None)  # No sequences found for the target species

# Print the mean and standard deviation for each species between Iran and South Africa
for species, (mean_dist, sd_dist) in species_distances_iran_south_africa.items():
    if mean_dist is not None and sd_dist is not None:
        print(f"Species: {species}")
        print(f"Mean genetic distance between Iran and South Africa: {mean_dist}")
        print(f"Standard deviation of genetic distance between Iran and South Africa: {sd_dist}")
        print()
    else:
        print(f"No sequences found for species '{species}' between Iran and South Africa.")
```
## Calculating the interregional divergence with the pairwise genetic distance od all 16 targeted species between India and Australia.
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations
import pandas as pd

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store genetic distances for each species between India and Australia
species_distances_india_australia = {}

# Read the TSV file into a pandas DataFrame with specified encoding
barcode_df = pd.read_csv(file_path, sep="\t", encoding='latin1')

# Filter rows based on country 
india_df = barcode_df[barcode_df['country'] == 'India']
australia_df = barcode_df[barcode_df['country'] == 'Australia']

# Iterate through each target species
for target_species in target_species_list:
    # Initialize lists to store sequences from India and Australia for the target species
    sequences_india = india_df.loc[india_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()
    sequences_australia = australia_df.loc[australia_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()

    # Calculate genetic distances between India and Australia sequences for the target species
    if sequences_india and sequences_australia:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Use global alignment
        aligner.match_score = 1  # Match score
        aligner.mismatch_score = -1  # Mismatch score
        aligner.open_gap_score = -1  # Open gap penalty
        aligner.extend_gap_score = -1  # Extend gap penalty

        total_dist = 0
        dists = []
        num_sequences = len(sequences_india)

        for seq_india, seq_australia in zip(sequences_india, sequences_australia):
            alignments = aligner.align(seq_india, seq_australia)
            top_alignment = alignments[0]
            aligned_seq_india = top_alignment[0]  # Get the aligned sequence from the alignment
            aligned_seq_australia = top_alignment[1]  # Get the aligned sequence from the alignment

            similarity = sum(a == b for a, b in zip(aligned_seq_india, aligned_seq_australia)) / len(aligned_seq_india)
            genetic_dist = 1 - similarity

            total_dist += genetic_dist
            dists.append(genetic_dist)

        mean_dist = total_dist / num_sequences
        sd_dist = (sum((dist - mean_dist) ** 2 for dist in dists) / len(dists)) ** 0.5

        # Store mean and standard deviation in the dictionary
        species_distances_india_australia[target_species] = (mean_dist, sd_dist)
    else:
        species_distances_india_australia[target_species] = (None, None)  # No sequences found for the target species

# Print the mean and standard deviation for each species between India and Australia
for species, (mean_dist, sd_dist) in species_distances_india_australia.items():
    if mean_dist is not None and sd_dist is not None:
        print(f"Species: {species}")
        print(f"Mean genetic distance between India and Australia: {mean_dist}")
        print(f"Standard deviation of genetic distance between India and Australia: {sd_dist}")
        print()
    else:
        print(f"No sequences found for species '{species}' between India and Australia.")
```
## Determing the interregional divergence with the average pairwise genetic distance of all the 16 targeted species between China and Australia.
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations
import pandas as pd

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store genetic distances for each species between China and Australia
species_distances_china_australia = {}

# Read the TSV file into a pandas DataFrame with specified encoding
barcode_df = pd.read_csv(file_path, sep="\t", encoding='latin1')

# Filter rows based on country 
china_df = barcode_df[barcode_df['country'] == 'China']
australia_df = barcode_df[barcode_df['country'] == 'Australia']

# Iterate through each target species
for target_species in target_species_list:
    # Initialize lists to store sequences from China and Australia for the target species
    sequences_china = china_df.loc[china_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()
    sequences_australia = australia_df.loc[australia_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()

    # Calculate genetic distances between China and Australia sequences for the target species
    if sequences_china and sequences_australia:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Use global alignment
        aligner.match_score = 1  # Match score
        aligner.mismatch_score = -1  # Mismatch score
        aligner.open_gap_score = -1  # Open gap penalty
        aligner.extend_gap_score = -1  # Extend gap penalty

        total_dist = 0
        dists = []
        num_sequences = len(sequences_china)

        for seq_china, seq_australia in zip(sequences_china, sequences_australia):
            alignments = aligner.align(seq_china, seq_australia)
            top_alignment = alignments[0]
            aligned_seq_china = top_alignment[0]  # Get the aligned sequence from the alignment
            aligned_seq_australia = top_alignment[1]  # Get the aligned sequence from the alignment

            similarity = sum(a == b for a, b in zip(aligned_seq_china, aligned_seq_australia)) / len(aligned_seq_china)
            genetic_dist = 1 - similarity

            total_dist += genetic_dist
            dists.append(genetic_dist)

        mean_dist = total_dist / num_sequences
        sd_dist = (sum((dist - mean_dist) ** 2 for dist in dists) / len(dists)) ** 0.5

        # Store mean and standard deviation in the dictionary
        species_distances_china_australia[target_species] = (mean_dist, sd_dist)
    else:
        species_distances_china_australia[target_species] = (None, None)  # No sequences found for the target species

# Print the mean and standard deviation for each species between China and Australia
for species, (mean_dist, sd_dist) in species_distances_china_australia.items():
    if mean_dist is not None and sd_dist is not None:
        print(f"Species: {species}")
        print(f"Mean genetic distance between China and Australia: {mean_dist}")
        print(f"Standard deviation of genetic distance between China and Australia: {sd_dist}")
        print()
    else:
        print(f"No sequences found for species '{species}' between China and Australia.")
```
## Interregional divergence with the average pair wise genetic distance of all targeted species between China and South Africa.
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations
import pandas as pd

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store genetic distances for each species between China and South Africa
species_distances_china_south_africa = {}

# Read the TSV file into a pandas DataFrame with specified encoding
barcode_df = pd.read_csv(file_path, sep="\t", encoding='latin1')

# Filter rows based on country 
china_df = barcode_df[barcode_df['country'] == 'China']
south_africa_df = barcode_df[barcode_df['country'] == 'South Africa']

# Iterate through each target species
for target_species in target_species_list:
    # Initialize lists to store sequences from China and South Africa for the target species
    sequences_china = china_df.loc[china_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()
    sequences_south_africa = south_africa_df.loc[south_africa_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()

    # Calculate genetic distances between China and South Africa sequences for the target species
    if sequences_china and sequences_south_africa:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Use global alignment
        aligner.match_score = 1  # Match score
        aligner.mismatch_score = -1  # Mismatch score
        aligner.open_gap_score = -1  # Open gap penalty
        aligner.extend_gap_score = -1  # Extend gap penalty

        total_dist = 0
        dists = []
        num_sequences = len(sequences_china)

        for seq_china, seq_safrica in zip(sequences_china, sequences_south_africa):
            alignments = aligner.align(seq_china, seq_safrica)
            top_alignment = alignments[0]
            aligned_seq_china = top_alignment[0]  # Get the aligned sequence from the alignment
            aligned_seq_safrica = top_alignment[1]  # Get the aligned sequence from the alignment

            similarity = sum(a == b for a, b in zip(aligned_seq_china, aligned_seq_safrica)) / len(aligned_seq_china)
            genetic_dist = 1 - similarity

            total_dist += genetic_dist
            dists.append(genetic_dist)

        mean_dist = total_dist / num_sequences
        sd_dist = (sum((dist - mean_dist) ** 2 for dist in dists) / len(dists)) ** 0.5

        # Store mean and standard deviation in the dictionary
        species_distances_china_south_africa[target_species] = (mean_dist, sd_dist)
    else:
        species_distances_china_south_africa[target_species] = (None, None)  # No sequences found for the target species

# Print the mean and standard deviation for each species between China and South Africa
for species, (mean_dist, sd_dist) in species_distances_china_south_africa.items():
    if mean_dist is not None and sd_dist is not None:
        print(f"Species: {species}")
        print(f"Mean genetic distance between China and South Africa: {mean_dist}")
        print(f"Standard deviation of genetic distance between China and South Africa: {sd_dist}")
        print()
    else:
        print(f"No sequences found for species '{species}' between China and South Africa.")
```
## Determining the divergence with the average genetic distance between Australia and South Africa.
```
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from itertools import combinations
import pandas as pd

file_path = "/Users/pujabajracharya/Desktop/Final_project/Data_file/cleaned_data/filtered_combined.txt"

# Define the target species names
target_species_list = [
    "Sphyraena putnamae",
    "Scomberomorus commerson",
    "Scomberoides tol",
    "Saurida tumbil",
    "Rhabdosargus sarba",
    "Rastrelliger kanagurta",
    "Psettodes erumei",
    "Platycephalus indicus",
    "Nemipterus japonicus",
    "Megalaspis cordyla",
    "Lethrinus nebulosus",
    "Epinephelus chlorostigma",
    "Decapterus russelli",
    "Chirocentrus dorab",
    "Chanos chanos",
    "Argyrops spinifer"
]

# Initialize a dictionary to store genetic distances for each species between Australia and South Africa
species_distances_aus_south_africa = {}

# Read the TSV file into a pandas DataFrame with specified encoding
barcode_df = pd.read_csv(file_path, sep="\t", encoding='latin1')

# Filter rows based on country 
aus_df = barcode_df[barcode_df['country'] == 'Australia']
south_africa_df = barcode_df[barcode_df['country'] == 'South Africa']

# Iterate through each target species
for target_species in target_species_list:
    # Initialize lists to store sequences from Australia and South Africa for the target species
    sequences_aus = aus_df.loc[aus_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()
    sequences_south_africa = south_africa_df.loc[south_africa_df['species_name'] == target_species, 'nucleotides'].apply(Seq).tolist()

    # Calculate genetic distances between Australia and South Africa sequences for the target species
    if sequences_aus and sequences_south_africa:
        aligner = PairwiseAligner()
        aligner.mode = 'global'  # Use global alignment
        aligner.match_score = 1  # Match score
        aligner.mismatch_score = -1  # Mismatch score
        aligner.open_gap_score = -1  # Open gap penalty
        aligner.extend_gap_score = -1  # Extend gap penalty

        total_dist = 0
        dists = []
        num_sequences = len(sequences_aus)

        for seq_aus, seq_safrica in zip(sequences_aus, sequences_south_africa):
            alignments = aligner.align(seq_aus, seq_safrica)
            top_alignment = alignments[0]
            aligned_seq_aus = top_alignment[0]  # Get the aligned sequence from the alignment
            aligned_seq_safrica = top_alignment[1]  # Get the aligned sequence from the alignment

            similarity = sum(a == b for a, b in zip(aligned_seq_aus, aligned_seq_safrica)) / len(aligned_seq_aus)
            genetic_dist = 1 - similarity

            total_dist += genetic_dist
            dists.append(genetic_dist)

        mean_dist = total_dist / num_sequences
        sd_dist = (sum((dist - mean_dist) ** 2 for dist in dists) / len(dists)) ** 0.5

        # Store mean and standard deviation in the dictionary
        species_distances_aus_south_africa[target_species] = (mean_dist, sd_dist)
    else:
        species_distances_aus_south_africa[target_species] = (None, None)  # No sequences found for the target species

# Print the mean and standard deviation for each species between Australia and South Africa
for species, (mean_dist, sd_dist) in species_distances_aus_south_africa.items():
    if mean_dist is not None and sd_dist is not None:
        print(f"Species: {species}")
        print(f"Mean genetic distance between Australia and South Africa: {mean_dist}")
        print(f"Standard deviation of genetic distance between Australia and South Africa: {sd_dist}")
        print()
    else:
        print(f"No sequences found for species '{species}' between Australia and South Africa.")

```
## Using Matplotlib to visualize the divergence data in a 2D barplot.
```
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

file_path = "/Users/pujabajracharya/Desktop/Final_project/genetic_dist.xlsx"
df = pd.read_excel(file_path)

# Check the data types of the columns
print(df.dtypes)

# Pivot the dataframe to have countries as columns
pivot_df = df.pivot(index='Species', columns='Country', values='Mean')

# Plotting
plt.figure(figsize=(12, 8))

# Create the bar plot
pivot_df.plot(kind='bar', stacked=True, figsize=(12, 8))

# Customize the plot
plt.title('Mean Genetic Distances Among Species')
plt.xlabel('Species')
plt.ylabel('Mean Genetic Distance')
plt.xticks(rotation=45, ha='right')
plt.legend(title='Country', bbox_to_anchor=(1.05, 1), loc='upper left')

# Show plot
plt.tight_layout()
plt.show()

```
## Reproducing the 3D divergence figure using the matplotlib and mpl_toolkits.mplot3d.
```
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Patch

# Load the Excel file
file_path = "/Users/pujabajracharya/Desktop/Final_project/genetic_dist.xlsx"
df = pd.read_excel(file_path)

# Pivot the dataframe to have countries as columns
pivot_df = df.pivot(index='Species', columns='Country', values='Mean')

# Plotting in 3D
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Define colors for each country
colors = {'global': 'r', 'SA': 'g', 'Ir/SA': 'b', 'Ir/In': 'y', 'Ir/Ch': 'c', 'Ir/Au': 'm', 'Ir': 'orange', 
          'In/Au': 'purple', 'In': 'pink', 'Ch/SA': 'brown', 'Ch/Au': 'olive', 'Ch': 'teal', 'Au/SA': 'lime', 'Au': 'magenta'}

# Create bars for each species
legend_handles = [Patch(facecolor=color, edgecolor='black', label=label) for label, color in colors.items()]
for i, (species, row) in enumerate(pivot_df.iterrows()):
    xpos = i
    ypos = range(len(row))
    zpos = [0] * len(row)
    dx = 0.5
    dy = 0.5
    dz = row.values
    color = [colors[col] for col in row.index]
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=color)

# Customize the plot
ax.set_xlabel('Species', fontsize=12, labelpad=60)  # Adjust fontsize and labelpad as needed
ax.set_ylabel('Country')
ax.set_zlabel('Mean Genetic Distance')
ax.set_xticks(range(len(pivot_df)))
ax.set_xticklabels(pivot_df.index, rotation=45, ha='right', fontsize=9, rotation_mode='anchor')  # Adjust rotation and spacing

# Create legend with custom handles and labels
ax.legend(handles=legend_handles, title='Country', bbox_to_anchor=(1.05, 1), loc='upper left')

# Show plot
plt.title('Global, Regional and interregional divergence between 16 selected fish species')
plt.tight_layout()

# Save the plot with extended bounding box to include the legend
plt.savefig('plot_with_legend.png', bbox_inches='tight')
plt.show()

```