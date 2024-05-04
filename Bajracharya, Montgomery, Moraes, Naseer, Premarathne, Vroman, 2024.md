
# Bajracharya, Montgomery, Moraes, Naseer, Premarathne, Vroman, 2024

# "Cytochrome c oxidase subunit 1 barcode data of fish of the Nayband National Park in the Persian Gulf and analysis using meta-data flag several cryptic species"

Paper citation: ASGHARIAN, H., SAHAFI, H.H., ARDALAN, A.A., SHEKARRIZ, S. and ELAHI, E. (2011), Cytochrome c oxidase subunit 1 barcode data of fish of the Nayband National Park in the Persian Gulf and analysis using meta-data flag several cryptic species. Molecular Ecology Resources, 11: 461-472. https://doi.org/10.1111/j.1755-0998.2011.02989.x

Abstract: We provide cytochrome c oxidase subunit 1 (COI) barcode sequences of fishes of the Nayband National Park, Persian Gulf,Iran. Industrial activities, ecological considerations and goals of The Fish Barcode of Life campaign make it crucial that fishspecies residing in the park be identified. To the best of our knowledge, this is the first report of barcoding data on fishesof the Persian Gulf. We examined 187 individuals representing 76 species, 56 genera and 32 families. The data flaggedpotentially cryptic species of Gerres filamentosus and Plectorhinchus schotaf. 16S rDNA data on these species are provided.Exclusion of these two potential cryptic species resulted in amean COI intraspecific distance of 0.18%, and a mean inter- tointraspecific divergence ratio of 66.7. There was no overlap between maximum Kimura 2-parameter distances among con-specifics (1.66%) and minimum distance among congeneric species (6.19%). Barcodes shared among species were not observed. Neighbour-joining analysis showed that most species formedcohesive sequence units with little variation.Finally, the comparison of 16 selected species from this studywith meta-data of conspecificsfrom Australia, India, Chinaand South Africa revealed high interregion divergences and potential existence of six cryptic species. Pairwise interregional comparisons were more informative than global divergence assessments with regard to detection of cryptic variation. Our analysis exemplifies optimal use of the expanding barcode data now becoming available.

# Biological Question

- Industrial activities, ecological considerations, and the Fish Barcode of Life campaign resulted in need to identify fish species present in Persian Gulf, Iran

- Goal: Create a catalogue of fish species living in Nayband Bay for species identification

- Some species are near threatened and vulnerable (IUCN), catalog will help with rehabilitation

- The Bay located in an energy zone that industrial waste is discharged

- 187 individuals examined; compare 16 selected species to meta-data from Australia, India, China, and South Africa (potential six cryptic species)

![image-2.png](attachment:image-2.png)

# Technical Details of Replication of Analyses

# Figure 2

# R

- In the paper they perform meta data analysis using bold system for sequences of 16 selected fish species and compare these sequences with the sequences from India,China, Australia, and South Africa. They retrieve the sequence information from the 5 projects using the Bold management and analysis system. However, it wasn’t clearly mentioned how we can retrieve that data. So, I retrieve the sequence information for figure 4 using the Bin list under the same project name NNPF.
- Combined all the 77 files and filtered with only the targeted 16 species using Unix command and R.
- Install Bio python using conda and import seq, pairwise aligner and combination using intertool module.
- Calculate the average global, regional, and interregional pairwise genetic distance for all 16 species.
- Then used matplotlib and mpl_toolkits. mplot3d to visualize the data in the 3D figure.

Figure 2 Rabsa Naseer
- Obtained relevant data from BOLD SYSTEMS
    - Utilized public data base and searched for sequences by project name "NPPF"
- Downloaded data as a TSV file in order to obtain relevant metadata
- Renamed sequences to match Figure 2 taxa in Excel
    - Format = species_name|sampleid|family_name
        - If the species_name was unavailable, the family_name was substituted, the paper did the same
    - Saved as a FASTA file using Notepad++
    
For R tree
- Began with reading in sequences as DNA strings 
    - Biostrings package
- Aligned using multiple sequence alignment (msa) package/function
    - Wrote alignment to new file
- Retrieved distrance matrix from alignment file
    - dist.hamming() function in phanghorn package
- Created and viewed neighbor joining tree based on matrix
    - nj function in ape package

# Command Line

For Command Line tree
- Performed MAFFT alignment using FASTA file
    - Used --auto option to automatically pick an alignment method based on data size
- Created tree using iqtree
    - Used GTR model
- Viewed and modified tree in FigTree
    - Swapped branches within clased in order to emulate Figure 2 results
    - Highlighted relevant taxa

# BOLD

- Search and open BOLD Systems in the browser
- Log in to the system or create an account if you are a first-time user
- Search for the NNPG project as the code
- Then select sequence analysis and select Taxon ID tree
- In pop-up window, select:
    - Tree type: multipage classic
    - Align sequences: None (use submitted alignment)
    - Select preferred labels, in this case Sequence/Process ID, Species and order
    - Ambiguous base/Gap handling: Pairwise deletion
    - Codon positions include: 1st, 2nd, 3rd
- Click: build tree
- From the next pop-up window, download the tree file and record the analysis description file

- See Figure_2_BOLD_description.pdf for attached images

# Figure 3

# MEGA

Steps to get the sequences for Figure 3:

- Obtain the sequences for Figure 3:
    - There were 14 sequences in Figure 3_ In the Supporting information from the manuscript, there is a file called “MEN_2989_sm_tS2.xls”
    - In this file was possible to obtain only 6 sequences at the column named “Genbank Num” by looking at the “Sample ID” described in the phylogenetic tree in Figure 3.
    - These 6 sequences were from Argyrops spinifer specimens from Iran.
    - Then we searched in NIH in search each gen bank number https://www.ncbi.nlm.nih.gov/. For example, for Sample ID “NPPF1048” we entered “HQ149794”.
    - For the 8 remaining sequences, we entered each Sample ID in gen Bank. For example, for Sample ID BW-A675, we entered “Argyrops spinifer BW-A675”
    - Then we saved the 14 sentences from FASTA format in a text file called “BCB 546 final project-sequences-Figure 3.
- Import the sequences to MEGA- Molecular Evolutionary Genetic Analysis
    - The version from the manuscript was MEGA 4.0
    - The current version that we used was MEGA 11
- The alignment was conducted by MUSCLE
- To build the phylogenetic tree: Neighbour-joining (NJ) tree

# Figure 4

# BOLD

- The tables were created using data found on BOLDSYSTEMS
- Using Workbench on BOLDSYSTEMS, enter NNPF under code
- Press Sequence Analysis and Distance Summary
- Use Kimura 2 Parameters
- Summarized in table:

![image.png](attachment:image.png)

- For Table 3, repeat the process from table 1 using codes AUSA, WLIND, FSCS, MM, TZFPC, and TZAIC
- Code MM does not retrieve the correct data. For MM, you have to use MEFM and MXII, these two data sets are combined using Project Options, Merge Projects. Then they are treated like the other files.
- The data is summarized in the table below:

![image.png](attachment:image.png)

- A suggestion was made that the final figure could be made in excel. I copied the data presented on table 2 into an excel file and then attempted to use excel to create the Figure_4_Excel.png

# R 

- Following information found in: FIgure_4_Markdown.md

- Use UNIX to filter data
- R script to filter the data with only the targeted 16 species and removing all the other rows
- Data inspection: Read the TSV file into a pandas DataFrame
- Look for the number of rows for species ‘Sphyraena putnamae’
- Similarly, look for the sequences of the species ‘Sphyraena putnamae’
- Calculate the genetic distance for for species ‘Sphyraena putnamae’ globally
- Calculate the mean and SD of genetic distance for all the 16 target species globally
- Determine MOTUs for each target species
- Determine the average genetic distance of all the 16 species in Iran
- Calculate average pair-wise genetic distance for all targeted species in India
- Look for the number of rows Argyrops spinifer in India
- Calculate the average pair-wise genetic distance for all the targeted species in China
- Determine the average pair-wise genetic distance for all the targeted species in Australia
- Calculate the average pair-wise genetic distance for all 16 targeted species in South Africa
- Determine the average pair-wise genetic distance of all targeted species between Iran and India
- Calculate the interregional average pair wise genetic distance of all targeted species between Iran and China
- Calculate the interregional divergence with average pairwise genetic distance between Iran and Australia
- Determine the interregional divergence with the pairwise genetic distance of all targeted species beyween Iran and South Africa
- Calculate the interregional divergence with the pairwise genetic distance od all 16 targeted species between India and Australia
- Determine the interregional divergence with the average pairwise genetic distance of all the 16 targeted species between China and Australia
- Interregional divergence with the average pair wise genetic distance of all targeted species between China and South Africa
- Determine the divergence with the average genetic distance between Australia and South Africa
- Use Matplotlib to visualize the divergence data in a 2D barplot
- Reproduce the 3D divergence figure using the matplotlib and mpl_toolkits.mplot3d

# Summarization

We were able to produce the figures to come degree, though several issues, outlined below, made it very challenging. Reproducing figure 2 using methods not described in the paper went well. However, we believe this was due to prior knowledge and experience of sequence alignment and tree building. The methods we used, though simple, would not be well suited for someone inexperienced in bioinformatics attempting to perform analyses described in the paper. Figure 4 that was generated using python was not same as in the figure. This might bedue to the differences in the sequence data used in the paper and in our project as it was unclear. It was also not clearly mentioned about which software they used to make the Figure 4. They also didn’t mention about global divergence whether they did that for only 5 countries or also with the sequences of the same species from other countries in the data file.

# Primary Issues

- In the paper, they talked about the sequence information from the 5 projects but didn't clearly mention where we can retrieve that information.
- It was also unclear about which software they used for the visualization of the divergence data.
- Aligning the format of the data – data files from various sources and not in one cohesive document
- The information was not clear and it was incomplete in the supplemental material
- Parameters were not easily accessible
- Limited permissions on the files created issues as well

