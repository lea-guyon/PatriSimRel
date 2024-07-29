# README

Pipeline simulating human genomes under scenarios with patrilocal residence and different descent systems (bilateral, patrilineal), and computing relatedness coefficients between males and between females, as well as diversity estimators on uniparental genetic markers.

## Installation
`git clone https://github.com/lea-guyon/PatriSimRel`

## Summary

To run the pipeline, modify the `PARAMETERS` section in the file `main.sh`. See "Parameters information" to have a description of each parameter. Then execute the .sh file with bash in a terminal. 

`bash main.sh`

This pipeline:
1) simulates autosomes, X chromosomes, Y chromosomes and mtDNA under the scenario specified by the parameters using the sofware SLiM (see "SLiM models"). 
2) converts SLiM outputs into VCF files with custom python scripts (see "Output VCF")
3) computes nucleotide diversity using custom python scripts (see "Nucleotide diversity")

NB: if you already ran the script once, it is recommended to set the `burnin` parameter to "false".

To run the scenarios presented in the paper, use the parameters values indicated in `scenarios_parameters.md`

## Requirements

### System requirements
The pipeline requires a standard computer with enough RAM to support the in-memory operations.

#### OS requirements
The pipeline is designed to be run on Linux. It has been tested on Ubuntu 20.04.

### SLiM
The pipeline requires SLiM v4.1. To install SLiM, follow the instructions given in the manual (https://messerlab.org/slim/). 

### Python dependencies
To install python (3.8) requirements enter:  
`pip3 install -r requirements.txt`

## Parameters information

`dir`: input path

`chr_size`: vector of chromosome sizes. First element: autosomes and X. Second element: Y chromosome. Third element: mtDNA.

`random_fission`: true or false. If true, the fission type is random.

`transmission`: "full" or "half". If "full", the growth rate of a splitting descent group is transmitted to the two resulting descent groups after the fission event. If "half" it is transmitted to one out of the two new groups.

`fission_threshold`: number of individuals from which a descent group splits

`pM`: float between 0 and 1. pM is the probability for the smallest group resulting from a split to move to another village.

`descent`: "unilineal" or "bilateral"

`descent_rule`: "patrilineal" or "matrilineal"

`nb_villages`: number of villages in the population

`nb_groups`: number of descent groups -> does not make sense for bilateral descent but usefull to normalize villages' size

`K`: initial number of individuals per descent group

`polygyny`: "F" or "T". If "T", males mate with several women (change the slim code to modify the rate of polygyny)

`mf`: proportion of females migrating to a new village at each generation.

`mm`: proportion of males migrating to a new village at each generation.

`sigma`: variance of the normal law used to draw relative fitnesses

`growth_rate`: growth rate of villages and outgroup, if 0 : population has a constant size

`sample_size`: number of individuals to sample in each village

`nbsimu`: number of replicates

`cores`: number of cores

`nameDir`: name of the output directory

## SLiM models
*Contains .slim files modeling populations with different descent (patrilineal, matrilineal, bilateral) and post-marital residence rules (patrilocal, matrilocal)*

- **burnin.slim**: simulate a panmictic population made of `K_total` individuals for 20,000 generations. This preliminary simulation is taken in entry of the models of villages with bilateral or unilineal descent. 

- **bilateral_descent_relatedness.slim**: simulates autosomes, X, Y chromosomes and mtDNA in `N` villages, each containing `K` individuals at the beginning of the simulation for 100 generations. A fraction of females `mf` and a fraction of males `mm` migrate between villages. Also computes relatedness coefficients between males and between females within villages or local groups.

- **unilineal_descent_relatedness.slim**: simulates autosomes, X, Y chromosomes and mtDNA in `N` villages, each containing `K` individuals divided into `nb_groups` descent groups at the beginning of the simulation for 100 generations. A fraction of females `mf` and a fraction of males `mm` migrate between villages. The descent rule is defined by `descent_rule`. Also computes relatedness coefficients between males and between females within villages or local groups.

- **patrilineal2bilateral_relatedness.slim**: simulates autosomes, X, Y chromosomes and mtDNA in `N` villages, each containing `K` individuals divided into `nb_groups` descent groups at the beginning of the simulation for 100 generations. A fraction of females `mf` and a fraction of males `mm` migrate between villages. The descent rule is defined by `descent_rule`. Then there is a transition to a bilateral system. This continues for another 100 generations. Also computes relatedness coefficients between males and between females within villages or local groups. This script can be executed by running **main_patrilineal2bilateral.sh**.

## Output VCF
*.py scripts taking .trees files (generated by simulating populations in SLiM) in entry to produce .vcf files for each chromosome and each village in every simulation*

In `Python_scripts` folder:
- **subset_trees_bilateral_descent_relatedness.py**
- **subset_trees_unilineal_descent_relatedness.py**

## Nucleotide diversity
*.py files used to compute diversity metrics (Pi) from .vcf files*

In `Python_scripts` folder:
- **metrics.py** : file containing functions to compute alleles frequencies and diversity metrics

- **Pi_villages.py** : compute nucleotidic diversity for models of villages

In `R_scripts` folder:

-**figures.R**: plots mean relatedness between males and between females and male and female effective sizes for each scenario.
