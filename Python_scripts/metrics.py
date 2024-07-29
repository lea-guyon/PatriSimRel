import random
import numpy as np

def alleles_count(vcf_file, pseudohap = False) :

    """
    Description : Compute the number of alleles 0 and the number of alleles 1
    for each SNP in a sampled population or village

    Arguments : vcf_file = VCF file for a given chromosome, a given replicate and a given generation
                chr_size = size of chromosomes
                pseudohap = boolean indicating if genomes are pseuhaploidized (default = False)

    Return : A = dictionnary associating each SNP with the number of alleles 0 in the sample
            B = dictionnary associating each SNP with the number of alleles 1 in the sample
            nb_chr = nb of chromosomes in the sample
    """
    A = {}
    B = {}
    nb_chr = 0
    nb_snps = 0

    with open(vcf_file, 'r') as f1:
        for line in f1:
            # skip header
            if line.startswith('#'):
                continue

            line = line.strip().split('\t')
            snp_info = line[9:]  # select only SNP information
            snp_id = line[1]

            if pseudohap == 'True':
                print('pseudo')
                coef_A = random.choices([0, 1], k=4)
                coef_B = [1 - i for i in coef_A]
                A_snp = float(snp_info.count('0|0') + sum(coef_A[i] * snp_info.count(j) for i, j in enumerate(['0|1', '1|0', '2|0', '0|2'])) + snp_info.count('0'))
                B_snp = float(snp_info.count('1|1') + snp_info.count('2|2') + snp_info.count('2|1') + snp_info.count('1|2') + sum(coef_B[i] * snp_info.count(j) for i, j in enumerate(['0|1', '1|0', '2|0', '0|2'])) + snp_info.count('1') + snp_info.count('2'))
            else:
                snp_info_str = '\t'.join(snp_info)
                A_snp = float(snp_info_str.count('0'))  # count number of 0
                B_snp = float(snp_info_str.count('1') + snp_info_str.count('2'))  # count number of 1 (and 2)

            A[snp_id] = A_snp
            B[snp_id] = B_snp

            if nb_chr == 0:
                nb_chr = A_snp + B_snp

            if 0 < A_snp < nb_chr:
                nb_snps += 1  # count number of polymorphic SNPs

    return(A, B, nb_chr, nb_snps)

def Pi(chr_size, A, B) :
    """
    Description : compute nucleotide diversity (ie expected heterozygosity)

    Arguments : chr_size = size (in nucleotides) of the modeled chromosomes
                A = dictionnary associating each SNP with the number of alleles 0 in the sample
                B = dictionnary associating each SNP with the number of alleles 1 in the sample

    Return : pi = mean nucleotide diversity over whole chromosomes
    """
    a_values = np.array(list(A.values()))
    b_values = np.array(list(B.values()))

    total_alleles = a_values + b_values
    valid_indices = total_alleles > 1

    a_valid = a_values[valid_indices]
    b_valid = b_values[valid_indices]
    total_valid = total_alleles[valid_indices]

    pi_sum = np.sum(a_valid * b_valid / (0.5 * (total_valid - 1) * total_valid))
    pi = pi_sum / float(chr_size)

    return(pi)