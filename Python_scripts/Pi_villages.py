import os
import numpy as np
from metrics import alleles_count, Pi
import glob
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Compute mean nucleotide diversity for each chromosome and each replicate')
    parser.add_argument('-s', '--source', dest='path_source', required=True, help='Source folder')
    parser.add_argument('-c', '--chromosome-size', dest='chr_size', default='1e6,1e6,1e6,1e4', help='Vector of autosome, X chromosome, Y chromosome and mtDNA sizes')
    parser.add_argument('-ph', '--pseudohap', dest='pseudohap', choices=('True', 'False'), default='False', help='Use a pseudo-haploid genome?')
    parser.add_argument('-o', '--output', dest='output', help='Output folder path')
    return parser.parse_args()

def write_header(file_path, header):
    with open(file_path, 'w') as f_out:
        print('\t'.join(header), file=f_out)

def write_line(file_path, line):
    with open(file_path, 'a') as f_out:
        print('\t'.join(map(str, line)), file=f_out)

def process_files(filenames, chr_size, pseudohap):
    pi_list = []
    for file in filenames:
        A1, B1, nb_chr, nb_snps = alleles_count(file, pseudohap=pseudohap)
        Pi_val = Pi(chr_size, A1, B1)
        pi_list.append(Pi_val)
    return np.array(pi_list)

def compute_means(rep, gen, chr_size, pseudohap, path_source, output):
    chr_size_Y, chr_size_mito = chr_size
    rep_path = os.path.join(path_source, str(rep))

    filenames_Y = sorted(glob.glob(os.path.join(rep_path, f'Sim_Y_*_{rep}_gen_{gen}_*.vcf')))
    filenames_Mito = sorted(glob.glob(os.path.join(rep_path, f'Sim_Mito_*_{rep}_gen_{gen}_*.vcf')))

    pi_Y, theta_Y, tajima_Y, snps_Y = process_files(filenames_Y, chr_size_Y, pseudohap)
    pi_Mito, theta_Mito, tajima_Mito, snps_Mito = process_files(filenames_Mito, chr_size_mito, pseudohap)

    metrics = [
        ('Pi_Y_mean_by_rep.txt', pi_Y, theta_Y, tajima_Y, snps_Y),
        ('Pi_Mito_mean_by_rep.txt', pi_Mito, theta_Mito, tajima_Mito, snps_Mito)
    ]

    os.chdir(output)
    for file_path, pi, theta, tajima, snps in metrics:
        line = [rep, gen, np.nanmean(pi), np.nanmean(theta), np.nanmean(tajima), np.nanmean(snps)]
        write_line(file_path, line)

def main():
    args = parse_args()
    path_source = args.path_source
    chr_size = [float(size) for size in args.chr_size.split(',')]
    pseudohap = str(args.pseudohap)
    output = args.output

    os.chdir(output)
    output_files = ['Pi_X_mean_by_rep.txt',
        'Pi_Y_mean_by_rep.txt', 
        'Pi_A_mean_by_rep.txt', 
        'Pi_Mito_mean_by_rep.txt'
    ]
    header = ['mf/m', 'Gen', 'Pi', 'Theta', 'tajima_D', 'Nb_SNPs']
    for out_file in output_files:
        write_line(out_file, header)
    
    generations = [0, 20, 40, 60, 80, 100]
    
    for rep in range(1, 101):
        print(rep)
        for gen in generations:
            compute_means(rep, gen, chr_size, pseudohap, path_source, output)

if __name__ == "__main__":
    main()