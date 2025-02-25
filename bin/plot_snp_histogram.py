#!/usr/bin/env python3
import sys
import numpy as np
from matplotlib import pyplot as plt


def load_matrix(MATRIX_FNAME: str):
    """
        Load all vs all SNP distance matrix,
        and return bottom half.
    """
    # Load file. Note first row has number of samples,
    # and last row is empty.
    file_content = open(MATRIX_FNAME, 'r').read().split('\n')[1:-1]
    # matrix_data = dict()
    snp_data = list()
    sample_id = list()
    for row in file_content:
        # sample_id = row.split('\t')[0]
        if len(row.split('\t')) > 1:
            sample_id.append(row.split('\t')[0].replace('barcode',
                                                        '').replace('reference',
                                                                    'Ref.').capitalize())
        data = [int(i) for i in row.split('\t')[1:]]
        # matrix_data[sample_id] = data
        snp_data.append(data)

    # Extract only the triangular/lower section
    # of the distance matrix (includes diagonal)
    idx = np.tril_indices(len(data), -1)
    return np.array(snp_data)[idx], np.matrix(snp_data), tuple(sample_id)


def plot_snp_histogram(matrix_data: np.array):
    """
    """
    fig = plt.figure(figsize=(5, 4))
    plt.hist(matrix_data, bins=50, rwidth=0.95)
    # Labels
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel('No. of SNVs', fontsize=14)
    plt.ylabel('No. of Isolates', fontsize=14)
    plt.savefig('/home/carlos/test1.pdf', format='pdf')
    return


def plot_distance_matrix(snp_data: np.array, sample_id: tuple):
    """
    """
    fig = plt.figure(figsize=(5, 4))
    # Since matrix is symmetric, extract triangular
    # lower section of the matrix (including diag.)
    idx = np.tril_indices(len(snp_data), -1)
    #snp_data[idx] = -1
    snp_alpha = np.ones_like(snp_data, dtype='float')
    snp_alpha[idx] = 0.25
    # Plot
    plt.imshow(snp_data, alpha=snp_alpha, vmin=0, vmax=51, cmap='PuRd')
    ax = plt.gca()  # Get axis
    # Major axis (for labels)
    plt.xticks(ticks=np.arange(len(sample_id)), labels=sample_id,
               rotation=90, fontsize=6)
    plt.yticks(ticks=np.arange(len(sample_id)), labels=sample_id,
               rotation=0, fontsize=6)
    # Minor axis (for grid)
    plt.xticks(ticks=np.arange(-0.5, len(sample_id) - 0.5), minor=True)
    plt.yticks(ticks=np.arange(-0.5, len(sample_id) - 0.5), minor=True)
    # Grid
    ax.grid(which='minor', color='black', linewidth=0.25, alpha=1.0)
    ax.tick_params(which='major', bottom=False, left=False)
    # Labels
    plt.xlabel('Sample ID', fontsize=14)
    plt.ylabel('Sample ID', fontsize=14)
    plt.savefig('/home/carlos/test2.pdf', format='pdf')
    return


def main(FNAME: str):
    """
    
    """
    # Load data
    mx_data, snp_data, sample_id = load_matrix(FNAME)

    # Plot
    plot_snp_histogram(mx_data)
    plot_distance_matrix(snp_data, sample_id)
    plt.show()

if __name__ == '__main__':
    main(sys.argv[1])
