#!/usr/bin/env python3

import gzip as gz
import numpy as np
import matplotlib.pyplot as plt 
import sys
import os


def _calculate_error(PHRED: int) -> float:
    """
        Phred score to error probability
    """
    return 10 ** (-PHRED/10)


def _calculate_score(ASCII_CODE: str) -> int:
    """
        Given a certain score, in ASCII format, return
        Phred-scale score.

        Because Q-score = ASCII code - 33, it
        coincidentally equals the index of the ASCII
        characters below with ! being 0, and ~ the
        highest (93)

        Note: Many parsers cap the highest at 60 
              (1 err per 1,000,000 nt)
    """
    Q_SCORE = str('!"#$%&') +  str("'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~")
    return _calculate_error(Q_SCORE.index(ASCII_CODE))


def _retrieve_scores(FASTQ: str) -> list:
    """
        Load up reads from FASTQ file, and retrieve the score (in ASCII format)
        assigned to each nucleotide within the read. Thus,
        len(score) == len(read)
    """
    from Bio import SeqIO

    RAW_DATA = gz.open(FASTQ, 'rt')  # Assumes fastq.gz
    reads = list()
    for read in SeqIO.parse(RAW_DATA, 'fastq'):
        reads.append(read.letter_annotations['phred_quality'])
    return reads


def plot_data(axis: plt.Axes, READS_QUAL: np.ndarray, THRESHOLD_QUAL: int):
    """
       Produce a box/scatter plot of data qualilty vs err.
    """
    # From QUAL to ACCURACY - NB: all data clustered around 99%
    #READS_ACC = 1 - _calculate_error(READS_QUAL)
    #THRESHOLD_ACC = 1 - _calculate_error(THRESHOLD_QUAL)

    # Box plot
    axis.boxplot(READS_QUAL[READS_QUAL >= THRESHOLD_QUAL],
                 positions=tuple([THRESHOLD_QUAL]), widths=1.5)
    # Scatter plot
    axis.scatter(np.ones_like(READS_QUAL[READS_QUAL >= THRESHOLD_QUAL]) * THRESHOLD_QUAL,
                 READS_QUAL[READS_QUAL >= THRESHOLD_QUAL], marker='o',
                 alpha=0.02, color='grey')
    return axis


def plot_inset(axis: plt.Axes, READS_QUAL: np.ndarray, THRESHOLD_QUAL: int):
    """
       Produce a box/scatter plot of data qualilty vs err.
    """
    # From QUAL to ACCURACY - NB: all data clustered around 99%
    #READS_ACC = 1 - _calculate_error(READS_QUAL)
    #THRESHOLD_ACC = 1 - _calculate_error(THRESHOLD_QUAL)

    # Bar plot
    axis.bar(THRESHOLD_QUAL, len(READS_QUAL[READS_QUAL >= THRESHOLD_QUAL]),
             width=0.8, color='black')
    return axis


def reads_mean_q(FASTQ: str) -> list:
    """
        For each read, compute the MEDIAN score following Phred-scale.
        MEAN is pulled by outliers.
    """
    READS = _retrieve_scores(FASTQ)

    reads_q = list()
    for read in READS:
        mean_err = np.mean([_calculate_error(nt) for nt in read])
        mean_q = -10 * np.log10(mean_err)
        # mean_q = np.mean(read)
        reads_q.append(mean_q)
    return np.array(reads_q)


def main(READS_DIR: str, FIGURE_DIR: str):
    """
        Produce a box/dot plot of number of reads with a Q-score of 10, 20,
        30, ... to help define the chosen threshold.

    """
    if READS_DIR[-1] != str('/'):
        READS_DIR += str('/')

    if FIGURE_DIR[-1] != str('/'):
        FIGURE_DIR += str('/')

    FASTQ_FILES = [f for f in os.listdir(READS_DIR) if f[0] != str('.') and
                   str('.') in f]

    for FILE in FASTQ_FILES:
        fig, ax = plt.subplots(nrows=1, ncols=1,
                               figsize=(5, 4.5))
        ax_inset = fig.add_axes([0.625, 0.135, 0.25, 0.25])
        ax.set_ylabel('Phred score per read')
        ax.set_xlabel('Phred score')

        ax_inset.set_xlabel('Phred score', fontsize=6)
        ax_inset.set_ylabel('Number of reads', fontsize=6)
        ax_inset.xaxis.set_label_position('top')
        ax_inset.tick_params(axis='both', labelsize=5, top=True, labeltop=True,
                             bottom=False, labelbottom=False)

        READS_QUAL = reads_mean_q(READS_DIR + FILE)
        for QUAL in list([10, 15, 20, 25, 30, 35, 40]):
            ax = plot_data(ax, READS_QUAL, QUAL)
            ax_inset = plot_inset(ax_inset, READS_QUAL, QUAL)

        if str('filtered') in READS_DIR:
            plt.savefig(FIGURE_DIR + str('filtered/') +
                        FILE.replace('.fastq.gz',
                        '_QC_per_read_filtered.pdf'), format='pdf',
                        bbox_inches='tight')
        else:
            plt.savefig(FIGURE_DIR + FILE.replace('.fastq.gz',
                        '_QC_per_read.pdf'), format='pdf',
                        bbox_inches='tight')
        plt.close('all')
    return


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
