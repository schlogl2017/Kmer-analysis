#!/usr/bin/env python
# coding: utf-8

# genome_kmers usage
# to count = python genome_kmers.py -p Data -d Results -c -k int
# frequency = python genome_kmers.py -p Data -d Results -k int -f -a(default ACGT)
# palindromes = python genome_kmers.py -p Data -d Results -pl -k 5


import os
import time
import argparse
import pickle
import pandas as pd
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
from fasta_reader import *
from alphabet import *


start_time = time.process_time()


parser = argparse.ArgumentParser(prog='Kmer analysis',
                                 usage='%(prog)s [options] path',
                                 description='Search for all k mers of length k from a genome.')
parser.add_argument('-p',
                    '--path',
                    metavar='path',
                    type=str,
                    required=True,
                    dest='path',
                    help='Path to the files')
parser.add_argument('-o',
                    '--out',
                    type=str,
                    dest='output',
                    help='name output files.')
parser.add_argument('-d',
                    '--dir_out',
                    required=True,
                    type=str,
                    dest='dir_out',
                    help='name output files.')
parser.add_argument('-k',
                    '--length',
                    type=int,
                    action="store",
                    dest='k',
                    help='length of the kmers')
parser.add_argument('-c',
                    '--count',
                    action="store_true",
                    dest='count',
                    help='Count kmers')
parser.add_argument('-f',
                    '--freq',
                    action="store_true",
                    dest='freq',
                    help='Frequency of the kmers')
parser.add_argument('-pl',
                    '--pal',
                    action="store_true",
                    dest='pal',
                    help='Check for palindromic sequences.')
parser.add_argument('--patter',
                    '-pat',
                    type=str,
                    action="store",
                    dest='pattern',
                    help='Count a pattern in a sequence.')
parser.add_argument('-sld',
                    '--slide',
                    type=str,
                    action="store",
                    dest='slide',
                    help='Count a pattern in a sequence.')
parser.add_argument('-plt',
                    '--plot',
                    type=str,
                    action="store",
                    dest='plot',
                    help='Count a pattern in a sequence.')
parser.add_argument('-a',
                    '--alph',
                    default=iupac_dna,
                    dest='alphabet',
                    help='Alphabet for the sequences.')
parser.add_argument('-w',
                    '--window',
                    type=int,
                    action="store",
                    dest='window',
                    help='Window size for the sequences chunks.')
parser.add_argument('-s',
                    '--step',
                    type=int,
                    action="store",
                    dest='step',
                    help='Step size for the sequences chunks.')
args = parser.parse_args()


def get_strand_complement(sequence):
    """Returns the complement strand of the genome."""
    seq = sequence.upper()
    change = str.maketrans('ACGT', 'TGCA')
    return seq.translate(change)


def get_reverse_complement(sequence):
    """Returns the reverse complement strand of the genome."""
    seq = sequence.upper()
    return get_strand_complement(seq)[::-1]


def count_pattern(sequence, pattern):
    '''Return the count the input pattern found in to a give string.'''
    return Counter([sequence[i:i+len(pattern)] for i in range(len(sequence)-len(pattern) + 1)
                    if sequence[i:i+len(pattern)] == pattern])


def kmer_count(sequence, alphabet, k=1):
    """Returns a dictionary with kmers and it counts."""
    seq = sequence.upper()
    seq_len = len(seq)
    kmers = [seq[i:i+k] for i in range(0, seq_len - k + 1)]
    filterd_kmers = [kmer for kmer in kmers if 
                     all(base in set(alphabet) for base in kmer)]
    return Counter(filterd_kmers)


def kmers_frequency(kmers_dict):
    num_kmers = len(kmers_dict.keys())
    kmer_freq = defaultdict(float)
    for kmer, count in kmers_dict.items():
        kmer_freq[kmer] = (count / num_kmers)
    return kmer_freq


def kmer_positions(sequence, alphabet, k):
    """ returns the position of all k-mers in sequence as a dictionary"""
    mer_position = defaultdict(list)
    for i in range(1,len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if all(base in set(alphabet) for base in kmer):
            mer_position[kmer] = mer_position.get(kmer,[]) + [i]
    # combine kmers with their reverse complements
    pair_position = defaultdict(list)
    for kmer, pos in mer_position.items():
        krev = get_reverse_complement(kmer)
        if kmer < krev:
            pair_position[kmer] = sorted(pos + mer_position.get(krev, []))
        elif krev < kmer:
            pair_position[krev] = sorted(mer_position.get(krev, []) + pos)
        else:
            pair_position[kmer] = pos
    return pair_position


def check_is_palindrome(sequence, kmer):
    """Returns True if the sequence is palindromic other wise False."""
    seq = sequence.upper()
    kmer = kmer.upper()
    return seq.find(kmer[::-1]) == 0


def get_palindromes(sequence, alphabet, k):
    """Returns the count of all the palindromic
    substrings of a genome."""
    mers = kmer_positions(sequence, alphabet, k)
    kmers = list(mers.keys())
    rev_kmers = [get_reverse_complement(kmer) for kmer in kmers]
    palindromes = defaultdict()
    for mer1, mer2 in zip(kmers, rev_kmers):
        if mer1 == mer2:
            if mer1 in mers:
                pos = mers[mer1]
                palindromes[mer1] = palindromes.get(mer1, []) + pos
                #palindromes.add((mers, pos))
    return palindromes


def kmers_clumps(sequence, alphabet, k, w, t):
    """ Find clumps of repeated k-mers in string.
    A clump occurs when t or more k-mers appear
    within a window of size w. A list of (k-mer, position, count)
    tuples is returned
    clumpList = kmers_clumps(seq, 9, 500, 6)
    print(len(clumpList))
    print([clumpList[i] for i in range(min(20,len(clumpList)))])
    """
    seq = sequence.upper()
    clumps = []
    kmers = kmer_positions(sequence, alphabet, k)
    for kmer, pos in kmers.items():
        for start in range(len(pos) - t):
            end = start + t - 1
            while ((pos[end] - pos[start]) <= w - k):
                end += 1
                if (end >= len(pos)):
                    break
            if end - start >= t and kmer not in itertools.chain(*clumps):
                clumps.append((kmer, pos[start], pos, end - start))
    return clumps


def get_files(dir_name):
    infiles = []
    for path, subdirs, files in os.walk(dir_name):
        for name in files:
            input_files = os.path.join(path, name)
            if input_files.endswith('fa.gz') or input_files.endswith('.fa') \
            or input_files.endswith('.fasta') or input_files.endswith('.fa.gz')\
            or input_files.endswith('.fna') or input_files.endswith('.fna.gz'):
                infiles.append(input_files)
    return infiles



def get_chunks(sequence, window, step=1):
    """Returns a chunk of length of window_size and the end of the window size"""
    k = len(sequence)
    for i in range(0, k - window + 1, step):
        end = i + window
        chunk = sequence[i:end]
        assert len(chunk) == window
        yield chunk, i, end


def get_kmer_count_slide_window(sequence, alphabet, window, step, k):
    slide_mer_count = defaultdict(Counter)
    for chunk, s, e in get_chunks(sequence, window, step):
        pos = '_'.join((str(s), str(e)))
        slide_mer_count[pos].update(kmer_count(chunk, alphabet, k))
    return slide_mer_count


dir_name = args.path
filenames = get_files(dir_name)
print(f'The number of files to be read are {len(filenames)}')
print(f'The files to be read are: {filenames}')
dir_out = args.dir_out
outfile = args.output
cnt_files = 0
print(f'Starting reading the fasta files!')
for filename in filenames:
    for name, sequence in parse_fasta_file(filename):
        name = '/'.join(str_punctuation_strip(name)[1:4:2])
        seq = sequence
        if args.count:
            kmer_cnts = kmer_count(seq, args.alphabet, args.k)
            pickle.dump( kmer_cnts, open(f"{dir_out}/kmer_cnts.pkl", "wb"))
            df = pd.DataFrame(list(kmer_cnts.items()), columns=['kmers', 'counts'])
            df.to_pickle(f"{dir_out}/kmer_cnts_data_frame.pkl", compression='gzip')
        elif args.freq:
            kmers_freq = kmers_frequency(kmer_count(seq, args.alphabet, args.k))
            #print(kmers_freq)
            pickle.dump(kmers_freq, open(f"{dir_out}/kmers_freq.pkl", "wb"))
            df = pd.DataFrame(list(kmers_freq.items()), columns=['kmers', 'freq'])
            df.to_pickle(f"{dir_out}/kmer_freqs_data_frame.pkl", compression='gzip')
        elif args.pal:
            kmer_palindromes = get_palindromes(seq, args.alphabet, args.k)
            print(kmer_palindromes)
            pickle.dump(kmer_palindromes, open(f"{dir_out}/kmers_palindromes.pkl", "wb"))
            df = pd.DataFrame(list(kmer_palindromes.items()), columns=['kmers', 'counts'])
            print(df.head(10))
            df.to_pickle(f"{dir_out}/kmers_palindromes_data_frame.pkl", compression='gzip')
        elif args.pattern:
            count = count_pattern(seq, args.pattern)
            print(f'The pattern: {args.pattern} was counted {count} in the sequence {name}')
        elif args.slide:
            cnt_mer_window = get_kmer_count_slide_window(sequence,
                                                         args.alphabet,
                                                         args.window,
                                                         args.step,
                                                         args.k)
            pickle.dump(cnt_mer_window, open(f"{dir_out}/kmers_cnt_slide_window.pkl", "wb"))
            df = pd.DataFrame(cnt_mer_window).T.fillna(0.0)
            df.to_pickle(f"{dir_out}/cnt_mer_window_data_frame.pkl", compression='gzip')
        elif args.plot:
            # save clumps plots
            clump_list = kmers_clumps(seq, args.alphabet, args.k, args.window, args.times)
            # Lets get the positions of all k-mers again
            kmers = kmer_positions(seq, args.alphabet, args.k)
            # find kmers appearing in the most clumps
            clumps = defaultdict(int)
            for kmer, start, pos_lst, size in clump_list:
                clumps[kmer] = clumps.get(kmer,0) + 1
            top_ten = [k for k in sorted(clumps, reverse=True, key=clumps.get)][:10]
            plt.figure(num=None, figsize=(16, 6), dpi=100, facecolor='w', edgecolor='k')
            for n, kmer in enumerate(top_ten):
                positions = kmers[kmer]
                plt.text(len(seq), n+0.4, kmer, fontsize=8)
                plt.plot(positions, [n + 0.5 for i in range(len(positions))], 'o', markersize=4.0)
            limit = plt.xlim((start, len(sequence)))
            plt.savefig(f"{dir_out}/{name}.pdf")
        cnt_files += 1


print(f'Readed files:{cnt_files}')


end = time.process_time() - start_time
print(f'This script take {end} seconds to finish!')
print('Done!')
