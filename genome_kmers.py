#!/usr/bin/env python
# coding: utf-8
import argparse
import os
import pickle
import time
import numpy as np
import pandas as pd
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
from fasta_parser import parse_fasta, str_punctuation_strip
from alphabet import alphabet_dna
plt.style.use('classic')


start_time = time.process_time()


def get_strand_complement(sequence):
    """Returns the complement strand of the genome."""
    seq = sequence.upper()
    change = str.maketrans('ACGT', 'TGCA')
    return seq.translate(change)


def get_reverse_complement(sequence):
    """Returns the reverse complement strand of the genome."""
    seq = sequence.upper()
    return get_strand_complement(seq)[::-1]

def base_counts(sequence):
    return Counter(base for base in sequence.upper())


def kmer_count(sequence, alphabet, k=1):
    """Returns a dictionary with kmers and it counts."""
    seq = sequence.upper()
    seq_len = len(seq)
    kmers = [seq[i:i+k] for i in range(0, seq_len - k + 1)]
    filterd_kmers = [kmer for kmer in kmers if
                     all(base in set(alphabet) for base in kmer)]
    return Counter(filterd_kmers)


def kmerPositions(sequence, alphabet, k):
    """ returns the position of all k-mers in sequence as a dictionary"""
    kmerPosition = {}
    for i in range(1,len(sequence)-k+1):
        kmer = sequence[i:i+k]
        if all(base in set(alphabet) for base in kmer):
            kmerPosition[kmer] = kmerPosition.get(kmer,[])+[i]
    # combine kmers with their reverse complements
    pairPosition = {}
    for kmer, posList in kmerPosition.items():
        krev = get_reverse_complement(kmer)
        if (kmer < krev):
            pairPosition[kmer] = sorted(posList + kmerPosition.get(krev, []))
        elif (krev < kmer):
            pairPosition[krev] = sorted(kmerPosition.get(krev, []) + posList)
        else:
            pairPosition[kmer] = posList
    return pairPosition


def check_is_palindrome(sequence, kmer):
    """Returns True if the sequence is palindromic other wise False."""
    seq = sequence.upper()
    kmer = kmer.upper()
    return seq.find(kmer[::-1]) == 0


def get_palindromes(sequence, alphabet, k):
    """Returns the count of all the palindromic
    substrings of a genome."""
    mers, _ = kmerPositions(sequence, k)
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


def kmers_clumps(sequence, alphabet, k, window, times):
    """ Find clumps of repeated k-mers in string. A clump occurs when t or
    more k-mers appear within a window of size w. A list of (k-mer, position, count)
    tuples is returned.
    clumpList = kmers_clumps(seq, 9, 500, 6)
    print(len(clumpList))
    print([clumpList[i] for i in range(min(20,len(clumpList)))])
    """
    seq = sequence.upper()
    clumps = []
    # get all kmers and it positions from a string prouced from a alphabet
    kmers = kmerPositions(seq, alphabet, k)
    for kmer, pos in kmers.items():
        for start in range(0, len(pos) - times):
            end = (start + times - 1)
            while (pos[end] - pos[start]) <= window - k:
                end += 1
                if end >= len(pos):
                    break
            if end - start >= times:
                clumps.append((kmer, pos[start], end - start))
    return clumps


# Lets get the positions of all k-mers again
kmers = kmerPositions(sequence, alphabet, k)

# find kmers appearing in the most clumps
clumps = defaultdict(int)
for kmer, start, size in clump_list:
    clumps[kmer] = clumps.get(kmer,0) + 1
top_ten = [k for k in sorted(clumps, reverse=True, key=clumps.get)][:10]

plt.figure(num=None, figsize=(16, 6), dpi=100, facecolor='w', edgecolor='k')
for n, kmer in enumerate(top_ten):
    positions = kmers[kmer]
    plt.text(len(sequence), n+0.4, kmer, fontsize=8)
    plt.plot(positions, [n + 0.5 for i in range(len(positions))], 'o', markersize=4.0)
limit = plt.xlim((start, len(sequence)))


def count_pattern(sequence, pattern):
    counter = 0
    pos = []
    for i in range(len(sequence) - len(pattern) + 1):
        if sequence[i:i+len(pattern)] == pattern:
            counter += 1
            pos.append(i)
    return counter, pos


def get_files(dir_name):
    infiles = []
    for path, subdirs, files in os.walk(dir_name):
        for name in files:
            input_files = os.path.join(path, name)
            if input_files.endswith('fa.gz') or input_files.endswith('.fa')\
                    or input_files.endswith('.fasta') or input_files.endswith('.fa.gz') \
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


def run_job(args):
    dir_out = args.dir_out
    outfile = args.output
    cnt_files = 0
    for filename in filenames:
        print('Starting reading the fasta files', filename)
    for name, sequence in parse_fasta_file(filename):
        name = '/'.join(str_punctuation_strip(name)[1:4:2])
        seq = sequence
        if args.count and args.alphabet:
            kmer_cnts = kmer_count(seq, args.alphabet, args.k)
            pickle.dump( kmer_cnts, open(f"{dir_out}/kmer_cnts.pickle", "wb"))
            df = pd.DataFrame(list(kmer_cnts.items()), columns=['kmers', 'counts'])
        elif args.frequency and args.alphabet:
            kmers_freq = kmers_frequency(kmer_count(seq, args.alphabet, args.k))
            pickle.dump(kmers_freq, open(f"{dir_out}/kmers_freq.pickle", "wb"))
            df = pd.DataFrame(list(kmers_freq.items()), columns=['kmers', 'freq'])
        elif args.pal:
            kmer_palindromes = get_palindromes(seq, args.k)
            pickle.dump(kmer_palindromes, open(f"{dir_out}kmers_palindromes.pickle", "wb"))
            df = pd.DataFrame(list(kmer_palindromes.items()), columns=['kmers', 'counts'])
        elif args.window and args.step:
            cnt_mer_window = get_kmer_count_slide_window(sequence, alphabet, window, step, k)
            pickle.dump(cnt_mer_window, open(f"{dir_out}/kmers_cnt_slide_window.pickle", "wb"))
            df = pd.DataFrame(d).T.fillna(0.0)
            clumps = kmers_clumps(sequence, k, window, times)




end = time.process_time() - start_time
print(f'This script take {end} seconds to finish')
print('Done')
