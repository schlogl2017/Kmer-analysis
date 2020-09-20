#!usr/bin/env python
# Read fasta files in a compressed gzip form or not
import itertools
import gzip


def is_header(line):
    """Check if the line starts with '>'."""
    return line[0] == '>'


def parse_fasta_file(filename):
    """It reads and returns a name and the sequence from a file (compressed or not)."""
    if filename.endswith('.gz'):
        opener = lambda filename: gzip.open(filename, 'rt')
    else:
        opener = lambda filename: open(filename, 'r')

    with opener(filename) as f:
        fasta_iter = (it[1] for it in itertools.groupby(f, is_header))
        for name in fasta_iter:
            name = name.__next__()[1:].strip()
            sequences = ''.join(seq.strip() for seq in fasta_iter.__next__())
            yield name, sequences


def str_punctuation_strip(word):
    punctuation = '!"#$%&\'()*+,-/:;<=>?@[\\]^`{|}~'
    for _ in word:
        for p in punctuation:
            word = word.replace(p, ' ')
    return word.strip().split()


def count_fasta_files(filename):
    if filename.endswith('.gz'):
        opener = lambda filename: gzip.open(filename, 'rt')
    else:
        opener = lambda filename: open(filename, 'r')

    with opener(filename) as f:
        return sum(g for g, _ in itertools.groupby(f, key=is_header))
