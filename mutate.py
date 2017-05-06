from random import random, choice
import sys
from itertools import groupby

#python mutate.py some.fasta 0.01 > mutated.fasta

#create a function that given a sequence, will create n mutations with rate p
alpha="GAUCRYWSMKHBVDN"

def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

def main(fasta_name, mutation_freq):

    for header, seq in fasta_iter(fasta_name):
        seq = list(seq)
        for i, s in enumerate(seq):
            val = random()
            if val < mutation_freq:
                # choose a random nucleotide that's different.
                seq[i] = choice([x for x in "ACTG" if x != s.upper()])
        print ">%s\n%s" % (header, "".join(seq))

def mutate(Mset, mutation_freq):
    new_set = set()
    for seq in Mset:
	new_set.add(seq)
	seq2 = list(seq)
	for i, s in enumerate(seq2):
	    val = random()
	    if val < mutation_freq:
		seq2[i] = choice([x for x in "ACUG" if x != s.upper()])
	new_set.add("".join(seq2))
    return new_set

def RandMotif(n, L):
    #returns a set of n random motifs of length L
    new_set = set()
    for ii in range(n):
	seq =[]
	for kk in range(L):
	    seq.append(choice(alpha))
	new_set.add("".join(seq))
    return new_set


if __name__ == "__main__":
    main(sys.argv[1], float(sys.argv[2]))
