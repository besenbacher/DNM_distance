import argparse
import sys
from collections import Counter

def main(args = None):
    """
    Run the main program.

    Arguments:
        args: Arguments passed from the command line.

    Returns:
        An exit code.
    """
    parser = argparse.ArgumentParser(
        prog='all_pair_distances',
        description='''
        Analyse distribution of all pairwise distances between de novo mutations.
    ''')

    #parser.add_argument('-v', '--verbose', action='store_true')
    #parser.add_argument('-l', '--log_file', type=argparse.FileType('w'))

    parser.add_argument('mutations', 
        type=argparse.FileType('r'), 
        help='A file with de novo mutations. First four columns should be: Chrom, pos, ref, alt, pn'
            'Other columns are ignored. Non-SNV variants are ignored.')

    args = parser.parse_args(args)
    chrom_poses = {}
    for line in args.mutations:
        chrom, pos, ref, alt, pn = line.split()
        if chrom not in chrom_poses:
            chrom_poses[chrom] = []
        chrom_poses[chrom].append((int(pos),pn))

    C = Counter()
    for chrom in chrom_poses:
        for i in range(len(chrom_poses[chrom])):
            for j in range(i+1,len(chrom_poses[chrom])):
                pos_i, pn_i = chrom_poses[chrom][i]
                pos_j, pn_j = chrom_poses[chrom][j]
                C[(abs(pos_i-pos_j), int(pn_i==pn_j))] += 1
                #print(abs(pos_i-pos_j), int(pn_i==pn_j))

    print("dist insame count")
    for x in C:
        dist, insame = x
        print(dist, insame, C[x])

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))