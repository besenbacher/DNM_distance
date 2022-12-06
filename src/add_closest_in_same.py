"""Module that contains the command line application."""

import argparse
import sys
#from DNM_distance.distance_functions import do_analysis
from distance_functions import do_analysis

#from DNM_distance.readers import *

def main(args = None):
    """
    Run the main program.

    Arguments:
        args: Arguments passed from the command line.

    Returns:
        An exit code.
    """
    parser = argparse.ArgumentParser(
        prog='DNM_distance',
        description='''
        Analyse distribution of distances between de novo mutations
    ''')

    #parser.add_argument('-v', '--verbose', action='store_true')
    #parser.add_argument('-l', '--log_file', type=argparse.FileType('w'))

    parser.add_argument('-f', '--infile', type=argparse.FileType('r'), default = '-',
        help='A file with de novo mutations. First four columns should be: Chrom, pos, ref, alt, pn'
        'Other columns are ignored. Should be sorted by chrom and pos (sort -k1,1 -k2,2n)')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'), default = '-')

    args = parser.parse_args(args)

    dists = []
    last_chrom = 'NA'
    last_pos = 'NA'
    last_line = ''
    D = {}
    cutoff = 1000000000000
    PN_INDEX = 4

    # if args.version:
    #     from DNM_distance import __version__
    #     print("version:", __version__)
    #     print()
    #     return 0

    for line in args.infile:
        L = line.split()
        if L[PN_INDEX] not in D:
            D[L[PN_INDEX]] = []
        D[L[PN_INDEX]].append([L[0], int(L[1]),line.strip(),cutoff,cutoff])

    for pn in D:
        lines = D[pn]
        for i in range(0, len(lines)-1):
            if lines[i][0] == lines[i+1][0]:
                d = lines[i+1][1]-lines[i][1]
                lines[i+1][3] = d
                lines[i][4] = d
        for L in lines:
            if L[3]<L[4]:
                print(L[2], L[3], file=args.outfile)
            else:
                if L[4] == cutoff:
                    print(L[2], "NA", file=args.outfile)
                else:
                    print(L[2], L[4], file=args.outfile)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

