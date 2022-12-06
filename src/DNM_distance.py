"""Module that contains the command line application."""

import argparse
import py2bit
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

    parser.add_argument('mutations', 
        type=str, 
        help='A file with de novo mutations. First four columns should be: Chrom, pos, ref, alt, pn. '
        'Other columns are ignored. Non-SNV variants are ignored.')
    parser.add_argument('ref_genome', 
        help='Reference genome in 2bit format', type=str)
    parser.add_argument('-n', '--n_random', default=2, type=int, 
        help='How many times should mutation positions be permuted')


    parser.add_argument('--bed', type=str,
        help='bed-file describing regions that should be counted. May be gzipped.')
    
    parser.add_argument('--obs_intpos', type=argparse.FileType('w'), default = None)
    parser.add_argument('--random_intpos', type=argparse.FileType('w'), default = None)
    #parser.add_argument('--obs_intpos', type=argparse.FileType('w'), default = None)
    parser.add_argument('--random_chrompos', type=argparse.FileType('w'), default = None)

    # parser.add_argument('--wig', type=str,
    #     help='wig-file describing regions that should be counted. May be gzipped. '
    #         'The context at a position will be weigthed by the value from the '
    #         'wig-file at that position. The output counts will thus be floats '
    #         'and not integers')
    # parser.add_argument('--all_autosomes', action="store_true",
    #     help='All parts of the autosomes will be counted')

    args = parser.parse_args(args)

    # if args.version:
    #     from DNM_distance import __version__
    #     print("version:", __version__)
    #     print()
    #     return 0

    if 'ref_genome' not in args:
        print('Error: ref_genome (as 2bit file) must be specified')
        print()
        parser.print_help()
        return 0

    tb = py2bit.open(args.ref_genome)

    do_analysis(args.mutations, args.bed, tb, args.n_random, args.obs_intpos, args.random_intpos, args.random_chrompos)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
