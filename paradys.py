import paradys
import argparse


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--patients', nargs='+', default=[], help='Specify list of patient IDs separated by space')
    parser.add_argument('--networks', type=str, help='Specify path to dysregulation networks file')
    parser.add_argument('--mutations', type=str, help='Specify path to mutation matrix file')
    parser.add_argument('--outdir', type=str, help='Specify path for output directory')
    parser.add_argument('--kappa', type=int, help='Specify parameter kappa')
    parser.add_argument('--d', type=float, help='Specify dumping factor for PageRank',default=0.85)
    parser.add_argument('--scores', action='store_true', help='Specify whether to return scores')
    parser.add_argument('--directed', action='store_true', help='Specify if networks are directed')
    return parser


if __name__ == '__main__':
    args = _get_parser().parse_args()
    paradys.process_patients(args.patients, args.kappa, args.d, args.scores, args.mutations, args.networks, args.directed, args.outdir)



