import paradys
import argparse


def _get_parser():
    parser = argparse.ArgumentParser(description="PARADYS tool for predicting genes driving dysregulation in cancer.")
    parser.add_argument('--patients', nargs='+', default=[], help='Specify list of to-analyze patient IDs separated by space. Can be omitted when "--all" flag is set.')
    parser.add_argument('--networks', type=str, help='Specify path to dysregulation networks file')
    parser.add_argument('--mutations', type=str, help='Specify path to mutation matrix file')
    parser.add_argument('--outdir', type=str, help='Specify path for output directory')
    parser.add_argument('--kappa', type=int, help='Specify parameter kappa')
    parser.add_argument('--d', type=float, help='Specify dumping factor for PageRank',default=0.85)
    parser.add_argument('--scores', action='store_true', help='Specify whether to compute and return impact scores. Default is false.')
    parser.add_argument('--undirected', action='store_true', help='Specify if line graph for impact scores calculationg should be undirected. Default is false.')
    parser.add_argument('--all', action='store_true', help='Specify if all patients should be analyzed. In this case, --patients can be omitted.')
    return parser


if __name__ == '__main__':
    args = _get_parser().parse_args()
    paradys.process_patients(args.patients, args.kappa, args.d, args.scores, args.mutations, args.networks, args.undirected, args.outdir, args.all)



