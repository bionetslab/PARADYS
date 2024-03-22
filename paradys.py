import paradys
import argparse


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('patients', nargs='+', default=[], help='specify list of patient IDs')
    parser.add_argument('networks', type=str, help='specify path to dysregulation networks')
    parser.add_argument('mutations', type=str, help='specify path to mutation matrix')
    parser.add_argument('kappa', type=int, help='specify parameter kappa')
    parser.add_argument('d', type=float, help='specify dumping factor',default=0.85)
    parser.add_argument('scores', type=bool, help='specify if return scores',default=False)
    parser.add_argument('directed', type=bool, help='specify if networks are directed',default=False)
    return parser


if __name__ == '__main__':
    args = _get_parser().parse_args()
    _, _ = paradys.process_patients(args.patients, args.kappa, args.d, args.scores, args.mutations, args.networks, args.directed)



