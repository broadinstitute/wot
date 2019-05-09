#!/usr/bin/env python
# -*- coding: utf-8 -*-

from wot.commands import *


def main():
    command_list = [convert_matrix, cells_by_gene_set, census, fates,
                    gene_set_scores, diff_exp, neighborhood_graph, optimal_transport,
                    optimal_transport_validation, trajectory, trajectory_divergence,
                    trajectory_trends, transition_table]
    parser = argparse.ArgumentParser(description='Run a wot command')
    command_list_strings = list(map(lambda x: x.__name__[len('wot.commands.'):], command_list))
    parser.add_argument('command', help='The wot command', choices=command_list_strings)
    parser.add_argument('command_args', help='The command arguments', nargs=argparse.REMAINDER)
    wot_args = parser.parse_args()
    command_name = wot_args.command
    command_args = wot_args.command_args
    cmd = command_list[command_list_strings.index(command_name)]
    sys.argv[0] = cmd.__file__
    cmd.main(command_args)


if __name__ == '__main__':
    main()
