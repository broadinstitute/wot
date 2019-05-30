#!/usr/bin/env python
# -*- coding: utf-8 -*-

from wot.commands import *


def main():
    command_list = [convert_matrix, cells_by_gene_set, census, diff_exp, fates,
                    gene_set_scores, neighborhood_graph, optimal_transport,
                    optimal_transport_validation, trajectory, trajectory_divergence,
                    trajectory_trends, transition_table]
    tool_parser = argparse.ArgumentParser(description='Run a wot command')
    command_list_strings = list(map(lambda x: x.__name__[len('wot.commands.'):], command_list))
    tool_parser.add_argument('command', help='The wot command', choices=command_list_strings)
    tool_parser.add_argument('command_args', help='The command arguments', nargs=argparse.REMAINDER)
    wot_args = tool_parser.parse_args()
    command_name = wot_args.command
    command_args = wot_args.command_args
    cmd = command_list[command_list_strings.index(command_name)]
    sys.argv[0] = cmd.__file__
    parser = cmd.create_parser()
    args = parser.parse_args(command_args)
    cmd.main(args)


if __name__ == '__main__':
    main()
