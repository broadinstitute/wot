#!/usr/bin/env python
# -*- coding: utf-8 -*-


import wot.io
import wot.ot


def main(argv):
    parser = wot.ot.OptimalTransportHelper.create_base_parser('Compute transport maps between pairs of time points')
    parser.add_argument('--format', help=wot.commands.FORMAT_HELP, default='loom', choices=wot.commands.FORMAT_CHOICES)

    args = parser.parse_args(argv)
    ot_helper = wot.ot.OptimalTransportHelper(args)

    params_writer = None
    # if args.solver is 'floating_epsilon':
    #     params_writer = open(args.out + '_params.txt', 'w')
    #     params_writer.write('t1' + '\t' + 't2' + '\t' + 'epsilon' + '\t' + 'lambda1' + '\t' + 'lambda2' +
    #                         '\n')

    def callback(cb_args):
        result = cb_args['result']
        # if args.solver is 'floating_epsilon':
        #     params_writer.write(
        #         str(cb_args['t0']) + '\t' + str(cb_args['t1']) + '\t' + str(result['epsilon']) + '\t' + str(
        #             result['lambda1']) + '\t' + str(
        #             result['lambda2']) + '\n')

        # save the tranport map

        if args.verbose:
            print('Saving transport map')

        filename = args.out + '_' + str(cb_args['t0']) + '_' + str(cb_args['t1'])

        row_meta = cb_args['df0'].copy()
        row_meta.drop(['cell_growth_rate', 'day'], axis=1, inplace=True)
        row_meta['g'] = cb_args['g']

        col_meta = cb_args['df1'].copy()
        col_meta.drop(['cell_growth_rate', 'day'], axis=1, inplace=True)
        wot.io.write_dataset(wot.Dataset(result['transport'], row_meta, col_meta), filename,
                             output_format=args.format)

    ot_helper.compute_transport_maps(callback)

    if params_writer is not None:
        params_writer.close()
