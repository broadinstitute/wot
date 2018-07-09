#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import wot.io


def main(argv):
    parser = argparse.ArgumentParser(
        description='Convert matrix data formats')
    parser.add_argument('--format', help=wot.commands.FORMAT_HELP, default='loom', choices=wot.commands.FORMAT_CHOICES)
    parser.add_argument('matrix', help='File(s) to convert', nargs='+')
    args = parser.parse_args(argv)
    files = args.matrix

    for f in files:
        name = wot.io.get_filename_and_extension(f)[0]
        wot.io.write_dataset(wot.io.read_dataset(f), name, output_format=args.format, txt_full=True)
