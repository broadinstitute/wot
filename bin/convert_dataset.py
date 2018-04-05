#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot.io

parser = argparse.ArgumentParser(
    description='Convert data formats')
parser.add_argument('--format', help='Output file format', default='loom')
parser.add_argument('file', help='File(s) to convert', nargs='+')
args = parser.parse_args()
files = args.file

for f in files:
    name = wot.io.get_file_basename_and_extension(f)[0]
    wot.io.write_dataset(wot.io.read_dataset(f), name, output_format=args.format, txt_full=True)
