#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import wot.io

parser = argparse.ArgumentParser(
    description='Convert data formats')
parser.add_argument('--file', help='File(s) to convert', required=True, action='append')
parser.add_argument('--format', help='Output file format', default='loom')

args = parser.parse_args()
files = args.file

for f in files:
    ds = wot.io.read_dataset(f)
    name = wot.io.get_file_basename_and_extension(f)[0]
    wot.io.write_dataset(ds, name, output_format=args.format)
