#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import wot


def main(argv):
    parser = argparse.ArgumentParser(
        'Generate fates for cell sets generated at the given time.')
    wot.commands.run_trajectory_or_fates(parser, argv, True)
