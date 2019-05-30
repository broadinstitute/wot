#!/usr/bin/env python
# -*- coding: utf-8 -*-

import wot


def create_parser():
    return wot.commands.get_trajectory_or_fates_parser(False)


def main(args):
    wot.commands.run_trajectory_or_fates(args, False)
