#!/usr/bin/env python3

#*****************************************************************************
#  Name: SVJedi-Tag
#  Description: Genotyping of SVs with linked-reads data
#  Copyright (C) 2025 INRIA
#  Author: Anne Guichard, Mélody Temperville
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as
#  published by the Free Software Foundation, either version 3 of the
#  License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************
import argparse
import sys
from . import run, predict_genotype, lr_stats, likelihood_calibration, construct_graph

def build_parser():
    parser = argparse.ArgumentParser(prog="svjedi-tag")
    subparsers = parser.add_subparsers(required=True)
    run.add_subparser(subparsers)
    construct_graph.add_subparser(subparsers)
    likelihood_calibration.add_subparser(subparsers)
    predict_genotype.add_subparser(subparsers)
    lr_stats.add_subparser(subparsers)
    return parser

def main():
    parser = build_parser()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()