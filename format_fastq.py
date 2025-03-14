#!/usr/bin/env python3

#*****************************************************************************
#  Name: SVJedi-GLR
#  Description: Genotyping of SVs with linked-reads data
#  Copyright (C) 2024 INRIA
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

"""
Module 'svjedi-glr.py': Format the linked-reads FASTQ file to keep barcode information.
"""

import sys
import argparse


#################
# Main function.
#################

def main(args):
    """
    Main method
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-q", 
        "--reads", 
        metavar="<queryReads>", 
        type=str, 
        required=True)

    parser.add_argument(
        "-o", 
        "--outputDir", 
        metavar="<outputDirectory>", 
        type=str,
        required=True)
    args = parser.parse_args()

    file = open(args.outputDir, "w")

    with open(args.reads) as sequenceFile:
        for line in sequenceFile:
            file.write(line.replace(" BX:Z:", "BX:Z:"))

if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")

    else:
        main(sys.argv[1:])