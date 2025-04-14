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

"""
Module 'svjedi-glr.py': Format the linked-reads FASTQ file to keep barcode information.
"""

import sys
import argparse
from xopen import xopen
import re


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

    file = xopen(args.outputDir, "wt")

    with xopen(args.reads, 'rt') as sequenceFile:
        for line in sequenceFile:
            if line.startswith("@") :
                header = line.split('\t')[0]
                res = re.search(r'BX:Z:(\S+)', line)
                barcode = res.group(1)

                header = "".join(header.split())
                barcode = "".join(barcode.split())
                file.write(f'{header+barcode}\n')
            else : 

                file.write(line)
                

if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")

    else:
        main(sys.argv[1:])