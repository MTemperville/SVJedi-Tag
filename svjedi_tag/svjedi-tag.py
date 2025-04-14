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
Module 'svjedi-tag.py': Pipeline of SVJedi-Tag.
"""

import sys
import argparse
import os
import re
import subprocess


#pylint: disable=line-too-long, disable=trailing-whitespace, disable=consider-using-f-string


def main(args):

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-v", 
        "--vcf", 
        metavar="<inputVCF>", 
        type=str,
        required=True)

    parser.add_argument(
        "-r", 
        "--ref", 
        metavar="<referenceGenome>", 
        type=str, 
        required=True)

    parser.add_argument(
        "-q", 
        "--reads", 
        metavar="<queryReads>", 
        type=str, 
        required=True)

    parser.add_argument(
        "-p", 
        "--prefix", 
        metavar="<outFilesPrefix>", 
        type=str, 
        required=True)

    parser.add_argument(
        "-t", 
        "--threads", 
        metavar="<threadNumber>", 
        type=int, 
        default=1)

    parser.add_argument(
        "-s",
        "--regionSize",
        metavar="<regionSize (default 10000)>",
        type=int,
        required=False,
        default=10000)

    parser.add_argument(
        "-a", 
        "--gaf", 
        metavar="<alignmentGAFFile>", 
        type=str,
        required=False)
    
    parser.add_argument(
        "-g", 
        "--gfa", 
        metavar="<Graphe File GFA>", 
        type=str,
        required=False)

    args = parser.parse_args()
    inVCF = args.vcf
    inREF = args.ref
    inFQ = args.reads
    outPrefix = args.prefix
    threads = args.threads
    regionSize = args.regionSize

    script_path = os.path.abspath(__file__)
    script_dir = os.path.dirname(script_path)

    if args.gaf or args.gfa :
        if not args.gaf : 
            print('### Gaf file (alignement) necessary please use parameter --gaf ###')
            
        if not args.gfa : 
            print('### Gfa file (graphe) necessary please use parameter --gfa ###')
            
        else:
            print("### Mapping of linked-reads onto graph already done ###")
            outGAF = os.path.abspath(args.gaf)
            outGFA = os.path.abspath(args.gfa)
            print("Alignment GAF file: " + str(outGAF))

            #### Analyze barcode signal & Genotype.
            print("### Analyze barcode signal & Genotype ###")
            outVCF = outPrefix + "_genotype.vcf"
            c6 = "python3 {}/predict_genotype.py -a {} -v {} -o {} -s {} -g {}".format(script_dir, outGAF, inVCF, outVCF,regionSize, outGFA)
            subprocess.run(c6, shell=True, check=True)
        
    else:
        #### Create variant graph.
        print("### Create variant graph ###")
        outGFA = outPrefix + ".gfa"
        c1 = "python3 {}/construct_graph.py -v {} -r {} -o {}".format(script_dir,inVCF, inREF, outGFA)
        subprocess.run(c1, shell=True, check=True)

        ### Index graph.
        print("### Index graph ###")
        c3 = "vg autoindex --workflow giraffe -g {} -p {}".format(outGFA, outPrefix)
        subprocess.run(c3, shell=True, check=True)

        ### Map linked-reads on graph.
        print ("### Map linked-reads on graph ###")
        outGBZ = outPrefix + ".giraffe.gbz"
        outMIN = outPrefix + ".min"
        outDIST = outPrefix + ".dist"
        outGAF = outPrefix + "_vgGiraffe.gaf"
        c4 = "vg giraffe -t {} -Z {} -m {} -d {} -f {} -i -o gaf --named-coordinates > {}".format(threads, outGBZ, outMIN, outDIST, inFQ, outGAF)
        subprocess.run(c4, shell=True, check=True)

        #### Analyze barcode signal & Genotype.
        print("### Analyze barcode signal & Genotype ###")
        outVCF = outPrefix + "_genotype.vcf"
        c6 = "python3 {}/predict_genotype.py -a {} -v {} -o {} -s {} -g {}".format(script_dir, outGAF, inVCF, outVCF,regionSize, outGFA)
        subprocess.run(c6, shell=True, check=True)


if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")

    else:
        main(sys.argv[1:])
