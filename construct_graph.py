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
Module 'construct_graph_classes.py': Create the GFA graph file.
"""

from __future__ import print_function
import argparse
import os
import sys
import pickle
from collections import OrderedDict
from Bio import SeqIO
from classes_creationGFA import Chrom, SV


#pylint: disable=line-too-long, disable=trailing-whitespace, disable=too-many-function-args


#################
# Main function.
#################

def main(args):
    """
    Main method
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-v",
        "--vcf",
        metavar="<inputVCF>",
        type=str,
        nargs=1,
        required=True)

    parser.add_argument(
        "-r",
        "--ref",
        metavar="<referenceGenome>",
        type=str,
        nargs=1,
        required=True)

    parser.add_argument(
        "-o",
        "--output",
        metavar="<outputGFAFile>",
        type=str)

    args = parser.parse_args()

    inputVCF = args.vcf[0]
    reference_fasta = args.ref[0]

    if args.output:
        graph_file_name = args.output
    else:
        graph_file_name = inputVCF.split("/")[-1].replace(".vcf", "_graph.gfa")


    ###################################################################
    #A. For each chrom of the reference genome, create a Chrom object.
    ###################################################################
    chromDict = {}
    with open(reference_fasta, "r", encoding='UTF-8') as refFile:
        for record in SeqIO.parse(refFile, "fasta"):
            chr_id, chr_seq = str(record.id), str(record.seq)
            chromDict[chr_id] = Chrom(chr_id, chr_seq)

    #####################################################
    #B. For each SV of the VCF file, create a SV object.
    #####################################################
    with open(inputVCF, "r", encoding='UTF-8') as vcfFile:
        l_discarded = []
        dict_ins_seq = {}
        sv_id = 0

        for line in vcfFile:
            if line.startswith('#'):                #headers lines
                continue

            else:
                if len(line.rstrip().split("\t")) >= 8 :
                    chromID, pos, __, __, alt, __, __, info, *__ = line.rstrip().split("\t")
                else : 
                    chromID, pos, __, __, alt, __, __, info = line.rstrip().split("\t")
                
                sv_id += 1

                pos = int(pos)
                sv_type = get_info(info, "SVTYPE")
                end = int(get_info(info, "END"))
                sv_coords = [pos, end]              #'pos' is 1-based and incl. ; 'end' is 1-based and excl.


                # DELetions.
                if sv_type == "DEL":
                    sv_format = format_nonBND_id(sv_type, pos, end)
                    sv_length = end - pos

                # INSertions.
                elif sv_type == "INS":
                    end = pos
                    sv_format = format_nonBND_id(sv_type, pos, end)
                    sv_length = int(get_info(info, "SVLEN"))

                    # INS seq not in ALT field, get seq from INFO field.
                    if alt.startswith("<"):
                        if any(["LEFT_SVINSSEQ=" in info, "RIGHT_SVINSSEQ=" in info]):
                            # left_insseq = get_info(info, "LEFT_SVINSSEQ")
                            # right_insseq = get_info(info, "RIGHT_SVINSSEQ")
                            l_discarded.append(line.rstrip())
                            continue

                        elif "SEQ=" in info:
                            dict_ins_seq[sv_id] = get_info(info, "SEQ")

                        else:
                            l_discarded.append(line.rstrip())
                            continue

                    # INS seq in ALT field but not in dict_sv_alt (no correction of pos required).
                    elif sv_id not in dict_ins_seq:
                        dict_ins_seq[sv_id] = alt.upper()

                # INVersions.
                elif sv_type == "INV":
                    sv_format = format_nonBND_id(sv_type, pos, end)         #sv_format: "INV-<pos>-<end>"
                    sv_length = end - pos

                # BND.
                elif sv_type == "BND":
                    sv_format = format_BND_id(pos, alt)
                    sv_length = end - pos                                   #TODO: determine sv_length

                # DUPlications and other SV types.
                else:
                    continue

                # Create a SV object.
                sv = SV(sv_id, sv_type, chromID, sv_coords, sv_format, sv_length)

                if sv.type != "BND":
                    # Add the SV breakpoints coords (pos, end) associated to their sv (objects) to the 'bkpts' dictionary of the corresponding Chrom object.
                    chromDict[chromID].addBreakpoint(pos, sv)               #bkpts are 1-based, with the first incl.    #TODO: to check if still working with entre-croisés SVs
                    chromDict[chromID].addBreakpoint(end, sv)

                    # Add the SV to the 'svs' list of the corresponding Chrom object.
                    chromDict[chromID].addSV(sv)

    ###############################
    #C. Create the GFA graph file.
    ###############################
    headerLines = []
    SLines = []
    PLines = []
    LLines = []
    leftEdges = {}
    rightEdges = {}

    with open(graph_file_name, "w", encoding='UTF-8') as graphFile:
        for chrID, chrObject in chromDict.items():

            #Write InfoSV line (containing the chrom ID and the list of SVs for the current Chrom object). 
            #chrObject.writeInfoSVLines(headerLines)

            #0. Process the case where there is no breakpoint (e.g. no SV).
            ###############################################################
            if len(chrObject.bkpts) == 0:
                
                # Create one GFA node.
                node_start = 0                              #0-based incl.
                node_end = int(chrObject.length)            #0-based excl.
                gfaNode_id = ":".join([str(chrID), str(node_start+1), str(node_end)])      #positions are 1-based to match with the VCF file, and incl.

                # Add the GFA node to the 'gfaNodes' dictionary of the current Chrom object and to the graph.
                chrObject.addGFANode(gfaNode_id)
                chrObject.writeSLines(SLines)

            # Process the case where there is at least one SV.
            else:

                # Sort the 'bkpts' dictionary of the corresponding Chrom object by 'bkpt_coord' (key).
                chrBkpts = OrderedDict(sorted(chrObject.bkpts.items()))
                chrBkpts_keys = list(chrBkpts.keys())

                #1. Create the GFA nodes and update sv.gfaNodes using the 'bkpts' dictionary of the current Chrom object.
                #########################################################################################################
                ## (excl. INS nodes and alternative links).
                for i in range(-1, len(chrBkpts_keys)):

                    # First node of region/graph component.
                    if i == -1:
                        node_start = 0                                          #0-based incl.
                        node_end = int(chrBkpts_keys[0]) - 1           #0-based excl.
                        svs = chrBkpts[int(node_end)+1]

                    # Last node of region/graph component.
                    elif i == len(chrBkpts_keys)-1:
                        node_start = int(chrBkpts_keys[i]) - 1         #0-based incl.
                        node_end = int(chrObject.length)                        #0-based excl.
                        svs = chrBkpts[int(node_start)+1]

                    # Middle node of region/graph component.
                    else:
                        node_start = int(chrBkpts_keys[i]) - 1         #0-based incl.
                        node_end = int(chrBkpts_keys[i+1]) - 1         #0-based excl.
                        svs = chrBkpts[int(node_start)+1] + chrBkpts[int(node_end)+1]

                    gfaNode_id = ":".join([str(chrID), str(node_start+1), str(node_end)])        #positions are 1-based to match with the VCF file, and incl.
                    associate_GFANode_To_SVsBkptsAdj(svs, node_start, node_end, gfaNode_id)      #TODO: to check for INS nodes

                    # Add the current GFA node to the 'gfaNodes' dictionary of the current Chrom object.
                    chrObject.addGFANode(gfaNode_id)

                # Add the GFA nodes to the graph and update the 'leftEdges' and 'rightEdges' dictionary when creating the L lines.
                chrObject.writeSLines(SLines)
                chrObject.writeRefPLine(PLines)
                chrObject.writeRefLLines(LLines, leftEdges, rightEdges)

            # Creation sv edges
            for sv in chrObject.svs:
                chrObject.writeAltLLines(LLines, sv, leftEdges, rightEdges, SLines, dict_ins_seq)
                #TODO: for BND



     #4. Write the lines to the GFA graph file.
        ##########################################
        for line in headerLines:
            graphFile.write(line)
        for line in SLines:
            graphFile.write(line)
        for line in PLines:
            graphFile.write(line)
        for line in LLines:
            graphFile.write(line)

    # Dump the dictionary 'ChromDict' in a pickle file.
    pickleFile = str(graph_file_name).rsplit(".gfa", maxsplit=1)[0] + "_chromDict.pickle"
    with open(pickleFile, "wb") as f1:
    #    pickle.dump(gfaNode2svRegionsDict, f1)
        pickle.dump(chromDict, f1)




#############
# Functions.
#############

def get_info(info, label):
    """Method to return the value corresponding to the requested label in the '#INFO' field of the VCF file."""
    #Info label first in info field
    if info.split(";")[0].startswith(label+"="):
        return info.split(label+"=")[1].split(";")[0]

    #Info label last in info field
    elif info.split(";")[-1].startswith(label+"="):
        return info.split(label.join([";", "="]))[1]

    else:
        return info.split(label.join([";", "="]))[1].split(";")[0]


def format_nonBND_id(sv_type, pos, end):
    """Method to format the deletions/insertions/inversions (DEL/INS/INV) SV: `<DEL|INS|INV>-<pos>-<end>`."""
    return '-'.join([sv_type, str(pos), str(end)])


def format_BND_id(pos, alt):
    """Method to format the BND SV: `BND-<alt>`."""
    if "[" in alt:
        parsed_alt = list(filter(bool, alt.split("[")))
        if ":" in parsed_alt[1]:
            s = parsed_alt[0]
        else:
            s = parsed_alt[1]

    elif "]" in alt:
        parsed_alt = list(filter(bool, alt.split("]")))
        if ":" in parsed_alt[1]:
            s = parsed_alt[0]
        else:
            s = parsed_alt[1]

    else:
        return "BND-format"

    if len(s) > 1:
        return "BND-format"

    alt = alt.replace(s, pos)
    return '-'.join(["BND", alt])


def associate_GFANode_To_SVsBkptsAdj(list_of_SV_objects, node_start, node_end, gfaNode_id):
    """Method to associate a GFA node to SVs breakpoints adjacencies."""
    for sv in list_of_SV_objects:

        #b1.left
        if sv.coords[0] == node_end + 1:            #sv.coords[0] == sv.pos
            sv.gfaNodes.append(gfaNode_id)

        # b1.right
        if (sv.coords[0] == node_start + 1) and (len(sv.gfaNodes) == 1):
            sv.gfaNodes.append(gfaNode_id)

        # b2.left
        if (sv.coords[1] == node_end+1) and (len(sv.gfaNodes) == 2):              #sv.coords[1] == sv.end
            sv.gfaNodes.append(gfaNode_id)

        # b2.right
        if sv.coords[1] == node_start+1:
            sv.gfaNodes.append(gfaNode_id)

##############################################
if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")

    else:
        main(sys.argv[1:])
