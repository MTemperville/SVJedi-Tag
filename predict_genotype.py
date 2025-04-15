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
Module 'predict_genotype.py': Process the barcodes signal to genotype the SVs.
"""

from __future__ import print_function
import argparse
import os
import re
import sys
import pickle
import math
import statistics                   #analysis
from decimal import Decimal

from pgGraphs.graph import Graph
from collections import deque
from pgGraphs.abstractions import Orientation

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
        "-a",
        "--gaf",
        metavar="<alignmentGAFFile>",
        type=str,
        nargs=1,
        required=True)

    parser.add_argument(
        "-g",
        "--gfa",
        metavar="<GrapheGFAFile>",
        type=str,
        nargs=1,
        required=True)

    parser.add_argument(
        "-v",
        "--vcf",
        metavar="<inputVCF>",
        type=str,
        nargs=1,
        required=True)
    
    parser.add_argument(
        "-s",
        "--regionSize",
        metavar="<regionSize>",
        type=int,
        required=True)

    parser.add_argument(
        "-o", 
        "--output", 
        metavar="<outputVCF>", 
        type=str,
        nargs=1,
        required=True)

    args = parser.parse_args()

    inputGFA = args.gfa[0]
    inputGAF = args.gaf[0]
    inputVCF = args.vcf[0]
    regionSize = args.regionSize
    outputVCF = args.output[0]


    # Load the dictionary 'chromDict' from pickle file.
    pickleFile = str(inputGAF).rsplit("_vgGiraffe.gaf", maxsplit=1)[0] + "_chromDict.pickle"
    with open(pickleFile, "rb") as pf:
        chromDict = pickle.load(pf)

    gfa_graph:Graph = Graph(
        gfa_file=inputGFA,
        with_sequence=False,
        with_reverse_edges=True,
        low_memory=False,
    )
    gfa_graph.compute_orientations()

    svsDict = {}
    gfaNode2svRegionsDict = {}
    for chr, chrObject in chromDict.items() :
        for sv in chrObject.svs:

            # Regions of interest are created from nodes directly adjacent to the inversion breakpoint.
            # If the node is smaller than the set region size, then it is created using a deep graph traversal.

            ## adjLeft.
            create_region(sv,  sv.gfaNodes[0], Orientation.REVERSE ,regionSize, "adjLeft", gfaNode2svRegionsDict, gfa_graph)              

            ## adjRight.
            create_region(sv,  sv.gfaNodes[-1], Orientation.FORWARD ,regionSize, "adjRight", gfaNode2svRegionsDict, gfa_graph)
            
            ## nodeSVbegin.
            create_region(sv,  sv.gfaNodes[1], Orientation.FORWARD ,regionSize, "nodeSVbegin", gfaNode2svRegionsDict,gfa_graph)

            ## nodeSVend.
            create_region(sv,  sv.gfaNodes[-2], Orientation.REVERSE ,regionSize, "nodeSVend", gfaNode2svRegionsDict, gfa_graph)

            svsDict[sv.id] = sv

    # Create a file containing the analysis results.
    analysis_file = str(inputGAF).rsplit("_vgGiraffe.gaf", maxsplit=1)[0] + "_analysis.txt"
    analysisFile = open(analysis_file, "w", encoding='UTF-8')
    analysisFile.write("\t".join(["SV", "Method", "Genotype", "Alt_Allelic_Freq", "NbBarc_tot", "NbAlns_tot", "NbBarc_allele0", "NbBarc_allele1", "NbBarc_alleleNA", "NbAlns_allele0", "NbAlns_allele1", "NbAlns_alleleNA", "NbBarc_adjLeft", "NbBarc_adjRight", "NbBarc_nodeSVbegin", "NbBarc_nodeSVend", "NbAlns_adjLeft", "NbAlns_adjRight", "NbAlns_nodeSVbegin", "NbAlns_nodeSVend", "MinOccAlns_adjLeft", "MinOccAlns_adjRight", "MinOccAlns_nodeSVbegin", "MinOccAlns_nodeSVend", "MaxOccAlns_adjLeft", "MaxOccAlns_adjRight", "MaxOccAlns_nodeSVbegin", "MaxOccAlns_nodeSVend", "MeanOccAlns_adjLeft", "MeanOccAlns_adjRight", "MeanOccAlns_nodeSVbegin", "MeanOccAlns_nodeSVend", "MedianOccAlns_adjLeft", "MedianOccAlns_adjRight", "MedianOccAlns_nodeSVbegin", "MedianOccAlns_nodeSVend"])+"\n")

    #########################
    #B. Process aln results.
    #########################
    with open(inputGAF, "r", encoding='UTF-8') as gafFile:
        for line in gafFile:
            
            readID, readLen, __, __, __, path, __, pos_start, pos_end, __, alnLen, mapq, *__ = line.split("\t")
            #readID, __, __, readLen, __, __, __, path, __, pos_start, pos_end, __, alnLen, mapq, *__ = line.split("\t")
            
            #1. Filters to keep only the valid alignments.
            ##############################################
            if path == "*":             #remove unmapped reads
                continue
            cov = int(alnLen) / int(readLen)
            if cov < 0.9:
                continue
            if int(mapq) < 30:
                continue

            #2. Get the barcode ID.
            #######################
            if "BX:Z:" in readID:
                barcodeID = "BX:Z:" + ''.join(readID.split("BX:Z:")[1]).split(" ")[0]
                readID = readID.split("BX:Z:")[0]
            else:
                barcodeID = ""


            #3. Get a set of SV regions where at least one vgNode of the path belongs to.
            #############################################################################
            
            list_way_node = extract_nodes(path)           
            for way, node in list_way_node:
                if way == "forward" :
                    if node[0] in gfaNode2svRegionsDict.keys() :
                        list_of_sv_regiontype_coords = gfaNode2svRegionsDict[node[0]]            
                    else :
                        continue

                    for sv_regiontype_coords in list_of_sv_regiontype_coords :
                        sv, region_type, coords_region, node_length = sv_regiontype_coords
                        if ((int(coords_region[0])< int(pos_start)) and (int(coords_region[1]) > int(pos_end))) :

                            if region_type == "adjLeft" :
                                sv.adjLeft.addBarcode(barcodeID)
                            elif region_type == "adjRight":
                                sv.adjRight.addBarcode(barcodeID)
                            elif region_type == "nodeSVbegin" :
                                sv.nodeSVbegin.addBarcode(barcodeID)
                            elif region_type == "nodeSVend" :
                                sv.nodeSVend.addBarcode(barcodeID)

                elif way == "backward" :
                    if node[0] in gfaNode2svRegionsDict.keys() :
                        list_of_sv_regiontype_coords = gfaNode2svRegionsDict[node[0]]
                    else :
                        continue
                    for sv_regiontype_coords in list_of_sv_regiontype_coords :
                        sv, region_type, coords_region, node_length = sv_regiontype_coords
                        start_aln = int(node_length)-int(pos_start)
                        end_aln = int(node_length)-int(pos_end)
                        #Comme on est en backward, l'annotation et comme si le 0 est la fin du noeud
                        if ((int(coords_region[0])< int(end_aln)) and (int(coords_region[1]) > int(start_aln))) :
                            #Comme on est en backward start et end inverse
                            if region_type == "adjLeft" :
                                sv.adjLeft.addBarcode(barcodeID)
                            elif region_type == "adjRight":
                                sv.adjRight.addBarcode(barcodeID)
                            elif region_type == "nodeSVbegin" :
                                sv.nodeSVbegin.addBarcode(barcodeID)
                            elif region_type == "nodeSVend" :
                                sv.nodeSVend.addBarcode(barcodeID)

                # TODO : take into account the information split-reads

    ###########################
    #C. Estimate the genotype.
    ###########################
    # start = time()
    with open(inputVCF, "r", encoding='UTF-8') as inVCF, open(outputVCF, "w", encoding='UTF-8') as outVCF:
        sv_id = 0
        for line in inVCF:
            if line.startswith("##"):
                outVCF.write(line)

            elif line.startswith("#C"):
                outVCF.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                outVCF.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Cumulated depth accross samples (sum)">\n')
                outVCF.write('##FORMAT=<ID=AD,Number=3,Type=Integer,Description="Depth of each allele by sample (allele0, allele1, alleleNA)">\n')
                outVCF.write('##FORMAT=<ID=Alt Allele Freq,Number=1,Type=Integer,Description="Alternative allelic frequency">\n')
                outVCF.write(line.rstrip("\n") + "\t" + "\t".join(["FORMAT", "SAMPLE"]) + "\n")

            else:
                # For each SV of the input VCF file, get the corresponding SV object.
                sv_id += 1

                sv = svsDict[sv_id]
                # Analyse the barcodes signal and estimate the genotype.
                if sv.type == "INV":

                    #NB. Corresponds to Method 3.

                    # Clean the set of barcodes to keep only the informative ones.
                    adjLeft_barcodesDict, adjRight_barcodesDict, nodeSVbegin_barcodesDict, nodeSVend_barcodesDict = cleanBarcodes(sv)
                    
                    # Estimate the genotype: on #alns with only alns of barcodes specific to one allele.
                    nbBarc_support_alleles, nbAlns_support_alleles = getSupportBarcodes(adjLeft_barcodesDict, adjRight_barcodesDict, nodeSVbegin_barcodesDict, nodeSVend_barcodesDict)
                    result_GT, allelic_frequency_allele1 = genotype(nbAlns_support_alleles)
                    
                    # Get statistics.
                    nbBarc_total, nbAlns_total, nbBarc_adjLeft, nbBarc_adjRight, nbBarc_nodeSVbegin, nbBarc_nodeSVend, nbAlns_adjLeft, nbAlns_adjRight, nbAlns_nodeSVbegin, nbAlns_nodeSVend = get_statistics(nbBarc_support_alleles, nbAlns_support_alleles, adjLeft_barcodesDict, adjRight_barcodesDict, nodeSVbegin_barcodesDict, nodeSVend_barcodesDict)

                    # Write the analysis results to the 'analysisFile'.
                    analysisFile.write("\t".join([str(sv.id), "3", str(result_GT), str(allelic_frequency_allele1), str(nbBarc_total), str(nbAlns_total), str(nbBarc_support_alleles[0]), str(nbBarc_support_alleles[1]), str(nbBarc_support_alleles[2]), str(nbAlns_support_alleles[0]), str(nbAlns_support_alleles[1]), str(nbAlns_support_alleles[2]), str(nbBarc_adjLeft), str(nbBarc_adjRight), str(nbBarc_nodeSVbegin), str(nbBarc_nodeSVend), str(nbAlns_adjLeft), str(nbAlns_adjRight), str(nbAlns_nodeSVbegin), str(nbAlns_nodeSVend)])+"\n")

                    # Output the genotype in the output VCF file.
                    numbers = ",".join(str(y) for y in nbAlns_support_alleles)
                    if len(line.split("\t")) <= 8:
                        new_line = (
                            line.rstrip("\n")
                            + "\t"
                            + "GT:DP:AD:AF"
                            + "\t"
                            + result_GT
                            + ":"
                            + str(round(sum(nbAlns_support_alleles), 3))
                            + ":"
                            + str(numbers)
                            + ":"
                            + str(allelic_frequency_allele1)
                        )
                        outVCF.write(new_line + "\n")
                    else:
                        line_without_genotype = line.split("\t")[0:8]
                        new_line = (
                            "\t".join(line_without_genotype)
                            + "\t"
                            + "GT:DP:AD:AF"
                            + "\t"
                            + result_GT
                            + ":"
                            + str(round(sum(nbAlns_support_alleles), 3))
                            + ":"
                            + str(numbers)
                            + ":"
                            + str(allelic_frequency_allele1)
                        )
                        outVCF.write(new_line + "\n")

        analysisFile.close()


#############
# Functions.
#############

def create_region(sv, node, orientation, region_size, region_type, gfaNode2svRegionsDict, gfa_graph):
    ''' Function to create regions based on node size and set region size '''

    # If the node is smaller than the set region size, then create the region using a deep graph traversal
    if length_node(node) < region_size :
        if region_type == 'nodeSVbegin' or region_type == 'nodeSVend' :
            if region_size > int(sv.length / 2):
                regionSize_nodeSV = int(sv.length / 2)
                dico_dfs_region = createRegion_DFS(node,orientation,regionSize_nodeSV,gfa_graph)
            else :
                dico_dfs_region = createRegion_DFS(node,orientation,region_size, gfa_graph)
        else :
            dico_dfs_region = createRegion_DFS(node,orientation,region_size, gfa_graph)
        
        dico_dfs_region = clean_region(dico_dfs_region)
        format_region(sv,dico_dfs_region,region_type,gfaNode2svRegionsDict)

    # Otherwise create a region on the node
    else :
        associate_GFANode_To_SVRegion(sv, node,region_type, region_size, gfaNode2svRegionsDict)


def createRegion_DFS(node, orientation, region_size, gfa_graph):
    ''' Function that allows a graph to be traversed in depth to create a region based on a fixed size. '''

    visited = {}    # We keep the node and the orientation in which it was traversed as a function of the distance covered to reach this node.
    stack = deque()
    region ={}

    stack.append((node,0,orientation))   # Add to stack: node, distance traveled to reach this node (first node therefore 0), orientation of entry into node
    visited[f"{node}-{Orientation.FORWARD}"] = 0    
    visited[f"{node}-{Orientation.REVERSE}"] = 0  

    while stack:
        node, dist, ori = stack.pop()   # dist = "distance covered to reach this node (i.e. sum of nodes covered to reach it)"
        dist_path = 0                   # dist_path = "distance covered to reach the "end" of the node (i.e. previous nodes + node size)"
        visited[f"{node}-{ori}"]=dist
        dist_path = dist + length_node(node)

        if dist_path < region_size :    # Region size = "fixed region size"
            region[node] = []
            region[node].append(('Full',[0,length_node(node)]))

            for info_succ in gfa_graph.segments[node]["out"][ori] : # gfa_graph.segments[node]["out"][ori] -> returns the successor nodes (with outgoing edges) of the node of interest, as well as their direction of entry into the node.
                succ = info_succ[0]
                ori_succ = info_succ[1]     # Ori_succ =direction of entry into the successor node (forward ou reverse)

                if f"{succ}-{ori_succ}" in visited.keys():  # If the node has already been visited in this orientation, we check whether the current route gives a smaller distance, and if so, we keep this distance (because we want to "go as far as possible").
                    if  visited[f"{succ}-{ori_succ}"] > dist_path : # if dist successor > dist_path
                        stack.append((succ,dist_path,ori_succ))
                        visited[f"{succ}-{ori_succ}"] = dist_path

                else :
                    stack.append((succ,dist_path,ori_succ)) 
                    visited[f"{succ}-{ori_succ}"] = dist_path

        elif dist_path > region_size :  # If we reach a size larger than the set region size, we take the part of the node that allows us to reach the set size, according to its reading direction.
            if node not in region.keys():
                region[node] = []
            if ori ==  Orientation.FORWARD:
                region[node].append(('CutF',[0,region_size-dist]))
            elif ori == Orientation.REVERSE:
                region[node].append(('CutR',[length_node(node)-(region_size-dist),length_node(node)]))
        else : #dist_path = region_size
            region[node] = []
            region[node].append(('Full',[0,region_size-dist]))
    return region


def length_node(node):
    ''' Function to retrieve node size from node identifier '''
    node_start = int(str(node).split(":")[1]) - 1            
    node_end = int(str(node).split(":")[2]) 
    node_length = (node_end-node_start)
    return node_length


def clean_region(region):
    ''' Function that processes nodes where several regions read in both orientations are present and checks that the regions do not overlap. '''
    for node in region.keys():
        if len(region[node]) > 1 :
            #print(region, "\n")
            size_forward = 0
            size_reverse = 0
            for i in range(2):
                if region[node][i][0] == "CutF" :
                    if size_forward < region[node][i][1][1] - region[node][i][1][0] :
                        size_forward = region[node][i][1][1] - region[node][i][1][0]
                if region[node][i][0] == "CutR" :
                    if size_reverse < region[node][i][1][1] - region[node][i][1][0] :
                        size_reverse = region[node][i][1][1] - region[node][i][1][0]

            if (size_forward + size_reverse) >= length_node(node) :
                region[node]=[('Full',[0,length_node(node)])]

            if size_reverse == 0 :
                region[node]=[('CutF',[0,size_forward])]
            elif size_forward == 0 :
                region[node]=[('CutR',[0,size_reverse])]
    return region

def format_region(sv,region_dico, region_type,gfaNode2svRegionsDict):
    ''' Function to switch to the expected format from dico_region post DFS '''
    for node in region_dico.keys():
        if node not in gfaNode2svRegionsDict.keys():
            gfaNode2svRegionsDict[node] = []

        node_lenght = length_node(node) 
        if len(region_dico[node]) > 1 :
            for i in range(2):
                coords = region_dico[node][i][1]
                if region_type == "adjLeft":
                    sv.adjLeft = sv.getAdjLeft(coords,node)
                elif region_type == "adjRight":
                    sv.adjRight = sv.getAdjRight(coords,node)
                elif region_type == "nodeSVbegin":
                    sv.nodeSVbegin = sv.getNodeSVbegin(coords,node)
                elif region_type == "nodeSVend":
                    sv.nodeSVend = sv.getNodeSVend(coords,node)

                gfaNode2svRegionsDict[node].append((sv, region_type,coords,node_lenght))
        else :
            coords = region_dico[node][0][1]
            if region_type == "adjLeft":
                sv.adjLeft = sv.getAdjLeft(coords,node)
            elif region_type == "adjRight":
                sv.adjRight = sv.getAdjRight(coords,node)
            elif region_type == "nodeSVbegin":
                sv.nodeSVbegin = sv.getNodeSVbegin(coords,node)
            elif region_type == "nodeSVend":
                sv.nodeSVend = sv.getNodeSVend(coords,node)

            gfaNode2svRegionsDict[node].append((sv, region_type,coords,node_lenght))

     
def associate_GFANode_To_SVRegion(sv_object, gfaNode, region_type, regionSize, gfaNode2svRegionsDict):
    """Method to associate a GFA node to a SV region."""
    
    node_start = int(str(gfaNode).split(":")[1]) - 1             #positions in 'gfaNode_id' are 1-based and incl.
    node_end = int(str(gfaNode).split(":")[2])  
    node_length = (node_end-node_start)     #'node_start' and 'node_end' are 0-based and incl./excl. resp.
    
    # adjLeft.
    if region_type == "adjLeft":
        coords = [(node_length-regionSize), node_length]
        sv_object.adjLeft = sv_object.getAdjLeft(coords,gfaNode) 


    # adjRight.
    elif region_type == "adjRight":
        coords = [0, regionSize]
        sv_object.adjRight = sv_object.getAdjRight(coords,gfaNode)

    # nodeSVbegin.
    elif region_type == "nodeSVbegin":
        if regionSize > int(sv_object.length / 2):
            regionSize_nodeSV = int(sv_object.length / 2)
        else:
            regionSize_nodeSV = regionSize
        coords = [0, regionSize_nodeSV]
        sv_object.nodeSVbegin = sv_object.getNodeSVbegin(coords,gfaNode)

    # nodeSVend.
    elif region_type == "nodeSVend":
        if regionSize > int(sv_object.length / 2):
            regionSize_nodeSV = int(sv_object.length / 2)
        else:
            regionSize_nodeSV = regionSize
        coords = [(node_length-regionSize_nodeSV), node_length]
        sv_object.nodeSVend = sv_object.getNodeSVend(coords,gfaNode)

    if gfaNode not in gfaNode2svRegionsDict:
        gfaNode2svRegionsDict[gfaNode] = [(sv_object, region_type, coords,node_length)]
    else:
        gfaNode2svRegionsDict[gfaNode].append((sv_object, region_type, coords,node_length))

def extract_nodes(path):       
    """Method to extract the nodes contained in a path from an alignment GAF file."""                                        
    list_way_node = []
 
    if path.startswith('>'):
        way_ALN = "forward"
    elif path.startswith('<'):
        way_ALN = "backward"

    node = [node for node in re.split(r'[<>]', path) if node]
    
    list_way_node.append((way_ALN, node))    
    
    return list_way_node

def genotype(nbAlnBarc_support_alleles):
    """Method to return the genotype of the SV using the 'allelic_frequency'"""

	# Allelic frequency.
	####################
    if nbAlnBarc_support_alleles[0] == 0 and nbAlnBarc_support_alleles[1] == 0 :
        result_GT = './.'
        allelic_frequency_allele1 = "NA"    
    else :
        allelic_frequency_allele1 = nbAlnBarc_support_alleles[1] / (nbAlnBarc_support_alleles[0] + nbAlnBarc_support_alleles[1])
        if allelic_frequency_allele1 >= 0.8:	# allelic_frequency_allele1 close to 1 --> supports allele 1.
            result_GT = "1/1"
        elif allelic_frequency_allele1 <= 0.2:	# allelic_frequency_allele1 close to 0 --> supports allele 0.
            result_GT = "0/0"
        else:                                   # allelic_frequency_allele1 close to 0.5 --> supports both alleles (heterozygous).
            result_GT = "0/1"

    return result_GT, allelic_frequency_allele1


def cleanBarcodes(sv_object):
    """Method to clean the set of barcodes mapping on the regions of the current SV."""
    # Remove barcodes mapping on both 'adjLeft' and 'adjRight'.

    if hasattr(sv_object.adjLeft, 'barcodesDict') and sv_object.adjLeft.barcodesDict:
        adjLeft_barcodesDict = sv_object.adjLeft.barcodesDict.copy()
    else : 
        adjLeft_barcodesDict = {}

    if hasattr(sv_object.adjRight, 'barcodesDict') and sv_object.adjRight.barcodesDict:
        adjRight_barcodesDict = sv_object.adjRight.barcodesDict.copy()
    else :
        adjRight_barcodesDict = {}


    if hasattr(sv_object.nodeSVbegin, 'barcodesDict') and sv_object.nodeSVbegin.barcodesDict :
        nodeSVbegin_barcodesDict = sv_object.nodeSVbegin.barcodesDict.copy()
    else :
        nodeSVbegin_barcodesDict = {}

    if hasattr(sv_object.nodeSVend, 'barcodesDict') and sv_object.nodeSVend.barcodesDict:
        nodeSVend_barcodesDict = sv_object.nodeSVend.barcodesDict.copy()
    else :
        nodeSVend_barcodesDict = {}

    for barcode in list(adjLeft_barcodesDict.keys()):               #pylint: disable=consider-using-dict-items
        if barcode in adjRight_barcodesDict:
            del adjLeft_barcodesDict[barcode]
            del adjRight_barcodesDict[barcode]

    return adjLeft_barcodesDict, adjRight_barcodesDict, nodeSVbegin_barcodesDict, nodeSVend_barcodesDict


def getSupportBarcodes(adjLeft_barcodesDict, adjRight_barcodesDict, nodeSVbegin_barcodesDict, nodeSVend_barcodesDict):
    """
    Method to get the number of alignments (and unique barcodes) that support the allele 0, and the number of alignments (and unique barcodes) that support the allele 1.
    NB: Allele 0 represented by junctions ('adjLeft' and 'nodeSVbegin') and ('nodeSVend' and 'adjRight') for INV.
    NB: Allele 1 represented by junctions ('adjLeft' and 'nodeSVend') and ('nodeSVbegin' and 'adjRight') for INV.
    """
    nbBarc_support_allele0 = 0
    nbBarc_support_allele1 = 0
    nbBarc_undetermined = 0
    nbAln_support_allele0 = 0
    nbAln_support_allele1 = 0
    nbAln_undetermined = 0

    #NB: Count only the barcodes specific to one allele.
    #NB: Count only the alns of barcodes specific to one allele.
    for barcode in list(adjLeft_barcodesDict.keys()):
        # Allele 0.
        if (barcode in nodeSVbegin_barcodesDict) and (barcode not in nodeSVend_barcodesDict):
            nbBarc_support_allele0 += 1
            nbAln_support_allele0 += adjLeft_barcodesDict[barcode] + nodeSVbegin_barcodesDict[barcode]
        # Allele 1.
        elif (barcode in nodeSVend_barcodesDict) and (barcode not in nodeSVbegin_barcodesDict):
            nbBarc_support_allele1 += 1
            nbAln_support_allele1 += adjLeft_barcodesDict[barcode] + nodeSVend_barcodesDict[barcode]
        # Undetermined.
        elif (barcode in nodeSVbegin_barcodesDict) and (barcode in nodeSVend_barcodesDict):
            nbBarc_undetermined += 1
            nbAln_undetermined += adjLeft_barcodesDict[barcode] + nodeSVbegin_barcodesDict[barcode] + nodeSVend_barcodesDict[barcode]
            
    for barcode in list(adjRight_barcodesDict.keys()):
        # Allele 0.
        if (barcode in nodeSVend_barcodesDict) and (barcode not in nodeSVbegin_barcodesDict):
            nbBarc_support_allele0 += 1
            nbAln_support_allele0 += adjRight_barcodesDict[barcode] + nodeSVend_barcodesDict[barcode]
        # Allele 1.
        elif (barcode in nodeSVbegin_barcodesDict) and (barcode not in nodeSVend_barcodesDict):
            nbBarc_support_allele1 += 1
            nbAln_support_allele1 += adjRight_barcodesDict[barcode] + nodeSVbegin_barcodesDict[barcode]
        # Undetermined.
        elif (barcode in nodeSVend_barcodesDict) and (barcode in nodeSVbegin_barcodesDict):
            nbBarc_undetermined += 1
            nbAln_undetermined += adjRight_barcodesDict[barcode] + nodeSVbegin_barcodesDict[barcode] + nodeSVend_barcodesDict[barcode]

    for barcode in list(nodeSVbegin_barcodesDict.keys()):
        # Undetermined.
        if (barcode in nodeSVend_barcodesDict) and (barcode not in adjLeft_barcodesDict) and (barcode not in adjRight_barcodesDict):
            nbBarc_undetermined += 1
            nbAln_undetermined += nodeSVbegin_barcodesDict[barcode] + nodeSVend_barcodesDict[barcode]
    
    return [nbBarc_support_allele0, nbBarc_support_allele1, nbBarc_undetermined], [nbAln_support_allele0, nbAln_support_allele1, nbAln_undetermined]


def get_statistics(nbBarc_support_alleles, nbAlns_support_alleles, adjLeft_barcodesDict, adjRight_barcodesDict, nodeSVbegin_barcodesDict, nodeSVend_barcodesDict):
    """Method to get statistics on the results of the genotype estimation."""
    nbBarc_total = sum(nbBarc_support_alleles)
    nbAlns_total = sum(nbAlns_support_alleles)
    nbBarc_adjLeft = len(list(adjLeft_barcodesDict.keys()))
    nbBarc_adjRight = len(list(adjRight_barcodesDict.keys()))
    nbBarc_nodeSVbegin = len(list(nodeSVbegin_barcodesDict.keys()))
    nbBarc_nodeSVend = len(list(nodeSVend_barcodesDict.keys()))
    nbAlns_adjLeft = sum(list(adjLeft_barcodesDict.values()))
    nbAlns_adjRight = sum(list(adjRight_barcodesDict.values()))
    nbAlns_nodeSVbegin = sum(list(nodeSVbegin_barcodesDict.values()))
    nbAlns_nodeSVend = sum(list(nodeSVend_barcodesDict.values()))

    return nbBarc_total, nbAlns_total, nbBarc_adjLeft, nbBarc_adjRight, nbBarc_nodeSVbegin, nbBarc_nodeSVend, nbAlns_adjLeft, nbAlns_adjRight, nbAlns_nodeSVbegin, nbAlns_nodeSVend #, minOccPerBarcode_adjLeft, minOccPerBarcode_adjRight, minOccPerBarcode_nodeSVbegin, minOccPerBarcode_nodeSVend, maxOccPerBarcode_adjLeft, maxOccPerBarcode_adjRight, maxOccPerBarcode_nodeSVbegin, maxOccPerBarcode_nodeSVend, meanOccPerBarcode_adjLeft, meanOccPerBarcode_adjRight, meanOccPerBarcode_nodeSVbegin, meanOccPerBarcode_nodeSVend, medianOccPerBarcode_adjLeft, medianOccPerBarcode_adjRight, medianOccPerBarcode_nodeSVbegin, medianOccPerBarcode_nodeSVend
 


##############################################
if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")

    else:
        main(sys.argv[1:])

