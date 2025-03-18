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
Module 'classes_creationGFA.py': Classes 'Chrom', 'SV'.
"""

# import os
# import re
# import sys
# from Bio import SeqIO

#pylint: disable=line-too-long, disable=trailing-whitespace, disable=consider-using-f-string

#----------------------------------------------------
# 'Chrom' class
#----------------------------------------------------
class Chrom:
    """
    Class defining a chromosome characterized by:
    - its ID
    - its sequence
    - its length
    - a dictionary of SV breakpoints present on this chromosome
        key: breakpoint coordinate
        value: [list of SV objects corresponding to this breakpoint]
    - a list of SVs (SV objects) present on this chromosome
    - a dictionary of GFA nodes characterizing this chromosome
        key: GFA node ID
        value: GFA node length
    """

    #Constructor
    def __init__(self, chr_id, sequence):
        self.id = chr_id
        self.sequence = sequence
        self.length = len(self.sequence)
        self.bkpts = {}
        self.svs = []
        self.gfaNodes = {}

    #Method "__str__"
    def __str__(self):
        return f"{self.id}"

    #Method "addBreakpoint"
    def addBreakpoint(self, bkpt_coord, sv_object):
        """Method to add a SV breakpoint present on this chromosome."""
        if int(bkpt_coord) not in self.bkpts:
            self.bkpts[int(bkpt_coord)] = [sv_object]
        else:
            self.bkpts[int(bkpt_coord)].append(sv_object)

    #Method "addSV"
    def addSV(self, sv_object):
        """Method to add a SV present on this chromosome."""
        self.svs.append(sv_object)

    #Method "writeInfoSVLines"
    def writeInfoSVLines(self, headerLines):
        """Method to append headerLines with a line containing the chrID and the list of SVs for this chromosome."""
        headerLines.append("#{}\t{}\n".format(self.id, ";".join([sv.format for sv in self.svs])))

    #Method "addGFANode"
    def addGFANode(self, gfaNode_id):
        """Method to add a GFA node characterizing this chromosome."""
        node_start = int(str(gfaNode_id).split(":")[1]) - 1             #positions in 'gfaNode_id' are 1-based and incl.
        node_end = int(str(gfaNode_id).split(":")[2])                   #'node_start' and 'node_end' are 0-based and incl./excl. resp.
        node_length = node_end - node_start
        self.gfaNodes[gfaNode_id] = node_length

    #Method "writeSLines"
    def writeSLines(self, SLines):
        """Method to append SLines with a segment line containing the GFA node and its sequence."""
        for gfaNode in self.gfaNodes:
            node_start = int(gfaNode.split(":")[1]) - 1                 #positions in 'gfaNode_id' are 1-based and incl.
            node_end = int(gfaNode.split(":")[2])                       #'node_start' and 'node_end' are 0-based and incl./excl. resp.
            SLines.append("\t".join(["S", gfaNode, self.sequence[node_start:node_end]]) + "\n")

    #Method "writeRefPLine"
    def writeRefPLine(self, PLines):
        """Method to append PLines with a path line containing the reference path."""
        path_nodes = "+,".join(self.gfaNodes.keys()) + "+"
        gfaNodes_values = [str(value) for value in self.gfaNodes.values()]
        path_lengths = "M,".join(gfaNodes_values) + "M"
        PLines.append("\t".join(["P", self.id, path_nodes, path_lengths]) + "\n")

    #Method "writeRefLLines"
    def writeRefLLines(self, LLines, leftEdges, rightEdges):
        """
        Method to append LLines with a link line containing a link between two GFA nodes with their respective orientations ;
        and to update the dictionary 'leftEdges' and 'rightEdges' containing the neighbors nodes (with their orientation) of each GFA node.
        """
        gfaNodes_keys = list(self.gfaNodes.keys())
        for i in range(len(gfaNodes_keys)-1):
            LLines.append("\t".join(["L", gfaNodes_keys[i], "+", gfaNodes_keys[i+1], "+", "0M"]) + "\n")
            #Ajout de ligne - A REMETTRE APREs TEST
            LLines.append("\t".join(["L", gfaNodes_keys[i+1], "-", gfaNodes_keys[i], "-", "0M"]) + "\n")
            
            if gfaNodes_keys[i+1] not in leftEdges:
                leftEdges[gfaNodes_keys[i+1]] = {"+": [], "-": []}
                rightEdges[gfaNodes_keys[i+1]] = {"+": [], "-": []}
            if gfaNodes_keys[i] not in rightEdges:
                rightEdges[gfaNodes_keys[i]] = {"+": [], "-": []}
                leftEdges[gfaNodes_keys[i]] = {"+": [], "-": []}

            leftEdges[gfaNodes_keys[i+1]]["+"].append((gfaNodes_keys[i], "+"))
            rightEdges[gfaNodes_keys[i]]["+"].append((gfaNodes_keys[i+1], "+"))
       
    #Method "writeAltLLines"
    def writeAltLLines(self, LLines, sv_object, leftEdges, rightEdges, SLines, dict_ins_seq):
        """
        Method to append LLines with alternative link lines according to the SV type ;
        and to update the dictionary 'leftEdges' and 'rightEdges' containing the neighbors nodes (with their orientation) of each GFA node.
        """
        b1_left, b1_right, b2_left, b2_right = sv_object.gfaNodes

        # INVersion.
        if sv_object.type == "INV":
            LLines.append("\t".join(["L", b1_left, "+", b2_left, "-", "0M"]) + "\n")
            LLines.append("\t".join(["L", b1_right, "-", b2_right, "+", "0M"]) + "\n")
            rightEdges[b1_left]["+"].append((b2_left, "-"))
            leftEdges[b2_left]["-"].append((b1_left, "+"))
            leftEdges[b2_right]["+"].append((b1_right, "-"))
            rightEdges[b1_right]["-"].append((b2_right, "+"))

            # #Ajout des lignes - A REMETTRE APRES TEST
            LLines.append("\t".join(["L", b2_left, "+", b1_left, "-", "0M"]) + "\n")
            LLines.append("\t".join(["L", b2_right, "-", b1_right, "+", "0M"]) + "\n")
            rightEdges[b2_left]["+"].append((b1_left, "-"))
            leftEdges[b1_left]["-"].append((b2_left, "+"))
            leftEdges[b1_right]["+"].append((b2_right, "-"))
            rightEdges[b2_right]["-"].append((b1_right, "+"))

        # DELetion.
        elif sv_object.type == "DEL":
            #TODO: check b1_right = b2_left, otherwise nested variant or issue
            LLines.append("\t".join(["L", b1_left, "+", b2_right, "+", "0M"]) + "\n")
            rightEdges[b1_left]["+"].append((b2_right, "+"))
            leftEdges[b2_right]["+"].append((b1_left, "+"))

        # INSertion.
        elif sv_object.type == "INS":

            # Coordinates of INS sequence end with "a" (for alternative).
            ins_node = ":".join([self.id, str(sv_object.coords[0])+"a"])
            SLines.append("\t".join(["S", ins_node, dict_ins_seq[sv_object.id]]) + "\n")
            sv_object.gfaNodes[1] = ins_node
            sv_object.gfaNodes[-2] = ins_node
            
            LLines.append("\t".join(["L", b1_left, "+", ins_node, "+", "0M"]) + "\n")
            LLines.append("\t".join(["L", ins_node, "+", b2_right, "+", "0M"]) + "\n")
            rightEdges[b1_left]["+"].append((ins_node, "-"))
            leftEdges[ins_node]["-"] = [(b1_left, "+")]
            leftEdges[b2_right]["+"].append((ins_node, "-"))
            rightEdges[ins_node]["-"] = [(b2_right, "+")]


#----------------------------------------------------
# 'SV' class
#----------------------------------------------------
class SV:
    """
    Class defining a SV characterized by:
    - its ID
    - its type
    - the chrom ID it is present on
    - its coords on the chrom (`[<pos>, <end>]`, 1-based with <pos> incl. and <end> excl.)
    - its format: `<sv_type>-<pos>-<end>`
    - its length (<end> - <pos>)
    - a list of GFA nodes IDs that are adjacents to the breakpoints of the SV: [b1.left, b1.right, b2.left, b2.right]
        with b1.left: GFA node ID adjacent to the left of the first bkpt
             b1.right: GFA node ID adjacent to the right of the first bkpt
             b1.left: GFA node ID adjacent to the left of the second bkpt
             b1.right: GFA node ID adjacent to the right of the second bkpt
    - the 4 regions corresponding to this SV:
        'adjLeft', 'nodeSVbegin', 'nodeSVend', 'adjRight'
    """
    svsList = []

    #Constructor
    def __init__(self, sv_id, sv_type, chrom, coords, sv_format, sv_length):
        self.id = sv_id
        self.type = sv_type
        self.chrom = chrom
        self.coords = coords
        self.format = sv_format
        self.length = sv_length
        self.gfaNodes = []
        self.adjLeft = None
        self.adjRight = None
        self.nodeSVbegin = None
        self.nodeSVend = None
        # self.adjLeft = self.Region(self.id, "adjLeft")
        # self.adjRight = self.Region(self.id, "adjRight")
        # self.nodeSVbegin = self.Region(self.id, "nodeSVbegin")
        # self.nodeSVend = self.Region(self.id, "nodeSVend")
        SV.svsList.append(self)
        
    #Method "__str__"
    def __str__(self):
        return f"{self.format}"

    def getAdjLeft(self, coords, node):
        """Method to return the adjLeft region of this SV."""
        if self.adjLeft is None:
            self.adjLeft = self.Region(self.id, "adjLeft", coords, node)
        return self.adjLeft

    def getAdjRight(self, coords, node):
        """Method to return the adjRight region of this SV."""
        if self.adjRight is None:
            self.adjRight = self.Region(self.id, "adjRight", coords, node)
        return self.adjRight
    
    def getNodeSVbegin(self, coords, node):
        """Method to return the nodeSVbegin region of this SV."""
        if self.nodeSVbegin is None:
            self.nodeSVbegin = self.Region(self.id, "nodeSVbegin", coords, node)
        return self.nodeSVbegin

    def getNodeSVend(self, coords, node):
        """Method to return the nodeSVend region of this SV."""
        if self.nodeSVend is None:
            self.nodeSVend = self.Region(self.id, "nodeSVend", coords, node)
        return self.nodeSVend

    def cleanBarcodes(self):
        """Method to clean the set of barcodes mapping on the regions of the current SV."""
        # # Remove barcodes mapping on both 'nodeSVbegin' and 'nodeSVend'.
        # ##--> actually, maybe keep because we can see a #ce if we count #reads / barcode
        # for barcode in list(self.nodeSVbegin.barcodesDict.keys()):          #pylint: disable=consider-using-dict-items
        #     if barcode in self.nodeSVend.barcodesDict:
        #         del self.nodeSVbegin.barcodesDict[barcode]
        #         del self.nodeSVend.barcodesDict[barcode]
                
        # Remove barcodes mapping on both 'adjLeft' and 'adjRight'.
        for barcode in list(self.adjLeft.barcodesDict.keys()):               #pylint: disable=consider-using-dict-items
            if barcode in self.adjRight.barcodesDict:
                del self.adjLeft.barcodesDict[barcode]
                del self.adjRight.barcodesDict[barcode]

        # # Remove barcodes from regions where it is seen only once.
        # ## adjLeft.
        # for barcode in list(self.adjLeft.barcodesDict.keys()):
        #     if self.adjLeft.barcodesDict[barcode] == 1:
        #         del self.adjLeft.barcodesDict[barcode]
        # ## adjRight.
        # for barcode in list(self.adjRight.barcodesDict.keys()):
        #     if self.adjRight.barcodesDict[barcode] == 1:
        #         del self.adjRight.barcodesDict[barcode]
        # ## nodeSVbegin.
        # for barcode in list(self.nodeSVbegin.barcodesDict.keys()):
        #     if self.nodeSVbegin.barcodesDict[barcode] == 1:
        #         del self.nodeSVbegin.barcodesDict[barcode]
        # ## nodeSVend.
        # for barcode in list(self.nodeSVend.barcodesDict.keys()):
        #     if self.nodeSVend.barcodesDict[barcode] == 1:
        #         del self.nodeSVend.barcodesDict[barcode]

    def getSupportBarcodes(self):
        """
        Method to get the number of alignments that support the allele 0, and the number of alignments that support the allele 1.
        NB: Allele 0 represented by junctions ('adjLeft' and 'nodeSVbegin') and ('nodeSVend' and 'adjRight') for INV.
        NB: Allele 1 represented by junctions ('adjLeft' and 'nodeSVend') and ('nodeSVbegin' and 'adjRight') for INV.
        """
        nbAln_support_allele0 = 0
        nbAln_support_allele1 = 0
        nbAln_undetermined = 0
        nbBarc_support_allele0 = 0
        nbBarc_support_allele1 = 0
        nbBarc_undetermined = 0
        #TODO: for other SV types (for now only for INV).

        for barcode, occ in self.adjLeft.barcodesDict.items():
            # Allele 0.
            if (barcode in self.nodeSVbegin.barcodesDict) and (barcode not in self.nodeSVend.barcodesDict):
                nbAln_support_allele0 += self.adjLeft.barcodesDict[barcode] + self.nodeSVbegin.barcodesDict[barcode]
                nbBarc_support_allele0 += 1
            # Allele 1.
            elif (barcode in self.nodeSVend.barcodesDict) and (barcode not in self.nodeSVbegin.barcodesDict):
                nbAln_support_allele1 += self.adjLeft.barcodesDict[barcode] + self.nodeSVend.barcodesDict[barcode]
                nbBarc_support_allele1 += 1
            # Undetermined.
            elif (barcode in self.nodeSVbegin.barcodesDict) and (barcode in self.nodeSVend.barcodesDict):
                nbAln_undetermined += self.adjLeft.barcodesDict[barcode] + self.nodeSVbegin.barcodesDict[barcode] + self.nodeSVend.barcodesDict[barcode]
                nbBarc_undetermined += 1

        for barcode, occ in self.adjRight.barcodesDict.items():
            # Allele 0.
            if (barcode in self.nodeSVend.barcodesDict) and (barcode not in self.nodeSVbegin.barcodesDict):
                nbAln_support_allele0 += self.adjRight.barcodesDict[barcode] + self.nodeSVend.barcodesDict[barcode]
                nbBarc_support_allele0 += 1
            # Allele 1.
            elif (barcode in self.nodeSVbegin.barcodesDict) and (barcode not in self.nodeSVend.barcodesDict):
                nbAln_support_allele1 += self.adjRight.barcodesDict[barcode] + self.nodeSVbegin.barcodesDict[barcode]
                nbBarc_support_allele1 += 1
            # Undetermined.
            elif (barcode in self.nodeSVend.barcodesDict) and (barcode in self.nodeSVbegin.barcodesDict):
                nbAln_undetermined += self.adjRight.barcodesDict[barcode] + self.nodeSVend.barcodesDict[barcode] + self.nodeSVbegin.barcodesDict[barcode]
                nbBarc_undetermined += 1

        return [nbAln_support_allele0, nbAln_support_allele1, nbAln_undetermined], [nbBarc_support_allele0, nbBarc_support_allele1, nbBarc_undetermined]

    #----------------------------------------------------
    # 'Region' class (inner class)
    #----------------------------------------------------
    class Region():
        """
        Class defining one of the 4 regions of an SV, characterized by:
        - the corresponding sv ID 
        - a dictionary of the #occurrences of each barcode present on this region
            key: barcode ID
            value: #occurrences of this barcode on this region
        #NB: The regions of an SV are:
            'adjLeft', 'nodeSVbegin', 'nodeSVend', 'adjRight'
        """
        def __init__(self, sv_id, region_type, coords, node):
            self.sv = sv_id
            self.type = region_type
            self.coords = coords
            self.node = node
            self.barcodesDict = {}

        def addBarcode(self, barcodeID):
            """Method to append barcodesDict with a barcodeID (key) and its #occurrences in this region (value)."""
            if barcodeID not in self.barcodesDict:
                self.barcodesDict[barcodeID] = 1
            else:
                self.barcodesDict[barcodeID] += 1

