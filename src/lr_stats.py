#!/usr/bin/env python3

#*****************************************************************************
#  Name: LR_Stats
#  Description: Genotyping of SVs with linked-reads data
#  Copyright (C) 2025 INRIA
#  Author: Mélody Temperville
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


import pysam
import sys
import statistics
#output
from prettytable import PrettyTable
import  matplotlib.pyplot as plt
import pandas as pd

#####################################################################################################
### CLASS

class Barcode:
    def __init__(self, reads=None, unmap=0, not_count_reads=0, mol=None):
        self.reads = reads if reads is not None else []
        self.mol = mol if mol is not None else []
        self.unmap = unmap
        self.not_count_reads = not_count_reads

    def add_read(self, new_reads) :
        self.reads.append(new_reads)

    def add_mol(self, new_mol) :
        self.mol.append(new_mol)

    def count_unmap(self) :
        self.unmap += 1 

    def count_reads(self):
        return len(self.reads)
    
    def count_mol(self):
        return len(self.mol)
    
    def count_not_count_reads(self, nb_reads): #number of reads in molecule with less than 3 reads
        self.not_count_reads += nb_reads

    def nb_read_post_deconvolution(self): #number of reads without the number of reads in molecule with less than 3 reads
        a = len(self.reads) - self.not_count_reads
        return a
    
#####################################################################################################


def add_subparser(subparsers):
    p = subparsers.add_parser("lr-stats", help="Analyze the characteristics of linked-read data")
    p.add_argument( "-b", "--bam", metavar="<sort_bam_file>", help= "Bam file out of a mapping linked-reads/reference and with BX tag", type=str, required=True)
    p.add_argument( "-s", "--molecule_max_size", metavar="<molecule_max_size>", help="Maximum size between two reads to share a same barcode and from the same molecule [default=100000]",type=int, required=False, default=100000)
    p.add_argument( "-G", "--graph_output", metavar="<graphe_output_path/name_file>", help="Path/Name_file for the output graph.png [default=LR_Stats_graph.png]", type=str, required=False, default="LRStats_graph.png")
    p.add_argument( "-o", "--output_table", metavar="<output_table_path/name_file>",help="Path/Name_file for the output table.csv [default=LR_Stats_table.csv] ", type=str, required=False, default="LRStats_table.csv")
    p.add_argument( "-g", "--genome_size", metavar="<genome_size>", help="Genome size required to calculate depth", type=int, required=False, default=0)
    p.add_argument( "-r", "--read_size", metavar="<read_size>", help="Read size required to calculate depth", type=int, required=False, default=150)
    p.set_defaults(func=main, _parser = p)


def main(args):
    if len(sys.argv) == 2:   # no args
        args._parser.print_help()
        sys.exit(1)

    input = args.bam
    mol_max_size = args.molecule_max_size 
    output_histo = args.graph_output
    out_table = args.output_table
    genome_size = args.genome_size
    read_size = args.read_size

    #####################################################################################################
    ### Parsing of bam file to create dictionary {barcode:(start,end)}
    #####################################################################################################

    barcodes = {} # {barcode:(start,end)}
    dico_Barcode = {} # {barcode:barcode_object}
    nb_read = 0
    
    bam_file = pysam.AlignmentFile(input, "rb")
    # Find genome size
    if genome_size == 0 :
        for sq in bam_file.header["SQ"]:
            genome_size = genome_size + sq["LN"]

    print("Genome_size:",genome_size, file = sys.stderr)

    for read in bam_file.fetch():
        nb_read +=1
        barcode = read.get_tag("BX") #Take barcode from the BX tag
        
        if barcode not in dico_Barcode :  #If not already create, creation {barcode:barcode_object}
            dico = False
            obj_bc = Barcode()
        else :  # If already create, edit barcode_object
            dico = True
            obj_bc = dico_Barcode[barcode]

        
        if read.is_unmapped:
            obj_bc.count_unmap()

        else :
            pos = (read.reference_start, read.reference_end)
            chrom = read.reference_name

            obj_bc.add_read(pos)
            if f'{chrom}-@-{barcode}' not in barcodes : 
                barcodes[f'{chrom}-@-{barcode}'] = [pos]
            else :
                barcodes[f'{chrom}-@-{barcode}'].append(pos)


        if dico == False : 
            dico_Barcode[barcode] = obj_bc

    bam_file.close()
  
    #####################################################################################################
    ### Deconvolution (find originate molecules)
    #####################################################################################################

    ### --------------------------------------- Function --------------------------------------- ###  
    def create_mol(num_first_read, num_last_read, mol_max_size, dico_mol, barcodes, barcode, c, not_mol, nb_read_not_count, dico_Barcode, nb_reads_mol):
        '''Function for finding the original molecules in linked-reads data based on barcode information and a maximum possible size for a molecule.'''

        dist_between_first_last = barcodes[barcode][num_last_read][1] - barcodes[barcode][num_first_read][0]  # distance between the start of the first read and the end of the last one (for a barocde)
        #Step 1 : To find the molecule extremity
        while dist_between_first_last > mol_max_size : #while the distance between the first and the  last (or the (n-x)) if lower than the max size for a molecule (paramter)
            num_last_read = num_last_read -1 # We take the n-1 as the last
            dist_between_first_last = barcodes[barcode][num_last_read][1] - barcodes[barcode][num_first_read][0] # And we calculate the distance again

        #Stap 2 : Add the molecule (all reads in the molecule) in the dictionary dico_mol {"barcode-numero_of_molecule" : [reads]}
        barcode_obj = dico_Barcode[barcode.split('-@-')[1]]
        temp = [] #list of reads in the molecule
        for i in range(num_first_read,num_last_read+1) : # +1 because of python count between 0 and x (ex : 0 a 5 = 6 (0,1,2,3,4,5) )
            temp.append(barcodes[barcode][i]) 
            nb_reads_mol +=1

        dico_mol[f"{barcode}-%-{c}"] =temp
        barcode_obj.add_mol(temp)
        c +=1 
        if len(temp) <= 2 :
            not_mol +=1 #molecule number with less than 3 reads
            nb_read_not_count = nb_read_not_count + len(temp) # number of reads lost
         
        return [num_last_read, dico_mol, c, not_mol, nb_read_not_count, nb_reads_mol]
        
    #######------------------------------------------------------------------------------------########


    dico_mol ={}
    nb_read_per_mol = []
    not_mol = 0
    nb_reads_mol = 0
    

    for barcode in barcodes.keys() :

        c=0 #number of mol
        nb_read_not_count = 0 #number of reads not count per barcode
         
        result = []
        
        ### Intialisation of the first read ###
        num_first_read = 0
        
        num_last_read = len(barcodes[barcode])-1 # -1 because of python (start with 0)
        result = create_mol(num_first_read, num_last_read, mol_max_size, dico_mol, barcodes, barcode, c, not_mol, nb_read_not_count, dico_Barcode, nb_reads_mol)

        # Create molecule with all the reads which share the same barcode 
        while result[0] != len(barcodes[barcode])-1 : # while the last read isn't see (-1 because of python)
            num_first_read = result[0]+1
            dico_mol = result[1]
            c = result[2]
            not_mol = result[3]
            nb_read_not_count = result[4]
            nb_reads_mol = result[5]
                
            result = create_mol(num_first_read, num_last_read, mol_max_size, dico_mol, barcodes, barcode, c, not_mol, nb_read_not_count, dico_Barcode, nb_reads_mol)

        nb_read_not_count = result[4]
        c = result[2] # utile ?
        nb_reads_mol = result[5] # utile ?
        # Save information of reads not count
        barcode_obj = dico_Barcode[barcode.split('-@-')[1]]
        barcode_obj.count_not_count_reads(nb_read_not_count) #number of reads on molecule with less than 3 reads



    #####################################################################################################
    ### Calculate statistics
    #####################################################################################################
    # Number of reads per barcode before deconvolution (and number of read unmap)
    list_nb_read_bc = [] 
    unmap_list = []
    list_nb_read_bc_postD = []
    nb_read_lost = []
    nb_mol_per_bc = []
    all_mol = []
    nb_read_per_post_mol = []
    nb_read_post_filter =[]

    for bc, obj in dico_Barcode.items() :

        if obj.count_reads() == 0 :
            continue
        
        list_nb_read_bc.append(obj.count_reads())
        unmap_list.append(obj.unmap)
        nb_read_post_filter.append(obj.nb_read_post_deconvolution())
        nb_mol_per_bc.append(obj.count_mol()) 

    nb_unmap = sum(unmap_list)

    mol_size_filt = [] #molecule size post deconvolution (with only molecule with more than 3 reads)
    bc_more3R = {}
    bc_mol = {}
    nb_read_mol_more3R =[]
    nb_mol_filt = 0
    read_bc_mol = {}
    read_bc_molfilt = {}
    ## Recuperation des stats ##
    for mol in dico_mol.keys(): 
    # Nombre de reads par molécule
        nb_read_per_mol.append(len(dico_mol[mol]))
        size = dico_mol[mol][-1][1] - dico_mol[mol][0][0]

        # Filtering molecule with more than 3 reads
        if len(dico_mol[mol]) >2 :
            nb_mol_filt += 1
            nb_read_mol_more3R.append(len(dico_mol[mol]))
            mol_size_filt.append(size)

            barcode = (mol.split('-%-')[0]).split('@')[1] #On veut barcode pas par chrom

            if barcode not in bc_more3R :
                bc_more3R[barcode] = 1 #clé nombre de barcode, valeur nombre de mol par bc
                read_bc_molfilt[barcode] = len(dico_mol[mol])
            else : 
                bc_more3R[barcode] += 1
                read_bc_molfilt[barcode] += len(dico_mol[mol])
                
        if  len(dico_mol[mol]) > 0 :
            all_mol.append(size)
            barcode = (mol.split('-%-')[0]).split('@')[1]

            if barcode not in bc_mol :
                bc_mol[barcode] = 1
                read_bc_mol[barcode] = len(dico_mol[mol])
            else :
                bc_mol[barcode] += 1
                read_bc_mol[barcode] += len(dico_mol[mol])


    #######################
    ##### STATISTICS ######
    #######################

    # Nombre de reads
    # nb_read : nombre de read dans fichier bam
    # nb_reads_mol : nombre de reads en l'ensemble de molecules
    nb_read_filtmore3R = sum(nb_read_post_filter)

    #Nombre de barcodes
    nb_bc_more3R = len(bc_more3R)
    nb_bc_mol = len(bc_mol) 

    #Liste du nombre de mol par barcode
    list_nb_mol_bc_allmol = list(bc_mol.values())
    list_nb_mol_bc_filtmol = list(bc_more3R.values())

    #List du nombre de reads par bc
    nb_read_per_bc_mol = list(read_bc_mol.values())
    nb_read_per_bc_molfilt = list(read_bc_molfilt.values())

    ### Number of reads per barcode -----------------------------------------###    
    # No deconvolution
    mean_nb_read_bc_raw = statistics.mean(list_nb_read_bc)
    median_nb_read_bc_raw = statistics.median(list_nb_read_bc)
    # Deconvolution no filtering
    mean_nb_read_per_bc = statistics.mean(nb_read_per_bc_mol)
    median_nb_read_per_bc = statistics.median(nb_read_per_bc_mol)
    # Deconvoultion filtering
    mean_nb_read_per_bc_filt = statistics.mean(nb_read_per_bc_molfilt)
    median_nb_read_per_bc_filt = statistics.median(nb_read_per_bc_molfilt)
    ### ------------------------------------------------------------------- ###
    ### Number of molecules per barcode ------------------------------------###
    # Deconvoultion no filtering
    median_nb_mol_per_bc = statistics.median(list_nb_mol_bc_allmol)
    mean_nb_mol_per_bc = statistics.mean(list_nb_mol_bc_allmol)

    # Deconvoultion filtering
    mean_nb_mol_per_bc_filtmol =statistics.mean(list_nb_mol_bc_filtmol)
    median_nb_mol_per_bc_filtmol =statistics.median(list_nb_mol_bc_filtmol)
    ### ------------------------------------------------------------------- ###
    ### Number of read per mol ---------------------------------------------###
    # Deconvolution no filtering
    median_nb_read_per_mol =statistics.median(nb_read_per_mol)
    mean_nb_read_per_mol =statistics.mean(nb_read_per_mol)
    # Deconvolution filtering
    median_nb_read_per_mol_filt =statistics.median(nb_read_mol_more3R)
    mean_nb_read_per_mol_filt =statistics.mean(nb_read_mol_more3R)
    ### ------------------------------------------------------------------- ###
    ### Molecule size ----------------------------------------------------- ###
    # Deconvolution no filtering
    median_all_mol_size = statistics.median(all_mol)
    mean_all_mol_size = statistics.mean(all_mol)
    # Deconvolution filtering
    median_mol_size_filt = statistics.median(mol_size_filt)
    mean_mol_size_filt = statistics.mean(mol_size_filt)
    ### ------------------------------------------------------------------- ###


    #####################################################################################################
    #####################################################################################################
    ### Coverage et deepth
    moy_reads_mol_covs = (mean_nb_read_per_mol * read_size) / mean_all_mol_size
    moy_reads_molfilt_covs = (mean_nb_read_per_mol_filt *read_size)/mean_mol_size_filt

    coverage_mol_genome = ( len(dico_mol) * mean_all_mol_size) / genome_size
    coverage_molfilt_genome = ( nb_mol_filt* mean_mol_size_filt) / genome_size
        
    # 3. Génome coverage by reads
    coverage_read_genome = (nb_read * read_size) / genome_size #4641652 

    cov_read_genome_mol =(nb_reads_mol * read_size) / genome_size
    cov_read_genome_molfilt =(nb_read_filtmore3R * read_size) / genome_size


############################################################
    #OUTPUT
    table = PrettyTable()
    table.field_names = [" ","Raw data", "Creating molecules no filtering","Creating molecules with more than 3 reads"]
    table.add_row(["Number of reads", nb_read , nb_reads_mol ,nb_read_filtmore3R])
    table.add_row(["Number of barcodes", len(dico_Barcode), nb_bc_mol, nb_bc_more3R ])
    table.add_row(["Number of unmapped reads",nb_unmap , "","" ])
    table.add_row(["Mean number of reads per barcode", round(mean_nb_read_bc_raw,1), round(mean_nb_read_per_bc,1),round(mean_nb_read_per_bc_filt,1) ])
    table.add_row(["Median number of reads per barcode",round(median_nb_read_bc_raw,1) ,round(median_nb_read_per_bc,1) , round(median_nb_read_per_bc_filt,1)])
    table.add_row(["Read depth per genome",round(coverage_read_genome,1) ,round(cov_read_genome_mol,1) ,round(cov_read_genome_molfilt,1) ])
    table.add_row(["---","---" ,"---" , "---"])
    table.add_row(["Number of molecules"," ", len(dico_mol),nb_mol_filt])
    table.add_row(["Mean molecule size","" ,round(mean_all_mol_size,1) ,round(mean_mol_size_filt,1) ])
    table.add_row(["Median molecule size","" , round(median_all_mol_size,1), round(median_mol_size_filt,1)])
    table.add_row(["Mean number of reads per molecule","" ,round(mean_nb_read_per_mol,1) ,round(mean_nb_read_per_mol_filt,1) ])   
    table.add_row(["Median number of reads per molecule", "",round(median_nb_read_per_mol,1) ,round(median_nb_read_per_mol_filt,1) ])
    table.add_row(["Mean number of molecules per barcode","" ,round(mean_nb_mol_per_bc,1) , round(mean_nb_mol_per_bc_filtmol,1)])
    table.add_row(["Median number of molecules per barcode",""  ,round(median_nb_mol_per_bc,1) , round(median_nb_mol_per_bc_filtmol,1)])
    table.add_row(["Read depth per molecule","" ,round(moy_reads_mol_covs,1) ,round(moy_reads_molfilt_covs,1) ])
    table.add_row(["Molecule depth per genome", "",round(coverage_mol_genome,1) ,round(coverage_molfilt_genome,1) ])
    # Afficher le tableau
    print(table)


    #### OUTPUT GRAPHE
    plt.style.use('seaborn-v0_8-paper')

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Molecule size distribution
    axes[0, 0].hist(mol_size_filt, bins=50, linewidth=0.5, edgecolor="white", color="teal")
    axes[0, 0].axvline(mean_mol_size_filt, color='black', linestyle='dotted', linewidth=1, label=f'Mean: {mean_mol_size_filt:.2f}')
    axes[0, 0].axvline(median_mol_size_filt, color='black', linestyle='solid', linewidth=1, label=f'Median: {median_mol_size_filt:.2f}')
    axes[0, 0].legend()
    axes[0, 0].set_title("Histogram of molecule size")
    axes[0, 0].set_xlabel("Size (pb)")
    axes[0, 0].set_ylabel("Number of molecules")

    # Number of reads per molecule vs molecule size
    axes[0, 1].scatter(mol_size_filt, nb_read_mol_more3R, color='olivedrab', s=4, alpha = 0.1 )
    axes[0, 1].set_title("Number of reads per molecule VS molecule size")
    axes[0, 1].set_xlabel("Size (pb)")
    axes[0, 1].set_ylabel("Number of reads")

    # Histogram of number of reads per molecule
    axes[1, 0].hist(nb_read_mol_more3R, bins=50, linewidth=0.5, edgecolor="white", color="coral")
    axes[1, 0].axvline(mean_nb_read_per_mol, color='black', linestyle='dotted', linewidth=1, label=f'Mean: {mean_nb_read_per_mol_filt:.2f}')
    axes[1, 0].axvline(median_nb_read_per_mol, color='black', linestyle='solid', linewidth=1, label=f'Median: {median_nb_read_per_mol_filt:.2f}')
    axes[1, 0].legend()
    axes[1, 0].set_title("Histogram of number of reads per molecule")
    axes[1, 0].set_xlabel("Number of reads")
    axes[1, 0].set_ylabel("Number of molecules")

    # Histogramme of number of reads per barcode
    axes[1, 1].hist(nb_read_per_bc_molfilt, bins=50, linewidth=0.5, edgecolor="white", color="goldenrod")
    axes[1, 1].axvline(mean_nb_read_per_bc_filt, color='black', linestyle='dotted', linewidth=1, label=f'Mean: {mean_nb_read_per_bc_filt:.2f}')
    axes[1, 1].axvline(median_nb_read_per_bc_filt, color='black', linestyle='solid', linewidth=1, label=f'Median: {median_nb_read_per_bc_filt:.2f}')
    axes[1, 1].legend()
    axes[1, 1].set_title("Histogram of number of reads per barcode post-deconvolution")
    axes[1, 1].set_xlabel("Number of reads")
    axes[1, 1].set_ylabel("Number of barcodes")



    plt.tight_layout()
    plt.show()
    plt.savefig(output_histo)

    print(f"The graphics have been saved as '{output_histo}'", file = sys.stderr)
    data = [
        ["Number of reads", nb_read, nb_reads_mol, nb_read_filtmore3R],
        ["Number of barcodes", len(dico_Barcode), nb_bc_mol, nb_bc_more3R],
        ["Number of unmapped reads", nb_unmap, "", ""],
        ["Mean number of reads per barcode",
        round(mean_nb_read_bc_raw, 1),
        round(mean_nb_read_per_bc, 1),
        round(mean_nb_read_per_bc_filt, 1)],
        ["Median number of reads per barcode",
        round(median_nb_read_bc_raw, 1),
        round(median_nb_read_per_bc, 1),
        round(median_nb_read_per_bc_filt, 1)],
        ["Read depth per genome",
        round(coverage_read_genome, 1),
        round(cov_read_genome_mol, 1),
        round(cov_read_genome_molfilt, 1)],
        ["---", "---", "---", "---"],
        ["Number of molecules", "", len(dico_mol), nb_mol_filt],
        ["Mean molecule size", "", round(mean_all_mol_size, 1), round(mean_mol_size_filt, 1)],
        ["Median molecule size", "", round(median_all_mol_size, 1), round(median_mol_size_filt, 1)],
        ["Mean number of reads per molecule", "",
        round(mean_nb_read_per_mol, 1),
        round(mean_nb_read_per_mol_filt, 1)],
        ["Median number of reads per molecule", "",
        round(median_nb_read_per_mol, 1),
        round(median_nb_read_per_mol_filt, 1)],
        ["Mean number of molecules per barcode", "",
        round(mean_nb_mol_per_bc, 1),
        round(mean_nb_mol_per_bc_filtmol, 1)],
        ["Median number of molecules per barcode", "",
        round(median_nb_mol_per_bc, 1),
        round(median_nb_mol_per_bc_filtmol, 1)],
        ["Read depth per molecule", "",
        round(moy_reads_mol_covs, 1),
        round(moy_reads_molfilt_covs, 1)],
        ["Molecule depth per genome", "",
        round(coverage_mol_genome, 1),
        round(coverage_molfilt_genome, 1)],
    ]

    # Création du DataFrame
    df = pd.DataFrame(
        data,
        columns=["Metric", "Raw data", "Creating molecules no filtering", "Creating molecules with more than 3 reads"]
    )

    # Sauvegarde en CSV
    df.to_csv(out_table, index=False)
    print(f"The table have been saved as '{out_table}'", file = sys.stderr)
