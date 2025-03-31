#!/usr/bin/env python3

#*****************************************************************************
#  Name: SVJedi-Tag
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
import argparse
import subprocess
import statistics
#output
from prettytable import PrettyTable
import  matplotlib.pyplot as plt
import pandas as pd

# input = sys.argv[1]
# mol_size = int(sys.argv[2])
# #mol_size = 100000

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
    
    def count_not_count_reads(self, nb_reads):
        self.not_count_reads += nb_reads

    def nb_read_post_deconvolution(self):
        return (len(self.reads) - self.not_count_reads)
    
#####################################################################################################


def main(args) :
    """ Main method """

    #####################################################################################################
    #Arguments
    #####################################################################################################
    parser = argparse.ArgumentParser()

    parser.add_argument( "-b", "--bam", metavar="<sort_bam_file>", help= "Bam file out of a mapping linked-reads/reference and with BX tag", type=str, required=True)
    parser.add_argument( "-s", "--molecule_max_size", metavar="<molecule_max_size>", help="Maximum size between two reads to share a same barcode and from the same molecule [default=100000]",type=int, required=False, default=100000)
    parser.add_argument( "-G", "--graph_output", metavar="<graphe_output_path/name_file>", help="Path/Name_file for the output graph.png [default=ACLR_graph.png]", type=str, required=False, default="ACLR_graph.png")
    parser.add_argument( "-o", "--output_table", metavar="<output_table_path/name_file>",help="Path/Name_file for the output table.csv [default=ACLR_table.csv] ", type=str, required=False, default="ACLR_table.csv")
    parser.add_argument( "-g", "--genome_size", metavar="<genome_size>", help="Genome size required to calculate depth", type=int, required=False, default=0)
    parser.add_argument( "-r", "--read_size", metavar="<read_size>", help="Read size required to calculate depth", type=int, required=False, default=150)


    args = parser.parse_args()
    input = args.bam
    mol_size = args.molecule_max_size 
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

    ### Function ###  
    def create_mol(num_first_read, num_last_read, mol_size, dico_mol, barcodes, barcode, c, not_mol, nb_read_not_count, dico_Barcode):

        dist_between_first_last = barcodes[barcode][num_last_read][1] - barcodes[barcode][num_first_read][0]  # La distance incluant les reads donc debut du premier - fin du dernier
        #Step 1 : Trouver les extrémités de la molécule 
        while dist_between_first_last > mol_size : #Si la distance entre le 1er read et le (n-x) est inférieur a molécule size alors ça s'arrete
            num_last_read = num_last_read -1 # On prend le n-1 comme dernier
            dist_between_first_last = barcodes[barcode][num_last_read][1] - barcodes[barcode][num_first_read][0]

        # Ajouter la molécule au dico mol (tous les readss compris entre les extrémités)
        #dico_mol[f"{barcode}-{c}"]=[] # Initialise clé

        barcode_obj = dico_Barcode[barcode.split('-@-')[1]]
        temp = []
        for i in range(num_first_read,num_last_read+1) :
            temp.append(barcodes[barcode][i]) # +1 car python et donc je dois compter de 0 a X (ex : 0 a 5 = 6 (0,1,2,3,4,5) )
            #dico_mol[f"{barcode}-{c}"].append(barcodes[barcode][i])
        if len(temp) > 2 :  #Ne pas prend en compte les molecules de moins de 3 reads
            dico_mol[f"{barcode}-{c}"] =temp
            barcode_obj.add_mol(temp)
            c +=1 # molecule number with more than 3 reads
        else :
            not_mol +=1 #molecule number woth less than 3 reads
            nb_read_not_count = nb_read_not_count + len(temp) # number of reads lost
         
        return [num_last_read, dico_mol, c, not_mol, nb_read_not_count]
        
    ###############


    dico_mol ={}
    #nb_mol_per_bc = []
    nb_read_per_mol = []
    post_mol_size = []
    not_mol=0


    for barcode in barcodes.keys() :

        c=0
        nb_read_not_count = 0
        result = []

        # Intialisation du prrmier read
        num_first_read = 0
        num_last_read = len(barcodes[barcode])-1 # -1 car python commence a 0
        result = create_mol(num_first_read, num_last_read, mol_size, dico_mol, barcodes, barcode, c, not_mol, nb_read_not_count, dico_Barcode)

        # Run create_mol sur l'ensemble du barcode temps que une molécule ne finit pas par le dernier read. 
        while result[0] != len(barcodes[barcode])-1 : #-1 car len compte pas a partir de 0
            num_first_read = result[0]+1
            dico_mol = result[1]
            c = result[2]
            not_mol = result[3]
            nb_read_not_count = result[4]
                
            result = create_mol(num_first_read, num_last_read, mol_size, dico_mol, barcodes, barcode, c, not_mol, nb_read_not_count, dico_Barcode)

        nb_read_not_count = result[4]
        c = result[2]
        # Save information of reads not count
        barcode_obj = dico_Barcode[barcode.split('-@-')[1]]
        barcode_obj.count_not_count_reads(nb_read_not_count)
        #nb_mol_per_bc.append(c) # liste des nombres de molécule par barcodes


    #####################################################################################################
    ### Calculate statistics
    #####################################################################################################
    # Number of reads per barcode before deconvolution (and number of read unmap)
    list_nb_read_bc = [] 
    unmap_list = []
    list_nb_read_bc_postD = []
    nb_read_lost = []
    nb_mol_per_bc = []

    for bc, obj in dico_Barcode.items() :
        list_nb_read_bc.append((obj.count_reads()))
        unmap_list.append(obj.unmap)
        list_nb_read_bc_postD.append(obj.nb_read_post_deconvolution())
        nb_read_lost.append(obj.not_count_reads)
        nb_mol_per_bc.append(obj.count_mol())

    # Supprimer les zeros pour ne prendre en compte que les barcodes qui n'ont pas de molecules
    list_nb_read_bc = [x for x in list_nb_read_bc if x > 0] 
    list_nb_read_bc_postD = [x for x in list_nb_read_bc_postD if x > 0]
    nb_mol_per_bc= [x for x in nb_mol_per_bc if x > 0]

    mean_nb_read_per_bc = statistics.mean(list_nb_read_bc)
    median_nb_read_per_bc = statistics.median(list_nb_read_bc)

    mean_nb_read_per_bc_postD = statistics.mean(list_nb_read_bc_postD)
    median_nb_read_per_bc_postD = statistics.median(list_nb_read_bc_postD)

    nb_unmap = sum(unmap_list)
    nb_reads_lost = sum(nb_read_lost)

    ## Recuperation des stats ##

    for mol in dico_mol.keys(): 
    # Nombre de reads par molécule
        nb_read_per_mol.append(len(dico_mol[mol]))
    # Taille des molécules
        size = dico_mol[mol][-1][1] - dico_mol[mol][0][0]
        post_mol_size.append(size) 

        
    median_nb_mol_per_bc = statistics.median(nb_mol_per_bc)
    mean_nb_mol_per_bc = statistics.mean(nb_mol_per_bc)


    median_nb_read_per_mol =statistics.median(nb_read_per_mol)
    mean_nb_read_per_mol =statistics.mean(nb_read_per_mol)

    median_mol_size = statistics.median(post_mol_size)
    mean_mol_size = statistics.mean(post_mol_size)


    #####################################################################################################
    #####################################################################################################
    ### Coverage et deepth
    moy_reads_mol_covs = (mean_nb_read_per_mol * read_size) / mean_mol_size

    if genome_size == 0 : 
        coverage_mol_genome = " "
        coverage_read_genome = " "
    else :
        #coverage_mol_genome = (10000 * len(molecule)) / genome_size
        coverage_mol_genome = ( len(dico_mol) * mean_mol_size) / genome_size
        
        # 3. Génome coverage by reads

        coverage_read_genome = (nb_read * read_size) / genome_size #4641652 


    per_read_lost = int((nb_reads_lost /nb_read) * 100)
    per_bc_lost = int((len(list_nb_read_bc_postD)/len(dico_Barcode))*100)
    #############################################################################################
    #OUTPUT
    table = PrettyTable()
    table.field_names = ["Pre-deconvolution", " "]

    table.add_row(["Number of reads", nb_read])
    table.add_row(["Number of unmapped reads", nb_unmap])
    table.add_row(["Number of barcodes", len(dico_Barcode)])
    table.add_row(["Read depth per genome", coverage_read_genome])
    table.add_row(["Mean number of reads per barcode (pre-deconvolution)", mean_nb_read_per_bc])
    table.add_row(["Median number of reads per barcode (pre-deconvolution)", median_nb_read_per_bc])
    table.add_row(["--------------------------------------------------------", "------------------------"])
    table.add_row(["Post-deconvolution", " "])
    table.add_row(["--------------------------------------------------------", "------------------------"])
    table.add_row(["Number of barcodes (post-deconvolution)", f"{len(list_nb_read_bc_postD)} ({per_bc_lost}% of barcodes lost)"])
    table.add_row(["Total number of molecules", len(dico_mol)])
    table.add_row(["Number of molecules with fewer than 3 reads", not_mol])
    table.add_row(["Number of reads lost", f"{nb_reads_lost} ({per_read_lost}% of reads)"])
    table.add_row(["Number of reads (post-deconvolution)", (nb_read - nb_reads_lost)])
    table.add_row(["  ", "  "])
    table.add_row(["Mean molecule size", mean_mol_size])
    table.add_row(["Median molecule size", median_mol_size])
    table.add_row(["  ", "  "])
    table.add_row(["Mean number of reads per barcode (post-deconvolution)", mean_nb_read_per_bc_postD])
    table.add_row(["Median number of reads per barcode (post-deconvolution)", median_nb_read_per_bc_postD])
    table.add_row(["  ", "  "])
    table.add_row(["Mean number of reads per molecule", mean_nb_read_per_mol])
    table.add_row(["Median number of reads per molecule", median_nb_read_per_mol])
    table.add_row(["  ", "  "])
    table.add_row(["Mean number of molecules per barcode", mean_nb_mol_per_bc])
    table.add_row(["Median number of molecules per barcode", median_nb_mol_per_bc])
    table.add_row(["  ", "  "])
    table.add_row(["Read depth per molecule", moy_reads_mol_covs])
    table.add_row(["Molecule depth per genome", coverage_mol_genome])


    # Afficher le tableau
    print(table)


    #### OUTPUT GRAPHE

    plt.style.use('seaborn-v0_8-paper')

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Molecule size distribution
    axes[0, 0].hist(post_mol_size, bins=50, linewidth=0.5, edgecolor="white", color="teal")
    axes[0, 0].axvline(mean_mol_size, color='black', linestyle='dotted', linewidth=1, label=f'Mean: {mean_mol_size:.2f}')
    axes[0, 0].axvline(median_mol_size, color='black', linestyle='solid', linewidth=1, label=f'Median: {median_mol_size:.2f}')
    axes[0, 0].legend()
    axes[0, 0].set_title("Histogram of molecule size")
    axes[0, 0].set_xlabel("Size (pb)")
    axes[0, 0].set_ylabel("Number of molecules")

    # Number of reads per molecule vs molecule size
    axes[0, 1].scatter(post_mol_size, nb_read_per_mol, color='olivedrab', s=10, alpha = 0.4 )
    axes[0, 1].set_title("Number of reads per molecule VS molecule size")
    axes[0, 1].set_xlabel("Size (pb)")
    axes[0, 1].set_ylabel("Number of reads")

    # Histogram of number of reads per molecule
    axes[1, 0].hist(nb_read_per_mol, bins=50, linewidth=0.5, edgecolor="white", color="coral")
    axes[1, 0].axvline(mean_nb_read_per_mol, color='black', linestyle='dotted', linewidth=1, label=f'Mean: {mean_nb_read_per_mol:.2f}')
    axes[1, 0].axvline(median_nb_read_per_mol, color='black', linestyle='solid', linewidth=1, label=f'Median: {median_nb_read_per_mol:.2f}')
    axes[1, 0].legend()
    axes[1, 0].set_title("Histogram of number of reads per molecule")
    axes[1, 0].set_xlabel("Number of reads")
    axes[1, 0].set_ylabel("Number of molecules")

    # # Quatrième graphique : Histogramme du nombre de reads par barcode
    # axes[1, 1].hist(list_nb_read_bc, bins=50, linewidth=0.5, edgecolor="white", color="goldenrod")
    # axes[1, 1].axvline(mean_nb_read_per_bc, color='black', linestyle='dotted', linewidth=1, label=f'Mean: {mean_nb_read_per_bc:.2f}')
    # axes[1, 1].axvline(median_nb_read_per_bc, color='black', linestyle='solid', linewidth=1, label=f'Median: {median_nb_read_per_bc:.2f}')
    # axes[1, 1].legend()
    # axes[1, 1].set_title("Histogram of number of reads per barcode preD")
    # axes[1, 1].set_xlabel("Number of reads")
    # axes[1, 1].set_ylabel("Number of barcodes")

    # # Quatrième graphique : Histogramme du nombre de reads par barcode
    axes[1, 1].hist(list_nb_read_bc_postD, bins=50, linewidth=0.5, edgecolor="white", color="goldenrod")
    axes[1, 1].axvline(mean_nb_read_per_bc_postD, color='black', linestyle='dotted', linewidth=1, label=f'Mean: {mean_nb_read_per_bc_postD:.2f}')
    axes[1, 1].axvline(median_nb_read_per_bc_postD, color='black', linestyle='solid', linewidth=1, label=f'Median: {median_nb_read_per_bc_postD:.2f}')
    axes[1, 1].legend()
    axes[1, 1].set_title("Histogram of number of reads per barcode post-deconvolution")
    axes[1, 1].set_xlabel("Number of reads")
    axes[1, 1].set_ylabel("Number of barcodes")



    plt.tight_layout()
    plt.show()
    plt.savefig(output_histo)

    print(f"The graphics have been saved as '{output_histo}'")


### OUTPUT TABLE 

    # Créer les données sous forme de dictionnaire
    data = {
        "Parameters": [
            "Number of reads",
            "Number of reads unmapped ",
            "Number of barcodes",
            "Read depth by genome",
            "Median number of reads per barcode pre-deconvolution",
            "Mean number of reads per barcode pre-deconvolution",
            "Number of barcodes post deconvolution",
            "Number of molecules total",
            "Number of molecules with less 3 reads",
            "Number of reads not count",
            "Number of reads post deconvolution",
            "Mean of molecule size",
            "Median of molecule size",
            "Mean number of reads per barcode post-deconvolution",
            "Median number of reads per barcode post-deconvolution",
            "Mean number of reads per molecule",
            "Median number of reads per molecule",
            "Mean number of molecule per barcode",
            "Median number of molecule per barcode",
            "Read depth by molecule",
            "Molecule depth by genome"
        ],
        "Values": [
            nb_read,
            nb_unmap,
            len(dico_Barcode),
            coverage_read_genome,
            mean_nb_read_per_bc,
            median_nb_read_per_bc,
            f"{len(list_nb_read_bc_postD)} ({per_bc_lost}% of barcodes lost)",
            len(dico_mol),
            not_mol,
            f"{nb_reads_lost} ({per_read_lost}% of reads)",
            (nb_read - nb_reads_lost),
            mean_mol_size,
            median_mol_size,
            mean_nb_read_per_bc_postD,
            median_nb_read_per_bc_postD,
            mean_nb_read_per_mol,
            median_nb_read_per_mol,
            mean_nb_mol_per_bc,
            median_nb_mol_per_bc,
            moy_reads_mol_covs,
            coverage_mol_genome
        ]
    }

    # Créer un DataFrame pandas
    df = pd.DataFrame(data)

    # Sauvegarder le DataFrame au format CSV
    df.to_csv(out_table, index=False)

    # Afficher le DataFrame pour vérifier
    
    print(f"The table have been saved as '{out_table}'")



# Run function main
if __name__ == "__main__":
    if sys.argv == 1:
        sys.exit("Error: missing arguments")

    else:
        main(sys.argv[1:])


