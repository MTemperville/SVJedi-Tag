#!/usr/bin/env python3

import sys
import re
import random
import statistics 
import matplotlib.pyplot as plt
import math
from decimal import *


def add_subparser(subparsers):
    p = subparsers.add_parser("likelihood-calibration", help="calibrate genotype likelihoods")
    p.add_argument( "--gaf", metavar="<alignment file>", help="", type=str, required=True)
    p.add_argument( "--gfa", metavar="<variation graph file>", help="", type=str, required=True)
    p.add_argument( "-w", metavar="<Windows size>", help="", type=int, required=True)
    p.add_argument( "-n", metavar="<Number of windows>", help="number of interest region/error rate", type=int, default=100)
    p.add_argument( "--regionssizesIn", help="regions sizes for region in window", type=int, default=10000)
    p.add_argument( "--regionssizesOut", help="regions sizes for region out window", type=int, default = 10000)
    p.add_argument( "-p", "--plot", metavar="path/plot.png", type=str, default = "Histogramme_error")
    p.set_defaults(func=main, _parser = p)

def main(args):
    if len(sys.argv) == 2:   # no args
        args._parser.print_help()
        sys.exit(1)

    inputGAF = args.gaf
    inputGFA = args.gfa
    windows_nb = args.n
    windows_Size = args.w  #Taille qui sépare les deux bk
    size_ro = args.regionssizesIn
    size_ri = args.regionssizesOut
    plot = args.plot 

    ### --------- Functions --------------- ###
    def create_region(bk1, bk2, size_ro, size_ri, windows_Size ):
        if size_ri > (windows_Size/2) :
            Coord_R2 = (bk1, (bk1 + (windows_Size/2)))
            Coord_R3 = ((bk2 - (windows_Size/2)), bk2)
        else :        
            Coord_R2 = (bk1, (bk1 + size_ri))
            Coord_R3 = ((bk2 - size_ri), bk2)
        Coord_R4 = (bk2, (bk2 + size_ro))
        Coord_R1 = ((bk1 - size_ro), bk1)

        return Coord_R1, Coord_R2, Coord_R3, Coord_R4
    

    def plot_creation(error_list, plot,windows_Size,windows_nb,size_ro, size_ri, mean_error,std_error,mean_diff,median_error):
        #sns.kdeplot(error_list, fill=True, color="blue")

        plt.scatter(error_list, range(len(error_list)),color='lightblue', s=50)

        #Axes
        # plt.xticks(np.arange(0.0, 1.1, 0.1), fontsize=5, rotation=0)
        # plt.yticks(fontsize=5, rotation=0)
        # plt.xlim(0, 1)

        ## Title
        plt.text(0.5, 1.08, "Density Curve of Error Rates", 
                ha="center", va="bottom", fontsize=12, fontweight="bold", transform=plt.gca().transAxes)
        plt.text(0.5, 1.04, f"For {windows_nb} Windows of Size {windows_Size}pb", 
                ha="center", va="bottom", fontsize=10, transform=plt.gca().transAxes)
        plt.text(0.5, 1.00, f"(Region sizes: {size_ro}pb out windows and {size_ri}pb in)", 
                ha="center", va="bottom", fontsize=8, style="italic", transform=plt.gca().transAxes)


        plt.xlabel("Error rate", fontsize=7)
        plt.ylabel("Windows", fontsize=7)
        plt.grid(True)

        plt.legend([f"Error rate density"], fontsize=6)
        textstr = f"Mean = {mean_error:.3f}\nStd = {std_error:.3f}\nMedian={median_error:.3f}\nLikeDiff = {mean_diff:.3f}"
        
        plt.text(0.7, 0.8, textstr, fontsize=7,
         transform=plt.gca().transAxes,  # coordonnées relatives (0 à 1)
         bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.3"))

        # sns.histplot(error_list, kde=True, color="skyblue")
        # plt.title("Distribution et densité des taux d'erreur")
        # plt.xlabel("Taux d'erreur")
        # plt.ylabel("Fréquence")

        plt.savefig(f"{plot}.png", dpi=300, bbox_inches='tight')


    def likelihood_diffBests(c1, c2, e): # calcule genotype likelihood pour les 3 genotypes et conclu
     
        lik0 = -(Decimal(c1*math.log10(1-e)) + Decimal(c2*math.log10(e)) ) # c1xlog10x(1−e))+ c2xlog10x(e)   (log10 pour facilité calcule)
        lik1 = -(Decimal((c1+c2)*math.log10(1/2))) #(c1+c2)xlog10(1/2)
        lik2 = -(Decimal(c2*math.log10(1-e)) + Decimal(c1*math.log10(e))) # c2xlog10x(1−e))+ c1xlog10x(e)
        # lik0 = 0/0 ; lik1 = 0/1 ; lik2 = 1/1
        L = [lik0, lik1, lik2]
        L.sort()
        diff = L[1] - L[0]

        return diff


    ### ------- Run --------- ### 
     
    # 1. Take nodes from GFA
    all_nodes = []
    with open(inputGFA, 'r') as gfa : 
        for line in gfa :
            if line.startswith('S'):
                node = line.split('\t')[1]
                chr = node.split(':')[0]
                if not chr.startswith('scaffold'):
                    all_nodes.append(node)
    #print(all_nodes)
    print("#Take node from GFA : Done", file = sys.stderr)
    #'''
    # 2. Create windows and regions
    #avoid_pos = [["LG1",8342182,33487673],["LG2",14083320,20869940],["LG3",7486933,13829649],["LG4",1088816,7995568], ["LG4",22421881,25145365],["LG4",30622035,31991919], ["LG5",15940464,32665323] ] 
    i = 0
    dico_windows = {}
    
    list_w_chr = []
    list_w_node = []
    
    while i != windows_nb :
        #Create window
        node = random.choice(all_nodes)
        if  (int(node.split(":")[1]) + size_ro) < (int(node.split(":")[2]) - size_ro):
            bk1 = random.randint((int(node.split(":")[1]) + size_ro), (int(node.split(":")[2]) - size_ro))
            bk2 = bk1 + windows_Size

            if bk2 > int(node.split(":")[2]) : 
                continue

            # for chr, start, end in avoid_pos:
            #     if node.split(':')[0] == chr and start <= bk1 <= end or node.split(':')[0] == chr and start <= bk2 <= end : 
            #         continue
            
            i += 1

            list_w_chr.append((f"{node.split(':')[0]}", bk1, bk2))
            list_w_node.append((node,bk1,bk2))

            #Creation region
            cR1, cR2, cR3, cR4 = create_region(bk1, bk2, size_ro, size_ri, windows_Size)
            if node not in dico_windows :
                #dico_windows[node]=[{"bk" : (bk1,bk2), "regions":[cR1, cR2, cR3, cR4],"R1":set(),"R2":set(),"R3":set(),"R4":set()}]
                dico_windows[node]=[{"bk" : (bk1,bk2), "regions":[cR1, cR2, cR3, cR4],"R1":{},"R2":{},"R3":{},"R4":{}}]
            else :
                #dico_windows[node].append({"bk" : (bk1,bk2), "regions":[cR1, cR2, cR3, cR4], "R1":set(),"R2":set(),"R3":set(),"R4":set()})
                dico_windows[node].append({"bk" : (bk1,bk2), "regions":[cR1, cR2, cR3, cR4], "R1":{},"R2":{},"R3":{},"R4":{}})
    print("#Create windows and regions : Done", file = sys.stderr)

    #print(f"CW\t{list_w_chr}")
    #print(f"NW\t{list_w_node}")
    #'''


    '''
    # Fenetre Fixe
    dico_windows = {}
    windows_list =  [('Chr1:2380175:2480173', 2413842, 2463842), ('Chr4:1:36243598', 31317396, 31367396), ('Chr4:1:36243598', 24203332, 24253332), ('Chr1:4452824:4552822', 4530802, 4580802), ('Chr3:1:38949715', 5644220, 5694220), ('Chr1:2380175:2480173', 2410799, 2460799), ('Chr2:1:38097547', 28604688, 28654688), ('Chr1:4383182:4452823', 4435938, 4485938), ('Chr1:2480174:4283182', 3302055, 3352055), ('Chr1:4283183:4383181', 4313058, 4363058), ('Chr4:1:36243598', 20573709, 20623709), ('Chr5:1:51160880', 5593000, 5643000), ('Chr4:1:36243598', 4408582, 4458582), ('Chr1:2380175:2480173', 2457650, 2507650), ('Chr4:1:36243598', 6996354, 7046354), ('Chr5:1:51160880', 16207849, 16257849), ('Chr1:1:2380174', 2175051, 2225051), ('Chr2:1:38097547', 29683274, 29733274), ('Chr4:1:36243598', 23907318, 23957318), ('Chr3:1:38949715', 2846822, 2896822), ('Chr1:4283183:4383181', 4365944, 4415944), ('Chr1:4383182:4452823', 4403363, 4453363), ('Chr6:1:5105880', 1698346, 1748346), ('Chr1:4283183:4383181', 4353826, 4403826), ('Chr1:4383182:4452823', 4436179, 4486179), ('Chr4:1:36243598', 27822055, 27872055), ('Chr3:1:38949715', 14490650, 14540650), ('Chr5:1:51160880', 38696006, 38746006), ('Chr6:1:5105880', 4758214, 4808214), ('Chr1:4283183:4383181', 4346735, 4396735), ('Chr1:2380175:2480173', 2399103, 2449103), ('Chr1:4383182:4452823', 4396203, 4446203), ('Chr6:1:5105880', 3698740, 3748740), ('Chr1:4552823:48873184', 31618484, 31668484), ('Chr1:2380175:2480173', 2420231, 2470231), ('Chr2:1:38097547', 15438710, 15488710), ('Chr1:4452824:4552822', 4542026, 4592026), ('Chr3:1:38949715', 18068347, 18118347), ('Chr4:1:36243598', 30561876, 30611876), ('Chr1:4552823:48873184', 7354370, 7404370), ('Chr1:2380175:2480173', 2394857, 2444857), ('Chr1:4552823:48873184', 9404683, 9454683), ('Chr1:4552823:48873184', 18327213, 18377213), ('Chr1:2380175:2480173', 2399029, 2449029), ('Chr1:4383182:4452823', 4404863, 4454863), ('Chr1:2380175:2480173', 2421937, 2471937), ('Chr1:4383182:4452823', 4404700, 4454700), ('Chr1:4283183:4383181', 4335727, 4385727), ('Chr1:4552823:48873184', 18061072, 18111072), ('Chr1:4452824:4552822', 4489977, 4539977), ('Chr5:1:51160880', 11401361, 11451361), ('Chr1:1:2380174', 1205720, 1255720), ('Chr4:1:36243598', 5138549, 5188549), ('Chr2:1:38097547', 8157189, 8207189), ('Chr1:4383182:4452823', 4409829, 4459829), ('Chr2:1:38097547', 31604119, 31654119), ('Chr1:4283183:4383181', 4336508, 4386508), ('Chr1:4383182:4452823', 4409705, 4459705), ('Chr1:4383182:4452823', 4442381, 4492381), ('Chr1:4283183:4383181', 4357352, 4407352), ('Chr1:4452824:4552822', 4487772, 4537772), ('Chr1:4452824:4552822', 4469757, 4519757), ('Chr1:1:2380174', 1002784, 1052784), ('Chr1:2480174:4283182', 2707357, 2757357), ('Chr1:4452824:4552822', 4491722, 4541722), ('Chr3:1:38949715', 15322672, 15372672), ('Chr1:4383182:4452823', 4411719, 4461719), ('Chr1:4452824:4552822', 4468614, 4518614), ('Chr3:1:38949715', 25722739, 25772739), ('Chr1:4452824:4552822', 4530167, 4580167), ('Chr1:4552823:48873184', 13238444, 13288444), ('Chr3:1:38949715', 30095459, 30145459), ('Chr1:4283183:4383181', 4305839, 4355839), ('Chr2:1:38097547', 34895703, 34945703), ('Chr1:4452824:4552822', 4507162, 4557162), ('Chr6:1:5105880', 3897149, 3947149), ('Chr1:1:2380174', 2291076, 2341076), ('Chr6:1:5105880', 655496, 705496), ('Chr1:4283183:4383181', 4347356, 4397356), ('Chr1:4452824:4552822', 4486433, 4536433), ('Chr5:1:51160880', 49118235, 49168235), ('Chr1:2380175:2480173', 2435457, 2485457), ('Chr1:4283183:4383181', 4337763, 4387763), ('Chr1:4452824:4552822', 4513661, 4563661), ('Chr1:1:2380174', 1075222, 1125222), ('Chr1:4283183:4383181', 4362590, 4412590), ('Chr1:2380175:2480173', 2442480, 2492480), ('Chr1:4552823:48873184', 10877797, 10927797), ('Chr6:1:5105880', 3497175, 3547175), ('Chr1:2380175:2480173', 2403408, 2453408), ('Chr1:4283183:4383181', 4316928, 4366928), ('Chr1:4283183:4383181', 4331773, 4381773), ('Chr1:4552823:48873184', 5678581, 5728581), ('Chr1:2480174:4283182', 3906792, 3956792), ('Chr1:1:2380174', 208291, 258291), ('Chr4:1:36243598', 895253, 945253), ('Chr3:1:38949715', 33269826, 33319826), ('Chr1:4452824:4552822', 4542382, 4592382), ('Chr4:1:36243598', 18623984, 18673984), ('Chr1:1:2380174', 856761, 906761)]

    for windows in windows_list :
        cR1, cR2, cR3, cR4 = create_region(windows[1], windows[2], size_ro, size_ri, windows_Size)
        if windows[0] not in dico_windows :
            #dico_windows[node]=[{"bk" : (bk1,bk2), "regions":[cR1, cR2, cR3, cR4],"R1":set(),"R2":set(),"R3":set(),"R4":set()}]
            dico_windows[windows[0]]=[{"bk" : (windows[1], windows[2]), "regions":[cR1, cR2, cR3, cR4],"R1":{},"R2":{},"R3":{},"R4":{}}]
        else :
            #dico_windows[node].append({"bk" : (bk1,bk2), "regions":[cR1, cR2, cR3, cR4], "R1":set(),"R2":set(),"R3":set(),"R4":set()})
            dico_windows[windows[0]].append({"bk" : (windows[1], windows[2]), "regions":[cR1, cR2, cR3, cR4], "R1":{},"R2":{},"R3":{},"R4":{}})
    '''

    # 3. List bx for each region
    n = 0
    x = 0
    with open(inputGAF, "r",encoding='UTF-8') as file : 
        for line in file :

            if line.startswith('@'):
                continue

            n += 1
            readID, readLen, __, __, __, path, __, pos_start, pos_end, __, alnLen, mapq, __, barcodeID, *__ = line.split("\t")
            

            # #Filters to keep only the valid alignments.
            if path == "*":
                x +=1             #remove unmapped reads
                continue
            # # cov = int(alnLen) / int(readLen)
            # # if cov < 0.9:
            # #     continue
            if int(mapq) < 20:
                continue

            #barcodeID = readID.split('/')[1]
            if path.startswith(">") :
                pos_node = re.split(r'[<>]',path)[1]
                pos_start = int(pos_start) + int(path.split(':')[1]) # Alignment coord : 0 to len(node) but node coord according to whole genome
                pos_end = int(pos_end) + int(path.split(':')[1])
            elif path.startswith("<") :
                pos_node = re.split(r'[<>]',path)[1]
                pos_start = int(path.split(":")[-1]) - int(pos_start) # Alignment coord : 0 to len(node) but node coord according to whole genome
                pos_end = int(path.split(":")[-1]) - int(pos_end)

            if pos_node not in all_nodes :
                continue

            #Filter to add bx in region list
            if pos_node in dico_windows :
                for windows in dico_windows[pos_node]:
                    r1, r2, r3, r4 = windows["regions"]
                    if r1[0]< pos_start and r4[1]> pos_end :
                        if r1[0]< pos_start and r1[1]> pos_end :
                            if barcodeID in windows["R1"].keys() :
                                windows["R1"][barcodeID] += 1
                            else :
                                windows["R1"][barcodeID] = 1

                        if r2[0]< pos_start and r2[1]> pos_end :
                            if barcodeID in windows["R2"].keys() :
                                windows["R2"][barcodeID] += 1
                            else :
                                windows["R2"][barcodeID] = 1

                        if r3[0]< pos_start and r3[1]> pos_end :
                            if barcodeID in windows["R3"].keys() :
                                windows["R3"][barcodeID] += 1
                            else :
                                windows["R3"][barcodeID] = 1

                        if r4[0]< pos_start and r4[1]> pos_end :
                            if barcodeID in windows["R4"].keys() :
                                windows["R4"][barcodeID] += 1
                            else :
                                windows["R4"][barcodeID] = 1

            # if pos_node in dico_windows :
            #     for windows in dico_windows[pos_node]:
            #         r1, r2, r3, r4 = windows["regions"]
            #         if r1[0]< pos_start and r1[1]> pos_end :
            #             windows["R1"].add(barcodeID)
            #         if r2[0]< pos_start and r2[1]> pos_end :
            #             windows["R2"].add(barcodeID)
            #         if r3[0]< pos_start and r3[1]> pos_end :
            #             windows["R3"].add(barcodeID)
            #         if r4[0]< pos_start and r4[1]> pos_end :
            #             windows["R4"].add(barcodeID)
    print("#Parse GAF file : Done", file = sys.stderr)


    # 4. Compare bx lists
    list_errors_rates = []
    count_GF = []
    filt_10bc =0
    w = windows_nb
    good_signal = 0
    false_signal = 0


    bc_good = 0
    bc_false = 0
    list_errors_rates_bc = []

    nb_bc_count = []
    nb_aln_count = []

    #print(dico_windows)
    for list_windows in dico_windows.values() :
        for windows in list_windows : 
            #print(windows)
            r1 = windows["R1"]
            r2 = windows["R2"]
            r3 = windows["R3"]
            r4 = windows["R4"]
            #print(len(r1),len(r2),len(r3),len(r4))
            if all(sum(windows[r].values()) > 2 for r in ["R1","R2","R3","R4"]):
                for barcode in windows["R1"].keys():
                    if barcode in windows["R2"] and barcode not in windows["R3"] and barcode not in windows["R4"] :
                        good_signal += windows["R1"][barcode] + windows["R2"][barcode]
                        bc_good +=1
                    if barcode in windows["R3"] and barcode not in windows["R2"] and barcode not in windows["R4"] :
                        false_signal += windows["R1"][barcode] + windows["R3"][barcode]
                        bc_false +=1
                for barcode in windows["R4"].keys():
                    if barcode in windows["R3"] and barcode not in windows["R1"] and barcode not in windows["R2"] :
                        good_signal += windows["R3"][barcode] + windows["R4"][barcode]
                        bc_good += 1
                    if barcode in windows["R2"] and barcode not in windows["R1"] and barcode not in windows["R3"] :
                        false_signal += windows["R4"][barcode] + windows["R2"][barcode]
                        bc_false += 1
                count_GF.append((good_signal,false_signal))

                # Calculate error rate
                N = good_signal+false_signal
                n_bc = bc_good +bc_false
                error_rate = false_signal / N if N >= 10 else None
                error_rate_bc = bc_false / n_bc if n_bc >= 10 else None
                #print(error_rate, N)
                if error_rate is not None:
                    if error_rate < 0.5 :
                        list_errors_rates.append(error_rate)
                        nb_aln_count.append((sum(r1.values()), sum(r2.values()), sum(r3.values()), sum(r4.values())))                      
                    else :
                        #print("Error rate to lower")
                        w = w - 1
                else :
                    #print("Moins de 10 reads total")
                    filt_10bc +=1
                    w = w - 1

                if error_rate_bc is not None:
                    if error_rate_bc < 0.5 :
                        list_errors_rates_bc.append(error_rate_bc)
                        nb_bc_count.append((len(r1), len(r2), len(r3), len(r4)))
            else :
                #print("Moins de deux reads par région")
                #print(len(r1),len(r2),len(r3),len(r4))
                w = w -1 

    # list_errors_rates = []
    # count_GF = []
    # filt_10bc =0
    # w = windows_nb
    # for list_windows in dico_windows.values() :
    #     for windows in list_windows : 
    #         r1 = windows["R1"]
    #         r2 = windows["R2"]
    #         r3 = windows["R3"]
    #         r4 = windows["R4"]
    #         if len(r1) > 2 and len(r2) > 2 and len(r3) > 2 and len(r4) > 2 :
    #             print(len(r1), len(r2), len(r3), len(r4))
    #             good_signal = len((r1 & r2) - r3 - r4) + len((r3 & r4) - r1 - r2)
    #             false_signal = len((r1 & r3) - r2 - r4) + len((r2 & r4) - r1 - r3)
    #             count_GF.append((good_signal,false_signal))

            #     # Calculate error rate
            #     N = good_signal+false_signal
            #     error_rate = false_signal / N if N >= 10 else None
            #     #print(error_rate, N)
            #     if error_rate is not None:
            #         if error_rate < 0.8 :
            #             list_errors_rates.append(error_rate)
            #         else :
            #             w = w - 1
            #     else :
            #         filt_10bc +=1
            #         w = w - 1

            # else :
            #     w = w -1 

    print(f"B\tWindows with less 10bc:",filt_10bc, file = sys.stderr)
    print(f"W\tTotal windows use to estimated:", w, file = sys.stderr)
    #print(f"BC\t{nb_bc_count}")
    #print(f"ALN\t{nb_aln_count}")

    #print(f"L\t{list_errors_rates}")
    if len(list_errors_rates) != 0 :
        mean_error = statistics.mean(list_errors_rates)
        median_error = statistics.median(list_errors_rates)
        std_error = statistics.stdev(list_errors_rates)

        print(f"E\tMean:{mean_error}\tStd:{std_error}\tMedian:{median_error}", file = sys.stderr)

        #5. Calculate the likelihood diff with the error estimate
        list_diff = []
        if mean_error != 0 :
            for (good_signal,false_signal) in count_GF:
                diff = likelihood_diffBests(good_signal,false_signal, mean_error)
                list_diff.append(diff)
            mean_diff = statistics.mean(list_diff)
        else :
            mean_diff = 'NA'

        
        #print(f"F\t{list_diff}")
        #print(f"D\tMean Diff Likelihood: {mean_diff}")


        #6. Create plot
        plot_creation(list_errors_rates, plot,windows_Size,w,size_ro, size_ri, mean_error,std_error,mean_diff,median_error)
    else :
        print("No windows with error rate above of 0.8 and/or a minimum number of 2 barcodes per regions and/or minimum 10 informatifs barcodes", file = sys.stderr)

    #print(f"BCL\t{list_errors_rates}")
    if len(list_errors_rates_bc) != 0 :
        mean_error_bc = statistics.mean(list_errors_rates_bc)
        median_error_bc = statistics.median(list_errors_rates_bc)
        std_error_bc = statistics.stdev(list_errors_rates_bc)

        #print(f"BCM\tMean:{mean_error_bc}\tStd:{std_error_bc}\tMedian:{median_error_bc}")


    #print("#All : Done")
    #print("nb reads tot", n, "nb unmap", x)
