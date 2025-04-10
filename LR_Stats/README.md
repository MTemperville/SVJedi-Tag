# LR_Stats

LR_Stats is a tool for analysing linked-read characteristics (number of barcodes, number of molecules, size of molecules, number of reads per molecule, etc.). LR_Stats uses the alignment positions of the linked-reads on a reference genome to reconstruct the original long DNA molecules. 

Input file: alignment_file_sorted.bam
Output file: statistique_tab.csv ; distribution_plot.png

Command : 
```
python LR_Stats.py -b [bam_file].bam 
```

You can have the help with :

```
python LR_Stats.py -h
```

```
usage: LR_Stats.py [-h] -b <sort_bam_file> [-s <molecule_max_size>] [-G <graphe_output_path/name_file>]
                   [-o <output_table_path/name_file>] [-g <genome_size>] [-r <read_size>]

options:
  -h, --help            show this help message and exit
  -b <sort_bam_file>, --bam <sort_bam_file>
                        Bam file out of a mapping linked-reads/reference and with BX tag
  -s <molecule_max_size>, --molecule_max_size <molecule_max_size>
                        Maximum size between two reads to share a same barcode and from the same molecule
                        [default=100000]
  -G <graphe_output_path/name_file>, --graph_output <graphe_output_path/name_file>
                        Path/Name_file for the output graph.png [default=ACLR_graph.png]
  -o <output_table_path/name_file>, --output_table <output_table_path/name_file>
                        Path/Name_file for the output table.csv [default=ACLR_table.csv]
  -g <genome_size>, --genome_size <genome_size>
                        Genome size required to calculate depth
  -r <read_size>, --read_size <read_size>
                        Read size required to calculate depth
(LREZ) [mtemperville@cl1n020:bam_file] $ 

```
