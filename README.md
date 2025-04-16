# SVJedi-Tag 
[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html) 

SVJedi-Tag is a tool for genotyping inversions using linked-read data. It is based on the analysis of the distribution of barcode signals on either sides of inversion breakpoints, with inversions being represented in a variation graph. The variation graph is built from a reference genome and a VCF file containing the inversions to be genotyped. [VG giraffe](https://github.com/vgteam/vg/wiki/Mapping-short-reads-with-Giraffe) is then used to map the sample reads onto the graph. Then, for each inversion, SVJedi-Tag analyzes the barcode signals of reads aligned on each side of the inversion breakpoints to estimate its allelic ratio and predict its genotype. 


---

## Installation

### Packages dependencies
- VG toolkits (version >1.54.0)
  
```
#VG installation command in a conda environment 
conda install bioconda::vg
```
- GFAGraphs (version 0.3.4) (Warning : GFAGraphs require python version 3.10 or higher and the packages networkx, tharos-pytools, matplotlib, mycolorpy)
```
pip install gfagraphs

# If necessary
pip install networkx
pip install matplotlib
pip install biopython
pip install tharos-pytools
```

### Download 
```
git clone https://github.com/MTemperville/SVJedi-Tag.git
```

### Command to test the correct installation : 
```
mkdir -p testAuto/output_test
python svjedi-tag.py -v testAuto/data/inversions_file.vcf -r testAuto/data/e_coli.fna -q testAuto/data/linked-reads.fastq -p testAuto/output_test/tag_test -s 10000 -t 8 
```

If the installation is successful, several files have been created in the directory `testAuto/output_test/`:
- the output file with predicted genotypes: `tag_test_genotype.vcf`
- intermediate files: `tag_test_analysis.txt`, `tag_test_chromDict.pickle`, `tag_test.dist`, `tag_test.gfa`, `tag_test.giraffe.gbz`, `tag_test.min`, `tag_test_vgGiraffe.gaf`.


### Automatic test
```
bash test.sh 
```
If the installation is successful, you will get the following output in the terminal :

```
### Create variant graph ###
### Index graph ###
[vg autoindex] Executing command: vg autoindex --workflow giraffe -g testAuto/output_test/tag_test.gfa -p testAuto/output_test/tag_test
[IndexRegistry]: Checking for haplotype lines in GFA.
[IndexRegistry]: Constructing VG graph from GFA input.
[IndexRegistry]: Constructing XG graph from VG graph.
[IndexRegistry]: Constructing a greedy path cover GBWT
[IndexRegistry]: Constructing GBZ using NamedNodeBackTranslation.
[IndexRegistry]: Constructing distance index for Giraffe.
[IndexRegistry]: Constructing minimizer index.
### Map linked-reads on graph ###
### Analyze barcode signal & Genotype ###
Done. Output genotypes in file testAuto/output_test/tag_test_genotype.vcf
### Testing SVJedi-Tag ###
Test-graph:PASS
Test-VCF:PASS
```

---

## Input Files 

* Linked-reads file (.fastq /.fq /.fastq.gz /.fa.gz)
* Variant calling file (.vcf)
* Reference genome file (.fasta/.fa)

#### File specificities : 
##### Linked-read file requirement

The sequence of each read should not contain the barcode sequence. The barcode information must be present in the header of each read in the following format (no space before the barcode): 
>headerBX:Z:XXXX

Please note that the full headers must not contain any tabulation. 

The `format_fastq.py` script can be used to format the file so that the headers correspond to the expected format.
```
python format_fastq.py -q <initial linked-read file> -o <formated linked-reads file>
```

##### VCF 
* The file must be a tabulated file with a standart vcf format : 
>#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO

*Example:* 
```
##fileformat=VCFv4.3
##fileDate=2025-03-05
##ALT=<ID=INV,Description="Inversion">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position on reference genome">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV">
##INFO=<ID=SVSIZE,Number=1,Type=Integer,Description="Size in nt of the SV">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
LG1	1173997	INV1	N	<INV>	.	.	END=2173996;SVTYPE=INV;SVSIZE=1000000
LG1	2287235	INV2	N	<INV>	.	.	END=3287234;SVTYPE=INV;SVSIZE=1000000
LG1	4393271	INV3	N	<INV>	.	.	END=5393270;SVTYPE=INV;SVSIZE=1000000
```

## Command line and parameters

```
python svjedi-glr.py -v <vcf_file.vcf> -r <reference_genome.fa> -q <linked-reads_file.fastq> -p <output_path/prefix> -s <region_size> -t <thread>
```

**Warning:** The linked-read file must be pre-processed to link BX tag to the header.

### Paramters
```
svjedi-tag.py [-h] -v <inputVCF> -r <referenceGenome> -q <queryReads> -p <outFilesPrefix> [-t <threadNumber>] [-s <regionSize (default 10000)>] [-a <alignmentGAFFile>] [-g <Graphe File GFA>]

options:
  -h, --help        Show this help message and exit
  -v, --vcf         Input VCF with structural variants to be genotyped
  -r, --ref         Reference genome 
  -q, --reads       Linked-reads file
  -p, --prefix      Output prefix 
  -t, --threads     Number of thread used
  -s, --regionSize  Region genotyping lenght (default 10000)
  -a, --gaf         GAF file
  -g, --gfa         GFA file

```

**Input VCF:**
The file must contain the inversions to be genotyped. It may also contain other types of structural variants that could have an impact on the genotyping of the inversions because of their position close to the breakpoints of the inversions.

**Region_size parameter:**
Genotyping is based on the analysis of barcodes in regions from either side of both inversion breakpoints. The **size** of these regions (in bp) is fixed by the `-s`(`--regionSize`) parameter.   
We recommend setting a region size similar to the average size of the large DNA molecules obtained with the linked-read data production protocol (default value = 10,000 bp).

## Output files

#### Output VCF 
SVJedi-Tag's output is a vcf file containing all the input inversions associated to their genotype. It contains a *SAMPLE* column containing the predicted genotypes and genotyping information.

- GT : Genotype, assuming a diploid individual: either `0/0`, `0/1`, or `1/1`. A `./.` value indicates a missing value due to an insufficient number of informative barcodes in the breakpoint regions.
- DP : Total number of distinct barcodes in the breakpoints regions
- AD : Number of barcodes supporting the absence of the inversion, Number of barcodes supporting the presence of the inversion, Number of non-informative barcodes.
- AF : Alternative allele ratio

```
Example : 
GT:DP:AD:AF	1/1:159:0,99,60:1.0 
```
#### GFA file
The GFA file contains the graph (in GFA format) constructed from the reference genome and the intput VCF file. It can be used to restart SVJedi-Tag by skipping the graph creation step.
#### GAF file 
The GFA file is the VG Giraffe output file (in GAF format), it contains the read alignments to the graph. It can be used to restart SVJedi-Tag by skipping the alignment step.

## Contact 
SVJedi-graph is a [Genscale](https://team.inria.fr/genscale/) tool developed by Mélody Temperville, Anne Guichard and Claire Lemaitre. For any bug report or feedback, please use the Github Issues form.
