# SVJedi-Tag 
[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html) 

SVJedi-Tag is a tool for genotyping inversions using linked-reads data. It is based on the creation of a graph from a reference genome and a VCF file containing the inversions to be genotyped. VG giraffe is then used to map reads onto the graph in order to analyse barcode signals specific to linked-reads and predict genotypes. 


---

## Installation

### Packages dependencies
- VG toolkits (version 1.54.0)
```
#VG installation command in a conda environment 
conda install bioconda::vg
```
- GFAGraphs (version 0.3.4) (Warning : GFAGraphs require python version 3.10 or higher and the packages networkx, tharos-pytools, matplotlib, mycolorpy)
```
pip install gfagraphs

# If necessary
pip install networkx
pip install tharos-pytools
pip install matplotlib
pip install mycolorpy
```

### Download 
```
git clone https://github.com/MTemperville/SVJedi-Tag.git
```

### Command to test the correct installation : 
```
python svjedi-tag.py -v testAuto/data/inversions_file.vcf -r testAuto/data/e_coli.fna -q testAuto/data/linked-reads.fastq -p testAuto/output_test/tag_test -s 10000 -t 8 
```

### Automatic test
```
bash test.sh 
```


---

## Files required 

* Linked-reads file (.fastq /.fq /.fastq.gz /.fa.gz)
* Variant calling file (.vcf)
* Reference genome file (.fasta/.fa)

#### File specificities : 
##### Linked-reads file requierement

The read header must contain the barcode and be in the following format: 
>headerBX:Z:XXXX

Please note that the header file must not contain tabulations.The *format_fastq.py* script is used to format the file so that the header corresponds to the expected format.
```
python format_fastq.py -q <linked-reads file> -o <linked-reads file pre-processing>
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

## Run tools

```
python svjedi-glr.py -v <vcf_file.vcf> -r <reference_genome.fa> -q <linked-reads_file.fastq> -p <output_path/prefix> -s <region_size> -t <thread>
```

**Warning:** The linked-reads file must be pre-processing to link BX tag to the header.

**Region_size parameter:**
*Genotyping is based on an analysis of barcodes in upstream and downstream regions and on the structure variant. The **region size** must be fixed. 
We recommend setting a region size greater than the size of the large DNA molecule used in the linked-reads data production protocol.*

## Contact 
SVJedi-graph is a Genscale tool developed by Mélody Temperville, Anne Guichard and Claire Lemaitre. For any bug report or feedback, please use the Github Issues form.
