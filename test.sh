# Run SVJedi-GLR

LR="testAuto/data/INV10-20kb_div400_ok.fastq"
VCF="testAuto/data/INV10-20kb_HMZAltRef+invclose.vcf"
REF="testAuto/data/e_coli.fna"
OUTPUT="testAuto/output_test/inv10-20_HMZ_DFS_div400"

date
python svjedi-glr.py -v $VCF -r $REF -q $LR -p $OUTPUT -s 10000 -t 8


# Automatic test

if diff "testAuto/output/inv10-20_HMZ_DFS_div400.gfa" "testAuto/output_test/inv10-20_HMZ_DFS_div400.gfa" > /dev/null; then
    echo "Same graphs"
else
    echo "Warning : Graphs are differents. Please check in construct_graph script."
fi

if diff "testAuto/output/inv10-20_HMZ_DFS_div400_genotype.vcf" "testAuto/output_test/inv10-20_HMZ_DFS_div400_genotype.vcf" > /dev/null; then
    echo "Same  VCF"
else
    echo "Warning : VCF are differents. Please check in genotype script"
fi

date