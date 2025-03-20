# Run SVJedi-GLR

LR="testAuto/data/linked-reads.fastq"
VCF="testAuto/data/inversions_file.vcf"
REF="testAuto/data/e_coli.fna"
OUTPUT="testAuto/output_test/tag_test"

date
python svjedi-tag.py -v $VCF -r $REF -q $LR -p $OUTPUT -s 10000 -t 8


# Automatic test

if diff "testAuto/output/tag_test.gfa" "testAuto/output_test/tag_test.gfa" > /dev/null; then
    echo "Same graphs"
else
    echo "Warning : Graphs are differents. Please check in construct_graph script."
fi

if diff "testAuto/output/tag_test_genotype.vcf" "testAuto/output_test/tag_test_genotype.vcf" > /dev/null; then
    echo "Same  VCF"
else
    echo "Warning : VCF are differents. Please check in genotype script"
fi

date