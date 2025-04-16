# Run SVJedi-GLR

LR="testAuto/data/linked-reads.fastq"
VCF="testAuto/data/inversions_file.vcf"
REF="testAuto/data/e_coli.fna"
OUTPUT="testAuto/output_test/tag_test"

mkdir -p testAuto/output_test

rm -f testAuto/output_test/*

date
python svjedi-tag.py -v $VCF -r $REF -q $LR -p $OUTPUT -s 10000 -t 8


# Automatic test

echo "### Testing SVJedi-Tag ###"

if diff "testAuto/output/tag_test.gfa" "testAuto/output_test/tag_test.gfa" > /dev/null; then
    echo "Test-graph:PASS"
else
    echo "Test-graph:FAILED - Graph files are different. Please check in construct_graph script."
fi

if diff "testAuto/output/tag_test_genotype.vcf" "testAuto/output_test/tag_test_genotype.vcf" > /dev/null; then
    echo "Test-VCF:PASS"
else
    echo "Test-VCF:FAILED  - VCF files are different. Please check in genotype script"
fi

date
