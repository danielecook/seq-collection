#!/bin/bash
test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

PARENT_DIR="`git rev-parse --show-toplevel`"
export PATH="${PATH}:${PARENT_DIR}"

set -o nounset
run test_json sc json "${PARENT_DIR}/tests/data/test.vcf.gz" X:17276844-17276844
assert_equal \"X\" "$(cat "$STDOUT_FILE" | jq '.CHROM')"
assert_equal 17276844 "$(cat "$STDOUT_FILE" | jq '.POS')"
assert_equal "\"PASS\"" "$(cat "$STDOUT_FILE" | jq '.FILTER[0]')"
assert_equal 999 "$(cat "$STDOUT_FILE" | jq '.QUAL')"
assert_equal \"T\" "$(cat "$STDOUT_FILE" | jq '.REF')"

run test_json_pretty sc json --pretty "${PARENT_DIR}/tests/data/test.vcf.gz" X:17276844-17276844
assert_equal 13 "$(cat $STDOUT_FILE | wc -l)"

# INFO
run single_info_item sc json --info="DP" "${PARENT_DIR}/tests/data/test.vcf.gz" X:17276844-17276844
assert_equal 9836 "$(cat $STDOUT_FILE | jq '.INFO.DP')"

run multi_info_item sc json --info="DP,MQ,DP4,HOB,INDEL" "${PARENT_DIR}/tests/data/test.vcf.gz" X:17276844-17276844
assert_equal 9836 "$(cat $STDOUT_FILE | jq '.INFO.DP')"
assert_equal 60 "$(cat $STDOUT_FILE | jq '.INFO.MQ')"
assert_equal 92 "$(cat $STDOUT_FILE | jq '.INFO.DP4[3]')"
assert_equal 0.5 "$(cat $STDOUT_FILE | jq '.INFO.HOB')"

# FORMAT
run format_dp sc json --format="DP" "${PARENT_DIR}/tests/data/test.vcf.gz" I:41947-41947 | jq '.FORMAT.DP|add'
assert_equal 2094 "$(cat $STDOUT_FILE | jq '.FORMAT.DP|add')"


# fq-meta
curl https://raw.githubusercontent.com/10XGenomics/supernova/master/tenkit/lib/python/tenkit/illumina_instrument.py | \
        grep -v 'import martian' | \
        sed 's/martian.exit/exit/g' > illumina_instrument.py

function test_fq() {
    gzip -c $1 > test.fq.gz
    python -c "import illumina_instrument as ia; print(ia.sequencer_detection_message(['test.fq.gz']))"
}

# Genome Analyzer II
run fq_meta sc fq-meta tests/fastq/illumina_1.fq
assert_equal "$(cat $STDOUT_FILE | cut -f 7)" \
             "$(test_fq tests/fastq/illumina_1.fq | egrep -o "Genome Analyzer IIx" | uniq)"
assert_equal "$(cat $STDOUT_FILE | cut -f 8 | cut -f 1 -d ':')" \
             "$(test_fq tests/fastq/illumina_1.fq | egrep -o "likely" | uniq)"

run fq_meta2 sc fq-meta tests/fastq/illumina_2.fq
assert_equal "$(cat $STDOUT_FILE | cut -f 7)" \
             "$(test_fq tests/fastq/illumina_2.fq | egrep -o "Genome Analyzer IIx" | uniq)"
assert_equal "$(cat $STDOUT_FILE | cut -f 8 | cut -f 1 -d ':')" \
             "$(test_fq tests/fastq/illumina_2.fq | egrep -o "likely" | uniq)"

# Unresolvable
run sc fq_meta tests/fastq/illumina_3.fq
assert_equal "$(cat $STDOUT_FILE | cut -f 7)" ""
assert_equal "$(cat $STDOUT_FILE | cut -f 8 | cut -f 1 -d ':')" ""