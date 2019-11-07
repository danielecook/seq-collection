#!/bin/bash
test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

PARENT_DIR=`git rev-parse --show-toplevel`
export PATH="${PATH}:${PARENT_DIR}"

#======#
# json #
#======#

set -o nounset
run test_json sc json "${PARENT_DIR}/tests/data/test.vcf.gz" X:17276844-17276844
assert_equal \"X\" "$(cat "$STDOUT_FILE" | jq '.CHROM')"
assert_equal 17276844 "$(cat "$STDOUT_FILE" | jq '.POS')"
assert_equal "\"PASS\"" "$(cat "$STDOUT_FILE" | jq '.FILTER[0]')"
assert_equal 999 "$(cat "$STDOUT_FILE" | jq '.QUAL')"
assert_equal \"T\" "$(cat "$STDOUT_FILE" | jq '.REF')"

run test_json_pretty sc json --pretty "${PARENT_DIR}/tests/data/test.vcf.gz" X:17276844-17276844
assert_equal 13 $(cat $STDOUT_FILE | wc -l)

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

#==========#
# fq-dedup #
#==========#

run fq_dedup sc fq-dedup tests/fastq/dup.fq
assert_equal 4 $(cat $STDOUT_FILE | grep '@' | wc -l)
assert_equal 4 $(cat $STDOUT_FILE | grep '@' | wc -l)

run fq_dedup_gz sc fq-dedup tests/fastq/dup.fq.gz
assert_equal 4 $(cat $STDOUT_FILE | grep '@' | wc -l)
assert_equal 4 $(cat $STDOUT_FILE | grep '@' | wc -l)

#=========#
# fq-meta #
#=========#
curl https://raw.githubusercontent.com/10XGenomics/supernova/master/tenkit/lib/python/tenkit/illumina_instrument.py | \
        grep -v 'import martian' | \
        sed 's/martian.exit/exit/g' | \
        sed 's/f.readline()/f.readline().decode("utf-8")/g' | \
        sed 's/exit/print/g' > "${PARENT_DIR}/illumina_instrument.py"

# Fix python script
echo "Fixing ${PARENT_DIR}/illumina_instrument.py"
2to3 --write --nobackups "${PARENT_DIR}/illumina_instrument.py"

function test_fq() {
    gzip -c $1 > test.fq.gz
    python -c "import illumina_instrument as ia; print(ia.sequencer_detection_message(['test.fq.gz']))"
}

# Genome Analyzer II
run fq_meta_genome_analyzer1 sc fq-meta tests/fastq/illumina_1.fq
assert_equal "$(cat $STDOUT_FILE | cut -f 2)" \
             "$(test_fq tests/fastq/illumina_1.fq | egrep -o "Genome Analyzer IIx" | uniq | sed 's/ //g')"
assert_equal "$(cat $STDOUT_FILE | cut -f 3 | cut -f 1 -d ':')" \
             "$(test_fq tests/fastq/illumina_1.fq | egrep -o "likely" | uniq)"

run fq_meta_genome_analyzer2 sc fq-meta tests/fastq/illumina_2.fq
assert_equal "$(cat $STDOUT_FILE | cut -f 2)" \
             "$(test_fq tests/fastq/illumina_2.fq | egrep -o "Genome Analyzer IIx" | uniq | sed 's/ //g')"
assert_equal "$(cat $STDOUT_FILE | cut -f 3 | cut -f 1 -d ':')" \
             "$(test_fq tests/fastq/illumina_2.fq | egrep -o "likely" | uniq)"

# Unresolvable; This FASTQ is from wikipedia...
run fq_meta3 sc fq-meta tests/fastq/illumina_3.fq
assert_equal "$(cat $STDOUT_FILE | cut -f 2)" ""
assert_equal "$(cat $STDOUT_FILE | cut -f 3 | cut -f 1 -d ':')" ""
assert_equal "$(test_fq tests/fastq/illumina_3.fq | egrep -o 'could not detect')" "could not detect"

run fq_meta4 sc fq-meta tests/fastq/illumina_4.fq
assert_equal "$(cat $STDOUT_FILE | cut -f 2)" ""
assert_equal "$(cat $STDOUT_FILE | cut -f 3 | cut -f 1 -d ':')" ""
assert_equal "$(test_fq tests/fastq/illumina_4.fq | egrep -o 'could not detect')" "could not detect"

run fq_meta_hiseq_2000_2500 sc fq-meta tests/fastq/illumina_2000_2500.fq
assert_equal "$(cat $STDOUT_FILE | cut -f 2)" "HiSeq2000/2500"
assert_equal "$(cat $STDOUT_FILE | cut -f 3)" "high:machine+flowcell"
assert_equal "$(cat $STDOUT_FILE | cut -f 2)" \
             "$(test_fq tests/fastq/illumina_2000_2500.fq | egrep -o "HiSeq2000/2500" | uniq)"
assert_equal "$(cat $STDOUT_FILE | cut -f 3 | cut -f 1 -d ':')" \
             "$(test_fq tests/fastq/illumina_2000_2500.fq | egrep -o 'high' )"

# HiSeq 4000
run fq_meta_hiseq_4000 sc fq-meta tests/fastq/illumina_3000_4000.fq
assert_equal "$(cat $STDOUT_FILE | cut -f 2)" "HiSeq3000/4000"
assert_equal "$(cat $STDOUT_FILE | cut -f 3)" "high:machine+flowcell"
assert_equal "$(cat $STDOUT_FILE | cut -f 2)" \
             "$(test_fq tests/fastq/illumina_3000_4000.fq | egrep -o "HiSeq3000/4000" | uniq)"
assert_equal "$(cat $STDOUT_FILE | cut -f 3 | cut -f 1 -d ':')" \
             "$(test_fq tests/fastq/illumina_3000_4000.fq | egrep -o 'high')"

# HiSeq X
run fq_meta_hiseq_x sc fq-meta tests/fastq/illumina_hiseq_x.fq
assert_equal "$(cat $STDOUT_FILE | cut -f 2)" "HiSeqX"
assert_equal "$(cat $STDOUT_FILE | cut -f 3)" "high:machine+flowcell"
assert_equal "$(cat $STDOUT_FILE | cut -f 3 | cut -f 1 -d ':')" \
             "$(test_fq tests/fastq/illumina_hiseq_x.fq | egrep -o 'high')"

# Novaseq
run fq_meta_novaseq sc fq-meta tests/fastq/novaseq.fq
assert_equal "$(cat $STDOUT_FILE | cut -f 2)" "NovaSeq"
assert_equal "$(cat $STDOUT_FILE | cut -f 3)" "high:machine+flowcell"

#=============#
# Insert size #
#=============#
run insert_size sc insert-size tests/data/test.bam
assert_equal "$(cat $STDOUT_FILE | cut -f 1 | tail -n 1)" "179"
assert_equal "$(cat $STDOUT_FILE | cut -f 2 | tail -n 1)" "176.5"
assert_equal "$(cat $STDOUT_FILE | cut -f 4 | tail -n 1)" "38"
assert_equal "$(cat $STDOUT_FILE | cut -f 5 | tail -n 1)" "358"
assert_equal "$(cat $STDOUT_FILE | cut -f 6 | tail -n 1)" "359"