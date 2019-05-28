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
