#!/bin/bash
test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtests/master/ssshtest

. ssshtest

set -o nounset

assert_equal 18 12

# slice
run slice_inf csv slice 5: tests/data/*.tsv
assert_exit_code 0
assert_no_stderr
assert_equal "18" $(cat $STDOUT_FILE | wc -l)
assert_in_stdout "Dodge Challenger"

run slice_low csv slice :3 tests/data/*.tsv
assert_exit_code 0
assert_no_stderr
assert_equal "12" $(cat $STDOUT_FILE | wc -l)
assert_in_stdout "Fiat 128"

run slice csv slice 1:3 tests/data/*.tsv
assert_exit_code 0
assert_no_stderr
assert_equal "9" $(cat $STDOUT_FILE | wc -l)

run slice csv slice huh:3 tests/data/*.tsv
assert_exit_code 1
assert_stderr
assert_in_stderr "Malformed range"

run slice_add_col csv slice -a 0:3 tests/data/*.tsv
assert_exit_code 0
assert_in_stdout "filename"
assert_equal "6" $(cat $STDOUT_FILE | cut -f 1 | uniq | wc -l)