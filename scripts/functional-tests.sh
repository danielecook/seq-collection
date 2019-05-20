#!/bin/bash
test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtests/master/ssshtest

. ssshtest

set -o nounset

# slice
run vcf2json sk
assert_exit_code 0
assert_no_stderr
assert_equal "18" $(cat $STDOUT_FILE | wc -l)